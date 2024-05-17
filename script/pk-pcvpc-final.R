# --- 
# title: pcVPC with the vpc package
# ---
# 
# # Scope
# 
# This document illustrates the mechanics of using the vpc package to 
# create a prediction-corrected vpc. 
# 
# We will create one vpc plot: pred-corrected vpc stratified by STUDY.
# 

# # Required packages
library(tidyverse)
library(glue)
library(yspec)
library(mrgsolve)
library(mrggsave)
library(vpc)
library(bbr)
library(bbr.bayes)
library(here)
library(posterior)

options(mrggsave.dir = here("deliv/figure"), mrg.script = "pk-pcvpc-final.R")
options(mrgsolve.project = here("script/model"))

mrg_vpc_theme = new_vpc_theme(list(
  sim_pi_fill = "steelblue3", sim_pi_alpha = 0.5,
  sim_median_fill = "grey60", sim_median_alpha = 0.5
))

spec <- ys_load(here("data/derived/pk.yml"))
lab <- ys_get_short_unit(spec, parens = TRUE, title_case = TRUE)

# 
# # Model
# 
runno <- 1100
MODEL_DIR <- here("model/pk")

# Read in complete dataset; don't forget to filter out other records that may
# have been ignored in the run
mod_bbr <- read_model(file.path(MODEL_DIR, runno))

# We also need the medians of the posterior parameter estimates. By default,
# the mrgsolve model is set up to read in only the estimates from the first
# chain.

# Extract posterior draws
draws <- read_fit_model(mod_bbr)
# Calculate medians
params <- draws %>% 
  subset_draws(variable = c("THETA", "OMEGA", "SIGMA")) %>% 
  summarize_draws(median) %>% 
  mutate(variable = str_remove_all(variable, "[:punct:]")) %>% 
  pivot_wider(names_from = variable, values_from = median)
theta <- params %>% select(starts_with("THETA"))
# This may need some tweaking for different block structures for OMEGA or SIGMA
omega <- params %>% 
  select(starts_with("OMEGA")) %>% 
  bmat()
sigma <- params %>% 
  select(starts_with("SIGMA")) %>% 
  bmat()

# ## Load the mrgsolve model
# 
# This should reflect `../model/pk/1100.ctl`
mod0 <- mread(glue("{runno}.mod"))
mod <- mod0 %>% 
  param(theta) %>% 
  omat(omega) %>% 
  smat(sigma)

# Generate `EPRED` with `nm_join_bayes` to use for prediction corrections.
# Since we only need `EPRED` (and the original data) we can skip IPRED, EWRES,
# and NPDE.
# Use `.superset = TRUE` to recover the full dataset, including BLQ records.
set.seed(1)
progressr::with_progress({
  data0 <- mod_bbr %>% 
    nm_join_bayes(
      mod_mrgsolve = mod,
      .superset = TRUE,
      ipred = FALSE,
      ewres_npde = FALSE
    )
})
saveRDS(data0, file.path(MODEL_DIR, runno, glue("diag-sims-epred-{runno}.rds")))
data0 <- readRDS(file.path(MODEL_DIR, runno, glue("diag-sims-epred-{runno}.rds")))
data <- filter(data0, is.na(C))


# Simulate `EPRED` for the records that are BLQ, otherwise we'll use the 
# value coming from NONMEM
data %>% filter(EVID == 0) %>% pull(EPRED) %>% anyNA() # TRUE
out <- mrgsim(zero_re(mod), data, digits = 5)
data <- mutate(data, EPRED = ifelse(BLQ > 0, out$Y, EPRED))
data %>% filter(EVID == 0) %>% pull(EPRED) %>% anyNA() # FALSE


# # Set up the simulation
# 
# Create a function to simulate out one replicate
sim <- function(rep, data, model) {
  mrgsim(
    model, 
    data = data,
    carry_out = "EVID,STUDYN,LDOS,DOSE,EPRED",
    recover = "STUDY,C,USUBJID,ACTARM,RF,Renal,Hepatic", 
    Req = "Y", 
    output = "df"
  ) %>%  mutate(irep = rep)
}

# Simluate data
isim <- seq(200)

set.seed(86486)
sims <- lapply(
  isim, sim, 
  data = data, 
  mod = mod
) %>% bind_rows()

# Filter both the observed and simulated data
# For the observed data, we only want actual observations that weren't BLQ
# For the simulated data, we take simulated observations that were above LQ
fdata <- filter(data, EVID==0, BLQ == 0)
fsims <- filter(sims, EVID==0, Y >= 10)

# # Create the plot
# 
# Pass observed and simulated data into vpc function
p1 <- vpc(
  obs = fdata,
  sim = fsims,
  pred_corr = TRUE,
  stratify = "STUDY",
  obs_cols = list(pred = "EPRED"),
  sim_cols = list(dv = "Y", sim = "irep", pred = "EPRED"), 
  log_y = TRUE,
  pi = c(0.05, 0.95),
  ci = c(0.025, 0.975), 
  vpc_theme = mrg_vpc_theme
) 

p1 <- 
  p1 +  
  theme_bw() + 
  xlab(lab$TIME) + 
  ylab("Prediction-corrected concentration (ng/mL)")

p1

mrggsave(p1, stem = "pk-vpc-{runno}-pred-corr", height = 6)
