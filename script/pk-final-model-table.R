### Libraries ----------------------------
library(tidyverse)
library(bbr)
library(bbr.bayes)
library(pmtables)
library(here)
library(posterior)
library(magrittr)


### Directories ----------------------------
tabDir <- here("deliv", "table", "report")
if(!file.exists(tabDir)) dir.create(tabDir)

thisScript <- "pk-final-model-table.R"
options(
  mrg.script = thisScript, 
  pmtables.dir = tabDir
)


set.seed(5238974)


# Helper functions ----------------------------
source(here("script", "functions-table.R"))

# Read in model output ----------------------------
run <- 1100
MODEL_DIR <- here("model/pk")

# Read in model object
mod_bbr <- read_model(file.path(MODEL_DIR, run))
# Extract posterior draws
draws <- read_fit_model(mod_bbr)

draws_param <- draws %>% 
  subset_draws(variable = c("THETA", "OMEGA", "SIGMA"))

# Calculate shrinkage from post hoc ETAs, or from the .shk files if the .iph
# files do not exist
shk0 <- shrinkage(mod_bbr)
omegas <- variables(draws) %>% 
  str_subset("^OMEGA") %>% 
  # diagonals only
  str_subset("\\[(\\d+),\\1\\]")
shk <- tibble(name = omegas, shrinkage = shk0)

# parameter estimates and MCMC diagnostics ------------------------------------

ptable <- get_ptable(draws_param)
ptable

# get parameter names  ---------------------------- ----------------------------
param_key <- yaml_as_df(here("script", "pk-parameter-key.yaml")) %>% 
  rename(name = .row)

# Extract PK parameters and generate values to be displayed for report table ----------------------------
param_df <- ptable %>% 
  left_join(shk, by = "name") %>% 
  # remove "[" and "]" (e.g., THETA[1] -> THETA1)
  mutate(name = str_remove_all(name, "[[:punct:]]")) %>% 
  inner_join(param_key, by = "name") %>% 
  checkTransforms() %>%               # check for associated THETAs (e.g. for logit transformations)
  defineRows() %>%                    # define series of T/F variables
  getValueSE() %>%                    # define which value and se are required
  formatValues() %>%                  # back transform as needed, round using' sig' and combine columns where needed
  mutate(parameter_names = name) %>% 
  formatGreekNames() %>%              # format the labels to display greek symbols
  getPanelName() %>%                  # Define panel names based on parameter type
  mutate(
    rhat = sig(rhat),
    across(c(ess_bulk, ess_tail), round)
  ) %>% 
  # select columns of interest
  dplyr::select(type, abb, greek, desc, value, ci, shrinkage,
                rhat, ess_bulk, ess_tail)

param_df 


# Define footnotes ---------------------------- ----------------------------
footAbbrev <- "Abbreviations: CrI = credible interval"
footAbbrev_Om <- "Abbreviations: CrI = credible interval; 
                        Corr = Correlation coefficient;
                        CV = coefficient of variation"
footDerive2 <- "CV\\% of log-normal omegas = sqrt(exp(estimate) - 1) $\\cdot$ 100"
footDerive3 <- "CV\\% of sigma = sqrt(estimate) $\\cdot$ 100"
footLog <- "Parameters estimated in the log-domain were back-transformed for clarity"
footMCMC <- "Abbreviations: 
  ESS = effective sample size;
  $\\hat{R}$ = Gelman-Rubin diagnostic"




##  FIXED EFFECTS tables ----------------------------

# Parameter estimates
fixed <- param_df %>% 
  filter(str_detect(type, "Struct") | 
           str_detect(type, "effect")) %>%
  select(-shrinkage, -rhat, -ess_bulk, -ess_tail) %>% 
  st_new() %>% 
  st_panel("type") %>% 
  st_center(desc = col_ragged(5.5),
            abb = "l") %>% 
  st_blank("abb", "greek", "desc") %>% 
  st_rename("Estimate" = "value", 
            "95\\% CrI" = "ci") %>% 
  st_files(output = "pk-param-final-fixed.tex") %>% 
  st_notes(footLog, footAbbrev) %>% 
  st_noteconf(type = "minipage", width = 1) %>% 
  stable() %T>% 
  # st2report(ntex = 2) %T>%
  stable_save()

# MCMC diagnostics
fixed_mcmc <- param_df %>% 
  filter(str_detect(type, "Struct") | 
           str_detect(type, "effect")) %>%
  select(-shrinkage, -value, -ci) %>% 
  st_new() %>% 
  st_panel("type") %>% 
  st_center(desc = col_ragged(5.5),
            abb = "l") %>% 
  st_blank("abb", "greek", "desc") %>% 
  st_rename("$\\hat{R}$" = "rhat", 
            "ESS bulk" = "ess_bulk",
            "ESS tail" = "ess_tail") %>% 
  st_files(output = "pk-param-final-fixed-mcmc.tex") %>% 
  st_notes(footMCMC) %>% 
  st_noteconf(type = "minipage", width = 1) %>% 
  stable() %T>% 
  # st2report(ntex = 2) %T>%
  stable_save()





##  RANDOM EFFECTS table ----------------------------

# Parameter estimates
random <- param_df %>% 
  filter(str_detect(greek, "Omega") | 
           str_detect(type, "Resid")) %>%
  select(-desc, -rhat, -ess_bulk, -ess_tail) %>%
  st_new %>% 
  st_panel("type") %>% 
  st_center(abb = "l") %>% 
  st_blank("abb", "greek") %>% 
  st_rename("Estimate" = "value", 
            "95\\% CrI" = "ci",
            "Shrinkage (\\%)" = "shrinkage") %>% 
  st_files(output = "pk-param-final-random.tex") %>% 
  st_notes(footAbbrev_Om, footDerive2, footDerive3) %>% 
  st_noteconf(type = "minipage", width = 1) %>% 
  stable() %T>% 
  # st2report(ntex = 2) %T>%
  stable_save()

# MCMC diagnostics
random_mcmc <- param_df %>% 
  filter(str_detect(greek, "Omega") | 
           str_detect(type, "Resid")) %>%
  select(-desc, -shrinkage, -value, -ci) %>% 
  st_new %>% 
  st_panel("type") %>% 
  st_center(abb = "l") %>% 
  st_blank("abb", "greek") %>% 
  st_rename("$\\hat{R}$" = "rhat", 
            "ESS bulk" = "ess_bulk",
            "ESS tail" = "ess_tail") %>% 
  st_files(output = "pk-param-final-random-mcmc.tex") %>% 
  st_notes(footMCMC) %>% 
  st_noteconf(type = "minipage", width = 1) %>% 
  stable() %T>% 
  # st2report(ntex = 2) %T>%
  stable_save()


##  Save tables out to pdf preview ----------------------------
# Check they fit within the report margins
st2report(
  list(fixed, random, fixed_mcmc, random_mcmc),
  ntex = 2,
  output_dir = tabDir,
  stem = "preview-final-param-table",
  show_pdf = interactive()
)
