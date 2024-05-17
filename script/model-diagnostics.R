### Purpose -------------------------------
## Produce diagnostic plots using Rmd templates
## (located in the diagnostics-templates folder)
## By defining model specific characteristics here, we can simply render the 
## Rmd templates to produce the diagnostics relevant to a specific model
## 
## This particular script focuses on producing diagnostics from two templates
## (one for model diagnostics and one for MCMC diagnostics)
## but it may be necessary for you to have other templates (depending how 
## different the models are).
## This script is simply meant to provide an example of how you can run the same
## diagnostic templates for multiple similar models. In this case, the templates
## focus on creating diagnostics that will be included in the final report but
## please note that these are just examples for this report. The template and
## choice of diagnostics should be modified to reflect the project specific
## questions being addressed in your project.
##
## Support for parameterized reports can be found 
## https://bookdown.org/yihui/rmarkdown/parameterized-reports.html


### Libraries ----------------------------
library(tidyverse)
library(glue)
library(bbr)
library(bbr.bayes)
library(here)
library(posterior)
library(mrgsolve)

### Script ----------------------------
thisScript <- "model-diagnostics.R"

### Model directory -------------------
model_dir <- here("model/pk")

source(here("script/functions-diagnostics.R"))

# Parallel options --------------------------------------------------------
options(future.fork.enable = TRUE)
future::plan(future::multicore, workers = parallelly::availableCores() - 1)

# diagnostic-templates > template-model-bayes.Rmd
# diagnostic-templates > template-mcmc.Rmd

# The final results included in the report were created with
# template-model-bayes.Rmd and template-mcmc.Rmd.
# This script was written with the report specifications in mind and generates
# the pdf figures appropriate to this report. Please go to the
# template-model-bayes.Rmd to see the fields that can be edited directly in this
# R script.


# The idea is that you build up the model diagnostic code shown below
# as you work though you project. We recommend keeping this code in case you
# need to re-run the diagnostic plots for some reason - it saves you re-typing
# all the model-specific details to re-create the diagnostics for multiple
# models.
# After an initial look at the diagnostics, run the `run-sims()`            
# function to generate the appropriate report diagnostics.

### run 1000 - Base model (Bayes) ----------------------------
model_name <- 1000

# Render the diagnostics report
# First we run the diagnostics using only output from NONMEM tables. Since the
# simulation-based diagnostics using the full posterior (i.e. using results from
# all chains) requires an mrgsolve model and could be time-consuming, it is
# reasonable to use only NONMEM output until the model is confirmed as final.
# To achieve this, we leave the defaults `sims_output_path = ""` so that
# simulation results aren't used and `nm_join_bayes_quick` is used instead.
rmarkdown::render(
  here("script/diagnostic-templates/template-model-bayes.Rmd"),
  # Only the run number parameter and the flag to disable running mrggsave are
  # set here.
  # All other parameters are set to default values: see template-model-bayes.Rmd
  # for defaults and descriptions.
  params = list(
    run = model_name,
    run_mrggsave = FALSE
  ),
  # the html that is created will be saved to the model directory
  output_dir = file.path(model_dir, model_name),
  output_file = glue("model-diagnostics-{model_name}.html")
)
# View the output
utils::browseURL(file.path(model_dir, model_name,
                           glue("model-diagnostics-{model_name}.html")))

# MCMC diagnostics
# Next we run the separate MCMC diagnostics to assess the MCMC sampling
rmarkdown::render(
  here("script/diagnostic-templates/template-mcmc.Rmd"),
  params = list(run = model_name),
  output_dir = file.path(model_dir, model_name),
  output_file = glue("mcmc-diagnostics-{model_name}.html")
)
# View the output
utils::browseURL(file.path(model_dir, model_name,
                           glue("mcmc-diagnostics-{model_name}.html")))

### run 1001 - Base model (NUTS) ----------------------------
model_name <- 1001

# After an initial look at the diagnostics, 1001 was determined to be suitable
# as the final base model. Hence we run simulations using 1000 samples from the
# full posterior distribution across all 4 chains. To do this, we have
# constructed an mrgsolve model (`1001.mod`) and now run the `nm_join_bayes()`
# function to generate the appropriate diagnostics.
# This function generates the following simulation-based diagnostics, and write
# them to an RDS file:
#   - EPRED (median/mean and percentiles)
#   - IPRED (median/mean and percentiles, if .iph files available)
#   - NPDE (using median/mean of posterior, or full posterior)
#   - EWRES (using median/mean of posterior, or full posterior)
# It also replaces NONMEM-output ETAs with medians of posterior ETAs from .iph
# files if available.
#
# Although these simulation-based diagnostics are available from NONMEM table
# files for each individual chain, for any final model run with BAYES/NUTS it is
# required that these simulations are run using the full posterior using the
# methods implemented in this code.
mod_ms <- mread(here(glue("script/model/{model_name}.mod")), outvars = "IPRED,Y")
mod <- read_model(file.path(model_dir, model_name))
set.seed(201)
system.time({
progressr::with_progress({
  sim <- nm_join_bayes(mod, mod_ms)
})
})
# Elapsed time run on 4 cores/32 GB RAM:
#    user  system elapsed
# 133.977  24.310  69.829 
sim
saveRDS(
  sim,
  file.path(model_dir, model_name, glue("diag-sims-{model_name}.rds"))
)

# Model diagnostics
rmarkdown::render(
  here("script/diagnostic-templates/template-model-bayes.Rmd"),
  params = list(
    run = model_name,
    # We now need to also specify the path of the simulation output file
    sims_output_path = file.path(model_dir, model_name,
                                 glue("diag-sims-{model_name}.rds"))
  ),
  output_dir = file.path(model_dir, model_name),
  output_file = glue("model-diagnostics-{model_name}.html")
)
utils::browseURL(file.path(model_dir, model_name,
                           glue("model-diagnostics-{model_name}.html")))

# MCMC diagnostics
rmarkdown::render(
  here("script/diagnostic-templates/template-mcmc.Rmd"),
  params = list(run = model_name),
  output_dir = file.path(model_dir, model_name),
  output_file = glue("mcmc-diagnostics-{model_name}.html")
)
# View the output
utils::browseURL(file.path(model_dir, model_name,
                           glue("mcmc-diagnostics-{model_name}.html")))


### run 1100 - Full covariate model ----------------------------
model_name <- 1100

# Generate columns for simulation-based diagnostics
mod_ms <- mread(here(glue("script/model/{model_name}.mod")), outvars = "IPRED,Y")
mod <- read_model(file.path(model_dir, model_name))
set.seed(201)
system.time({
progressr::with_progress({
  sim <- nm_join_bayes(mod, mod_ms)
})
})
# Elapsed time run on 4 cores/32 GB RAM:
#   user  system elapsed
# 52.423  19.376  56.630 
sim
saveRDS(
  sim,
  file.path(model_dir, model_name, glue("diag-sims-{model_name}.rds"))
)

# Model diagnostics
rmarkdown::render(
  here("script/diagnostic-templates/template-model-bayes.Rmd"),
  params = list(
    run = model_name,
    sims_output_path = file.path(model_dir, model_name,
                                 glue("diag-sims-{model_name}.rds"))
  ),
  output_dir = file.path(model_dir, model_name),
  output_file = glue("model-diagnostics-{model_name}.html")
)
utils::browseURL(file.path(model_dir, model_name,
                           glue("model-diagnostics-{model_name}.html")))

# MCMC diagnostics
rmarkdown::render(
  here("script/diagnostic-templates/template-mcmc.Rmd"),
  params = list(run = model_name),
  output_dir = file.path(model_dir, model_name),
  output_file = glue("mcmc-diagnostics-{model_name}.html")
)
# View the output
utils::browseURL(file.path(model_dir, model_name,
                           glue("mcmc-diagnostics-{model_name}.html")))
