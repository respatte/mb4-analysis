# LIBRARY IMPORTS ==================================================================================
library(brms)
library(bridgesampling)
library(tidyverse)
library(scales)
library(beepr)
library(future)
library(future.apply)
plan(multiprocess, workers = 4) # Adapt to the number of cores you want to use

# GENERATE DATA ====================================================================================
# Define simulation-specific variables
n_labs <- 20
infants_by_lab <- rnorm(n_labs, 35, 20) %>% round() %>% pmax.int(16) # Mean and SD from MB1, min 16
n_infants <- sum(infants_by_lab)
age_min <- 165
age_max <- 320
pr_helper <- 0.64     # From Margoni & Surian 2018 (estimated true effect size)
pr_helper_sd <- 0.03  # From Margoni & Surian 2018 (confidence interval)

# Generate repeat datasets
n_sims <- 15 # Running multiple sims can lead to brms errors with save_all_pars, resulting in NAs
data_sims <- future_replicate(n_sims,
                              tibble(lab_id = mapply(rep,
                                                     paste0("lab",1:n_labs),
                                                     infants_by_lab) %>%
                                       unlist() %>%
                                       as_factor(),
                                     age_days = sample(age_min:age_max,
                                                       n_infants,
                                                       replace = T),
                                     z_age_days = scale(age_days),
                                     chose_helper = rbinom(n_infants, 1,
                                                           rnorm(n_infants,
                                                                 mean = pr_helper,
                                                                 sd = pr_helper_sd))),
                              simplify = F)
# Define informative priors for Bayesian models (sampling fails otherwise)
# Same Intercept whether raw or scaled age
# Same mildly informative prior on age_days whether raw or scaled age 
priors.full <- c(set_prior("normal(.5753641, .1)", # From Margoni & Surian, through logit
                           class = "Intercept"),
                 set_prior("normal(0, .5)",
                           class = "b"),
                 set_prior("student_t(3, 0, 2)",
                           class = "sd"))
priors.nointercept <- c(set_prior("normal(0, .5)",
                                  class = "b"),
                        set_prior("student_t(3, 0, 2)",
                                  class = "sd"))

# RUN STATISTICS: RAW AGE ==========================================================================
save_path <- "simulations_results/raw_age/"
# Run bayesian models and bridge-sample
run_models <- T
if(run_models){
  ## Full model
  ### Run models
  mods.raw_age.full <- lapply(data_sims,
                              function(df){
                                tryCatch(
                                  {
                                    brm(chose_helper ~ age_days +
                                          (age_days | lab_id),
                                        data = df,
                                        family = bernoulli,
                                        prior = priors.full,
                                        iter = 10000,
                                        future = T,
                                        save_all_pars = T,
                                        control = list(adapt_delta = .9999,
                                                       max_treedepth = 20))
                                  },
                                  error = function(cond){return(NA)})
                              })
  ### Bridge-sample posteriors
  bridges.raw_age.full <- lapply(mods.raw_age.full,
                                 function(m){
                                   tryCatch(bridge_sampler(m, silent = T),
                                            error = function(cond){return(NA)})
                                 })
  ## No Intercept
  ### Run models
  mods.raw_age.nointercept <- lapply(data_sims,
                                     function(df){
                                       tryCatch(
                                         {
                                           brm(chose_helper ~ 0 + age_days +
                                                 (0 + age_days | lab_id),
                                               data = df,
                                               family = bernoulli,
                                               prior = priors.nointercept,
                                               iter = 10000,
                                               future = T,
                                               save_all_pars = T,
                                               control = list(adapt_delta = .9999,
                                                              max_treedepth = 20))
                                         },
                                         error = function(cond){return(NA)})
                                     })
  ### Bridge-sample posteriors
  bridges.raw_age.nointercept <-lapply(mods.raw_age.nointercept,
                                       function(m){
                                         tryCatch(bridge_sampler(m, silent = T),
                                                  error = function(cond){return(NA)})
                                       })
  ## Save all
  lapply(1:n_sims,
         function(i){
           saveRDS(mods.raw_age.full[[i]],
                   paste0(save_path, "model_full_",i,".rds"))
           saveRDS(mods.raw_age.nointercept[[i]],
                   paste0(save_path, "model_nointercept_",i,".rds"))
           saveRDS(bridges.raw_age.full[[i]],
                   paste0(save_path, "bridge_full_",i,".rds"))
           saveRDS(bridges.raw_age.nointercept[[i]],
                   paste0(save_path, "bridge_nointercept_",i,".rds"))
         })
}else{
  ## Read all
  mods.raw_age.full <- lapply(1:n_sims,
                              function(i){
                                readRDS(paste0(save_path, "model_full_",i,".rds"))
                              })
  mods.raw_age.nointercept <- lapply(1:n_sims,
                                     function(i){
                                       readRDS(paste0(save_path, "model_nointercept_",i,".rds"))
                                     })
  bridges.raw_age.full <- lapply(1:n_sims,
                                 function(i){
                                   readRDS(paste0(save_path, "bridge_full_",i,".rds"))
                                 })
  bridges.raw_age.nointercept <- lapply(1:n_sims,
                                        function(i){
                                          readRDS(paste0(save_path, "bridge_nointercept_",i,".rds"))
                                        })
}
beep(8)

# Get Bayes factors
bf.raw_age <- lapply(1:n_sims,
                     function(i){
                       tryCatch(
                         {bayes_factor(bridges.raw_age.full[[i]],
                                       bridges.raw_age.nointercept[[i]])
                         },
                         error = function(cond){return(NA)})
                     })

# RUN STATISTICS: SCALED AGE =======================================================================
save_path <- "simulations_results/scaled_age/"
# Run bayesian models and bridge-sample
run_models <- T
if(run_models){
  ## Full model
  ### Run models
  mods.scaled_age.full <- lapply(data_sims,
                                 function(df){
                                   tryCatch(
                                     {
                                       brm(chose_helper ~ z_age_days +
                                             (z_age_days | lab_id),
                                           data = df,
                                           family = bernoulli,
                                           prior = priors.full,
                                           iter = 10000,
                                           save_all_pars = T,
                                           control = list(adapt_delta = .9999,
                                                          max_treedepth = 20))
                                     },
                                     error = function(cond){return(NA)})
                                 })
  ### Bridge-sample posteriors
  bridges.scaled_age.full <- lapply(mods.scaled_age.full,
                                    function(m){
                                      tryCatch(bridge_sampler(m, silent = T),
                                               error = function(cond){return(NA)})
                                    })
  ## No Intercept
  ### Run models
  mod.scaled_age.nointercept <- lapply(data_sims,
                                       function(df){
                                         tryCatch(
                                           {
                                             brm(chose_helper ~ 0 + z_age_days +
                                                   (0 + z_age_days | lab_id),
                                                 data = df,
                                                 family = bernoulli,
                                                 prior = priors.nointercept,
                                                 iter = 10000,
                                                 save_all_pars = T,
                                                 control = list(adapt_delta = .9999,
                                                                max_treedepth = 20))
                                           },
                                           error = function(cond){return(NA)})
                                       })
  ### Bridge-sample posteriors
  bridge.scaled_age.nointercept <- lapply(mods.scaled_age.nointercept,
                                          function(m){
                                            tryCatch(bridge_sampler(m, silent = T),
                                                     error = function(cond){return(NA)})
                                          })
  ## Save all
  saveRDS(mod.scaled_age.full, paste0(save_path, "model_full.rds"))
  saveRDS(mod.scaled_age.nointercept, paste0(save_path, "model_nointercept.rds"))
  saveRDS(bridge.scaled_age.full, paste0(save_path, "bridge_full.rds"))
  saveRDS(bridge.scaled_age.nointercept, paste0(save_path, "bridge_nointercept.rds"))
}else{
  ## Read all
  mod.scaled_age.full <- readRDS(paste0(save_path, "model_full.rds"))
  mod.scaled_age.nointercept <- readRDS(paste0(save_path, "model_nointercept.rds"))
  bridge.scaled_age.full <- readRDS(paste0(save_path, "bridge_full.rds"))
  bridge.scaled_age.nointercept <- readRDS(paste0(save_path, "bridge_nointercept.rds"))
}
beep(8)

# Get Bayes factors
bf.scaled_age <- lapply(1:n_sims,
                        function(i){
                          tryCatch(
                            {
                              bayes_factor(bridges.scaled_age.full[[i]],
                                           bridges.scaled_age.nointercept[[i]])
                            },
                            error = function(cond){return(NA)})
                        })
