# NOTES ----------------------------------------------------


# PRELIMINARIES ----------------------------------------------------

#rm(list=ls())

# This script uses renv to preserve the R environment specs (e.g., package versions.)
library(renv)
# run this if you want to reproduce results using the R environment we had:
# renv::restore()

toLoad = c("crayon",
           "dplyr",
           "foreach",
           "doParallel",
           "metafor",
           "robumeta",
           "data.table",
           "purrr",
           "metRology",
           "fansi",
           "MetaUtility",
           "ICC",
           "cfdecomp",
           "tidyr",
           "tibble",
           "testthat",
           "rstan", # note: to reinstall this one, need to use high-mem session
           "optimx",
           "weightr",
           "phacking",
           "here")  # note: to reinstall this one, need ml load jags

# to install everything
# lapply(toLoad, install.packages)

lapply( toLoad,
        require,
        character.only = TRUE)



# prevent masking
select = dplyr::select

# run this only if you want to update the R environment specs
# setwd(here())
# renv::snapshot()

# ~~ User-specified global vars -------------------------
# no sci notation
options(scipen=999)

# control which results should be redone and/or overwritten
# but note that not all fns respect this setting
overwrite.res = TRUE


# ~~ Set directories -------------------------
( code.dir = str_replace( string = here(),
                          pattern = "Applied example",
                          replacement = "Simulation study/Code \\(git\\)") )
# test it
setwd(code.dir)

# ( data.dir = str_replace( string = here(),
#                           pattern = "Code \\(git\\)",
#                           replacement = "Results/Working dataset") )
# 
# ( results.dir = str_replace( string = here(),
#                              pattern = "Code \\(git\\)",
#                              replacement = "Results/Working results") )

# check that they're specified correctly
setwd(data.dir)
setwd(results.dir)

# below are the only absolute paths
# write results directly to directory containing TeX manuscript in Overleaf so stats can be piped directly into text
# this is an absolute path because it must live in Dropbox, outside the project directory, in order to sync with Overleaf
# to reproduce results, just set this to any directory on your local machine
# results will be written to a csv file in that location
overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/JTE (Jeffreys tau estimation) Overleaf/R_objects/figures"
overleaf.dir.stats = "/Users/mmathur/Dropbox/Apps/Overleaf/JTE (Jeffreys tau estimation) Overleaf/R_objects"


setwd(code.dir)
source("analyze_sims_helper_JTE.R")
source("helper_JTE.R")  # for lprior(), etc.
source("init_stan_model_JTE.R")


# PRELIMINARIES ----------------------------------------------------



library(multibiasmeta)
d = meta_meat  # only 34 clusters, but might work
# look at smaller dataset
d = meta_meat[1:5,]

# other metas
setwd("/Users/mmathur/Dropbox/Personal computer/Reference sheets/Library of prepped example meta-analyses/MRM")
d = fread("hu_data_prepped.csv")
d$sei = sqrt(d$vi)
d = d[1:8,]


setwd("/Users/mmathur/Dropbox/Personal computer/Reference sheets/Library of prepped example meta-analyses/SAPB")
d = fread("boehm_prepped.csv")
d$sei = sqrt(d$vi)
d = d[1:8,]



#### AWR subset of high-quality studies - works pretty well
setwd("/Users/mmathur/Dropbox/Personal computer/Reference sheets/Library of prepped example meta-analyses/MRM")
d = fread("mathur_awr_prepped.csv")


# from AWR code
d = d %>% mutate( hi.qual = (randomized == TRUE) & 
                    qual.exch %in% c("a.Low", "b.Medium") &
                    qual.sdb %in% c("a.Low", "b.Medium") &  # this is the killer
                    !is.na(qual.missing) & qual.missing < 15 )   # reducing this to 5 doesn't change

d = d %>% filter(hi.qual == TRUE)

d$authoryear
d$sei = sqrt(d$vi)

# pre-aggregate the clusters
d = d %>% group_by(authoryear) %>%
  mutate( yi_agg = rma.uni(yi = yi, sei = sei, method = "FE")$b,
          vi_agg = rma.uni(yi = yi, sei = sei, method = "FE")$se^2 )

d2 = d %>% filter( !duplicated(authoryear) ) 

.dat = as.data.frame( d2 %>% select(yi_agg, vi_agg) )
names(.dat) = c( "authoryear", "yi", "vi")
.dat$sei = sqrt(.dat$vi)
####


### Insomnia meta
# see prep script
.dat = d %>% filter(group=="Posttreatment")
.dat = d %>% filter(group=="Early follow-up")  # makes a big difference to width
.dat = d %>% filter(group=="Late follow-up")  # makes a big difference, and becomes significant




rep.res = data.frame()

rep.res = run_method_safe(method.label = c("REML"),
                          method.fn = function() {
                            mod = rma( yi = .dat$yi,
                                       vi = .dat$vi,
                                       method = "REML",
                                       knha = TRUE )
                            
                            report_meta(mod, .mod.type = "rma")
                          },
                          .rep.res = rep.res )


srr(rep.res)



## Jeffreys
rep.res = run_method_safe(method.label = c("jeffreys-pmean",
                                           "jeffreys-pmed",
                                           "jeffreys-max-lp-iterate"),
                          method.fn = function() estimate_jeffreys(.yi = as.numeric(.dat$yi),
                                                                   .sei = .dat$sei,
                                                                   
                                                                   .Mu.start = 0,
                                                                   # can't handle start value of 0:
                                                                   .Tt.start = 0.2,
                                                                   .stan.adapt_delta = 0.995,
                                                                   .stan.maxtreedepth = 25), .rep.res = rep.res )


# start values for finding posterior mode analytically
maxlp = rep.res[ rep.res$method == "jeffreys-max-lp-iterate", ]
Mhat.MaxLP = maxlp$Mhat
Shat.MaxLP = maxlp$Shat

srr(rep.res)



# find posterior mode analytically
rep.res = run_method_safe(method.label = c("Jeffreys"),
                          method.fn = function() {
                            
                            # as in the future pkg
                            mle_fit <- mle_params(Mhat.MaxLP, Shat.MaxLP, .dat$yi, .dat$sei)
                            modes <- c(mle_fit@coef[["mu"]], mle_fit@coef[["tau"]])
                            optim_converged <- mle_fit@details$convergence == 0
                            
                            
                            return( list( stats = data.frame( 
                              
                              Mhat = modes[1],
                              Shat = modes[2],
                              
                              # all inference is again from MCMC
                              MhatSE = maxlp$MhatSE,
                              ShatSE = maxlp$ShatSE,
                              MLo = maxlp$MLo,
                              MHi =maxlp$MHi,
                              SLo = maxlp$SLo,
                              SHi = maxlp$SHi,
                              
                              stan.warned = maxlp$stan.warned,
                              stan.warning = maxlp$stan.warning,
                              MhatRhat = maxlp$MhatRhat,
                              ShatRhat = maxlp$ShatRhat,
                              
                              OptimConverged = optim_converged) ) )
                            
                          }, .rep.res = rep.res )


srr(rep.res)


