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
           "here",
           "stringr")  # note: to reinstall this one, need ml load jags

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

# need helper fns from simulation study
( code.dir = str_replace( string = here(),
                          pattern = "Applied example",
                          replacement = "Simulation study/Code \\(git\\)") )
# test it
setwd(code.dir)

data.dir = here()
results.dir = here("Results from R")

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


# READ PREPPED DATA ----------------------------------------------------

setwd(data.dir)
d = fread("data_insomnia.csv")


# ANALYZE EACH SUBSET ----------------------------------------------------

if (exists("rs")) rm(rs)

# fit REML and Jeffreys for each group (i.e., section in their forest plot)
for ( .group in unique(d$group) ) {
  
  .dat = d %>% filter(group==.group)
  
  rep.res = data.frame()
  
  
  rep.res = run_method_safe(method.label = c("DL"),
                            method.fn = function() {
                              mod = rma( yi = .dat$yi,
                                         vi = .dat$vi,
                                         method = "DL",
                                         knha = TRUE )
                              
                              report_meta(mod, .mod.type = "rma")
                            },
                            .rep.res = rep.res )
  
  
  srr(rep.res)
  
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
  
  rep.res$group = .group
  
  if ( .group == unique(d$group)[1] ) rs = rep.res else rs = rbind(rs, rep.res)
  
}

# remove the extra Jeffreys methods
rsp = rs %>% filter(method %in% c("REML", "DL", "Jeffreys") )



setwd(results.dir)
fwrite(rsp, "insomnia_forest_results.csv")


# sanity check: do REML estimates agree with published forest plots?
# close, but not exact
rsp %>% filter( method == "REML" ) %>%
  select( group, Mhat, MLo, MHi )

# maybe they used DL instead?
rsp %>% filter( method == "DL" ) %>%
  select( group, Mhat, MLo, MHi )




# FOREST PLOT ----------------------------------------------------

# set y-axis order
correct.order = c( "SOL; posttreatment (k = 16)",
                   "SOL; early follow-up (k = 4)",
                   "SOL; late follow-up (k = 3)",
                   
                   "WASO; posttreatment (k = 14)",
                   "WASO; early follow-up (k = 3)",
                   "WASO; late follow-up (k = 3)",
                   
                   "TST; posttreatment (k = 16)",
                   "TST; early follow-up (k = 4)",
                   "TST; late follow-up (k = 4)",
                   
                   "SE%; posttreatment (k = 16)",
                   "SE%; early follow-up (k = 4)",
                   "SE%; late follow-up (k = 4)" )
  
rsp$group = factor( rsp$group, levels = rev(correct.order) )
levels(rsp$group)



.colors = c("orange", "black")

p = ggplot( data = rsp,
            aes( x = group,
                 y = Mhat, 
                 ymin = MLo, 
                 ymax = MHi,
                 color = method) ) +
  
  # reference line at null
  geom_hline(yintercept = 0,
             lwd = .8,
             color = "gray") +
  
  
  geom_errorbar( aes(ymax = MHi,
                      ymin = MLo),
                  width = 0,
                  lwd = 0.8,
                  position = position_dodge(width = 0.5) ) +
  

  geom_point(size=3,
             position=position_dodge(width = 0.5) ) +
  
  # manually provided colors
  scale_colour_manual(values = .colors ) +
  

  scale_x_discrete( name = "Subset" ) +
  #scale_y_continuous(name="Odds ratio", limits = c(0.5, 5)) +
  coord_flip() +
  
  
  ylab( "Pooled estimate with 95% CI" ) +

  labs(color  = "Method") +

  
  theme_bw() +
  
  theme( text = element_text(face = "bold"),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         legend.position = "bottom" )


p


#bm: make sure REML results all agree with the original paper
# also work on plot order :)
# report: For what percent of metas < 10 studies was Jeffreys more precise? and for what percent of larger ones?
# look into how I could report p-values









# # TRASH - other possible examples
# 
# library(multibiasmeta)
# d = meta_meat  # only 34 clusters, but might work
# # look at smaller dataset
# d = meta_meat[1:5,]
# 
# # other metas
# setwd("/Users/mmathur/Dropbox/Personal computer/Reference sheets/Library of prepped example meta-analyses/MRM")
# d = fread("hu_data_prepped.csv")
# d$sei = sqrt(d$vi)
# d = d[1:8,]
# 
# 
# setwd("/Users/mmathur/Dropbox/Personal computer/Reference sheets/Library of prepped example meta-analyses/SAPB")
# d = fread("boehm_prepped.csv")
# d$sei = sqrt(d$vi)
# d = d[1:8,]
# 
# 
# 
# #### AWR subset of high-quality studies - works pretty well
# setwd("/Users/mmathur/Dropbox/Personal computer/Reference sheets/Library of prepped example meta-analyses/MRM")
# d = fread("mathur_awr_prepped.csv")
# 
# 
# # from AWR code
# d = d %>% mutate( hi.qual = (randomized == TRUE) & 
#                     qual.exch %in% c("a.Low", "b.Medium") &
#                     qual.sdb %in% c("a.Low", "b.Medium") &  # this is the killer
#                     !is.na(qual.missing) & qual.missing < 15 )   # reducing this to 5 doesn't change
# 
# d = d %>% filter(hi.qual == TRUE)
# 
# d$authoryear
# d$sei = sqrt(d$vi)
# 
# # pre-aggregate the clusters
# d = d %>% group_by(authoryear) %>%
#   mutate( yi_agg = rma.uni(yi = yi, sei = sei, method = "FE")$b,
#           vi_agg = rma.uni(yi = yi, sei = sei, method = "FE")$se^2 )
# 
# d2 = d %>% filter( !duplicated(authoryear) ) 
# 
# .dat = as.data.frame( d2 %>% select(yi_agg, vi_agg) )
# names(.dat) = c( "authoryear", "yi", "vi")
# .dat$sei = sqrt(.dat$vi)
# ####

