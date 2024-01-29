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

# should analyses be run from scratch?
# FALSE if just redoing plots
rerun.analyses = TRUE


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



# READ PREPPED DATA ----------------------------------------------------

setwd(data.dir)
d = fread("data_insomnia.csv")


# ANALYZE EACH SUBSET ----------------------------------------------------

if ( rerun.analyses == TRUE ) {
  
  setwd(code.dir)
  source("init_stan_model_JTE.R")
  
  if (exists("rs")) rm(rs)
  
  # fit Jeffreys, DL, and REML for each group (i.e., section in their forest plot)
  for ( .group in unique(d$group) ) {
    
    # test only
    #.group = "SOL; early follow-up (k = 4)"
    
    .dat = d %>% filter(group==.group)
    
    rep.res = data.frame()
    
    ### DL
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
    
    ### REML
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
    
    ### MLE + profile interval
    rep.res = run_method_safe(method.label = c("MLE-profile"),
                              method.fn = function() {
                                
                                
                                nll_fun <- function(mu, tau) get_nll(mu, tau, .dat$yi, .dat$sei)
                                my_mle = stats4::mle(minuslogl = nll_fun,
                                                     start = list(mu = 0, tau = 0.1),
                                                     method = "L-BFGS-B")
                                
                                # this fn will complain if optimizer in my_mle hasn't found the true max
                                #  which is a good thing
                                cis = stats4::confint(my_mle)
                                
                                return( list( stats = data.frame( 
                                  
                                  Mhat = as.numeric( attr(my_mle, "coef")["mu"] ),
                                  Shat = as.numeric( attr(my_mle, "coef")["tau"] ),
                                  
                                  MhatSE = NA,
                                  ShatSE = NA,
                                  
                                  MLo = cis["mu", 1],
                                  MHi = cis["mu", 2],
                                  
                                  SLo = cis["tau", 1],
                                  SHi = cis["tau", 2] ) ) )
                                
                              },
                              .rep.res = rep.res )
    
    
    srr(rep.res)

   

    
    ### Jeffreys
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
    
    # this can happen if there is a failure
    if ("method.1" %in% names(rep.res)) rep.res = rep.res %>% select(-method.1)
    
    
    if ( .group == unique(d$group)[1] ) rs = rep.res else rs = rbind(rs, rep.res)
    
} # end loop over .group
  
  # POST-PROCESSING
  
  # remove the extra Jeffreys methods
  rsp = rs %>% filter(method %in% c( "MLE-profile", "REML", "DL", "Jeffreys") )
  
  # ratio of CI width of Jeffreys vs. the winner among the other two methods
  rsp = rsp %>% group_by(group) %>%
    mutate( MhatWidth = MHi - MLo, 
            MhatTestReject = sign(MHi) == sign(MLo),
            CI_ratio = min( MhatWidth[ method != "Jeffreys" ] ) / MhatWidth[ method == "Jeffreys" ] )
  
  # OnlyJeffreysRejects = MhatTestReject[ method != "Jeffreys" ] == 0 &
  #   MhatTestReject[ method == "Jeffreys" ] == 1)
  
  # extract k, which is the numeric part of the group variable
  rsp$k = as.numeric( str_extract(rsp$group, "\\d+") )
  
  View(rsp %>% select(method, Mhat, MLo, MHi, MhatTestReject, CI_ratio))
  
  
  setwd(results.dir)
  fwrite(rsp, "insomnia_forest_results.csv")
  
}  # end "if (rerun.analyses == TRUE)"




# EXPLORE OTHER METHODS  -------------------------------------------------

.group = "TST; early follow-up (k = 4)"
d2 = d %>% filter( group == .group )


# ~ MLE-profile  -------------------------------------------------

mle_params2 <- function(mu_start, tau_start, yi, sei) {
  nll_fun <- function(mu, tau) get_nll(mu, tau, yi, sei)
  stats4::mle(minuslogl = nll_fun,
              start = list(mu = mu_start, tau = tau_start),
              # using this method because in the single k=4 example below, 
              #  seems more likely than other methods to converge
              method = "L-BFGS-B")
}

my_mle = mle_params2(mu_start = 0,
                     tau_start = 0.1, 
                     yi = d2$yi,
                     sei = d2$sei)
my_mle
# warning that original solution hadn't converged:
cis = stats4::confint(my_mle) 


od = get_optimx_dataframe(.yi = d2$yi, 
                          .sei = d2$sei,
                          .mu.start = 0,
                          .tau.start = 0.1)


# which method won?
od2 = od %>% select( all_of( namesWith(dat = od, pattern = "nll") ) )
min(od2)


### Try with metaLik package
library(metaLik)
mod2 = metaLik( yi ~ 1,
                data = d2,
                sigma2 = vi)

summ = summary(mod2)

confint( summary(mod2) )
inf = test.metaLik(mod2, param=1)

# try their example meta-analysis
data(education)
m <- metaLik(y~1, data=education, sigma2=sigma2)
summary(m)
confint(m)
test.metaLik(m, param=1)

profile(m)

#bm: how did the Annals authors get CIs from this package??




### Why does my MLE-profile not agree with metafor for tau specifically?
# is it because of REML vs. ML?

# can we just use metafor?
# no: it won't profile mu
# https://wviechtb.github.io/metafor/reference/profile.rma.html
# some underlying fns: https://github.com/cran/metafor/blob/master/R/profile.rma.uni.r
# uses Q-profile intervals
# from here (https://wviechtb.github.io/metafor/reference/confint.rma.html):
# For objects of class "rma.uni" obtained with the rma.uni function, a confidence interval for the amount of (residual) heterogeneity (i.e., 𝜏2
# ) can be obtained by setting random=TRUE (which is the default). The interval is obtained iteratively either via the Q-profile method or via the generalized Q-statistic method (Hartung and Knapp, 2005; Viechtbauer, 2007; Jackson, 2013; Jackson et al., 2014). The latter is automatically used when the model was fitted with method="GENQ" or method="GENQM", the former is used in all other cases.
#**this DOES change if original model is fit using REML
mod = rma( yi = d2$yi,
           vi = d2$vi,
           method = "ML",
           knha = FALSE )
( x = confint(mod) )
( tau.lb = x$random["tau","ci.lb"] )
( tau.ub = x$random["tau","ci.ub"] )
data.frame( tau.lb, tau.ub )


# my profile CI, as in doParallel
nll_fun <- function(mu, tau) get_nll(mu, tau, d2$yi, d2$sei)
my_mle = stats4::mle(minuslogl = nll_fun,
                     start = list(mu = 0, tau = 0.1),
                     method = "L-BFGS-B")

# this fn will complain if optimizer in my_mle hasn't found the true max
#  which is a good thing
cis = stats4::confint(my_mle)
cis
#??? why is the CI for tau symmetric around 0? doesn't make sense.

prof = stats4::profile(my_mle) 
plot( prof )


x = as.data.frame( attr(prof, "profile")$tau )
# okay, so the issue is just that it's trying to profile values of tau that are <0
plot( x$par.vals[, "tau"],
      x$z )



# ~ Exact -------------------------------------------------

### try rma.exact
library(rma.exact)
rma.exact(yi = d2$yi,
          vi = d2$vi)

# ~ bayesmeta  -------------------------------------------------

### Bayesmeta with Jeffreys prior on tau alone
library(bayesmeta)

# sanity check: should match MLE
#  because when mu.prior isn't specified, defaults to uniform
# yes, matches :)
m = bayesmeta(y = d2$yi,
              sigma = d2$sei,
              tau.prior = "uniform")

m = bayesmeta(y = d2$yi,
              sigma = d2$sei,
              tau.prior = "Jeffreys")
m
m$MAP["joint", "mu"]
as.numeric( m$post.interval(tau.level=0.95) )

m$post.interval(mu.level=0.95)
m$post.interval(mu.level=0.95, method = "central")



# ~ Compare to other methods  -------------------------------------------------

rs %>% filter(group == .group)


# FOREST PLOT ----------------------------------------------------

# retrieve existing analysis results
if ( rerun.analyses == FALSE ) {
  setwd(results.dir)
  rsp = fread("insomnia_forest_results.csv")
}


# set y-axis order
unique(rsp$group)
correct.order = c( "SOL; posttreatment (k = 16)",
                   "SOL; early follow-up (k = 4)",
                   "SOL; late follow-up (k = 3)",
                   
                   "WASO; posttreatment (k = 14)",
                   "WASO; early follow-up (k = 3)",
                   "WASO; late follow-up (k = 3)",
                   
                   "TST; posttreatment (k = 16)",
                   "TST; early follow-up (k = 4)",
                   "TST; late follow-up (k = 4)",
                   
                   "SE%; posttreatment (k = 17)",
                   "SE%; early follow-up (k = 5)",
                   "SE%; late follow-up (k = 4)" )
  
rsp$group = factor( rsp$group, levels = rev(correct.order) )
levels(rsp$group)

# reorder methods
correct.order = c("MLE-profile", "DL", "REML", "Jeffreys")
rsp$method = factor(rsp$method, levels = correct.order)
levels(rsp$method)

# same colors as in analyze_sims_helper.R / prior_plot_one_k for prettiness

.colors = c("black",
            "#0E96F0",
            "#0F5A8C",
            "#F2340E")


# find good x-axis limits
min(rsp$MLo)
max(rsp$MHi)
xmin = -100
xmax = 80

p = ggplot( data = rsp,
            aes( y = group,
                 x = Mhat, 
                 xmin = MLo, 
                 xmax = MHi,
                 color = method) ) +
  
  # reference line at null
  geom_vline(xintercept = 0,
             lwd = .8,
             color = "gray") +
  
  
  geom_errorbarh( aes(xmax = MHi,
                     xmin = MLo),
                 height = 0,
                 lwd = 0.8,
                 position = position_dodge(width = 0.5) ) +
  
  
  geom_point(size=3,
             position=position_dodge(width = 0.5) ) +
  
  # manually provided colors
  scale_colour_manual(values = .colors,
                      guide = guide_legend(reverse = TRUE)) +
  
  
  scale_y_discrete( name = "Study subset (outcome; follow-up duration)" ) +
  scale_x_continuous( limits = c(xmin - 2, xmax + 2),
                      breaks = seq(xmin, xmax, 20) ) +
  
  
  xlab( "Pooled estimate with 95% CI" ) +
  
  labs(color  = "Method") +
  
  
  theme_bw(base_size = 16) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=20),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" )

p



my_ggsave(name = "insomnia_forest.pdf",
          .plot = p,
          .width = 10,
          .height = 13,
          .results.dir = results.dir,
          .overleaf.dir = overleaf.dir.figs)


# ONE-OFF STATS FOR PAPER  -------------------------------------------------


# retrieve existing analysis results
if ( rerun.analyses == FALSE ) {
  setwd(results.dir)
  rsp = fread("insomnia_forest_results.csv")
}

# ~ CI width comparisons  -------------------------------------------------
#@definitely check these and the underlying CI_ratio calculation
update_result_csv( name = "Mean perc narrower Jeffreys vs winning other method",
                   value = round( 100 * ( mean(rsp$CI_ratio) - 1 ) ),
                   print = TRUE )

update_result_csv( name = "Mean perc narrower Jeffreys vs winning other method - small metas",
                   value = round( 100 * ( mean(rsp$CI_ratio[ rsp$k <= 5] ) - 1 ) ),
                   print = TRUE )

update_result_csv( name = "Mean perc wider Jeffreys vs winning other method - larger metas",
                   value = round( 100 * ( mean( 1 / rsp$CI_ratio[ rsp$k > 5] ) - 1 ) ),
                   print = TRUE )

# ~ MhatTestReject comparisons -------------------------------------------------
update_result_csv( name = "Mean MhatTestReject Jeffreys",
                   value = round( 100 * mean( rsp$MhatTestReject[ rsp$method == "Jeffreys"] ) ),
                   print = TRUE )

update_result_csv( name = "Mean MhatTestReject REML",
                   value = round( 100 * mean( rsp$MhatTestReject[ rsp$method == "REML"] ) ),
                   print = TRUE )

update_result_csv( name = "Mean MhatTestReject DL",
                   value = round( 100 * mean( rsp$MhatTestReject[ rsp$method == "DL"] ) ),
                   print = TRUE )



update_result_csv( name = "Mean MhatTestReject Jeffreys - small metas",
                   value = round( 100 * mean( rsp$MhatTestReject[ rsp$method == "Jeffreys" & rsp$k <= 5] ) ),
                   print = TRUE )

update_result_csv( name = "Mean MhatTestReject REML - small metas",
                   value = round( 100 * mean( rsp$MhatTestReject[ rsp$method == "REML" & rsp$k <= 5] ) ),
                   print = TRUE )

update_result_csv( name = "Mean MhatTestReject DL - small metas",
                   value = round( 100 * mean( rsp$MhatTestReject[ rsp$method == "DL"& rsp$k <= 5] ) ),
                   print = TRUE )


update_result_csv( name = "Mean MhatTestReject Jeffreys - larger metas",
                   value = round( 100 * mean( rsp$MhatTestReject[ rsp$method == "Jeffreys" & rsp$k > 5] ) ),
                   print = TRUE )

update_result_csv( name = "Mean MhatTestReject REML - larger metas",
                   value = round( 100 * mean( rsp$MhatTestReject[ rsp$method == "REML" & rsp$k > 5] ) ),
                   print = TRUE )

update_result_csv( name = "Mean MhatTestReject DL - larger metas",
                   value = round( 100 * mean( rsp$MhatTestReject[ rsp$method == "DL"& rsp$k > 5] ) ),
                   print = TRUE )

# # basic info about meta-analysis
# update_result_csv( name = "Insomnia meta k",
#                    value = nrow(d),
#                    print = TRUE )


# ~ Rhat's -------------------------------------------------

update_result_csv( name = "Max MhatRhat - applied",
                   value = max( rsp$MhatRhat[ rsp$method == "Jeffreys"] ),
                   print = TRUE )

update_result_csv( name = "Max ShatRhat - applied",
                   value = max( rsp$ShatRhat[ rsp$method == "Jeffreys"] ),
                   print = TRUE )



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

