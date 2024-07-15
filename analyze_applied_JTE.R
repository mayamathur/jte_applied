
# PRELIMINARIES ----------------------------------------------------

# rm(list=ls())

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
           "stringr",
           "bayesmeta",
           "rma.exact")  

# to install everything
lapply( toLoad,
        require,
        character.only = TRUE)

# prevent masking
select = dplyr::select

# run this only if you want to update the R environment specs
# setwd(here())
# renv::snapshot()

# no sci notation
options(scipen=999)


# ~~ User-specified global vars -------------------------

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
d = fread("data_zito.csv")


# ANALYZE EACH SUBSET ----------------------------------------------------

if ( rerun.analyses == TRUE ) {
  
  if (exists("rs")) rm(rs)
  
  # fit Jeffreys, DL, and REML for each group (i.e., section in their forest plot)
  for ( .group in unique(d$group) ) {
    
    # test only
    # .group = "CV death and myocardial infarction (k = 2)"
    
    .dat = d %>% filter(group==.group)
    
    rep.res = data.frame()
    
    ### DL
    rep.res = run_method_safe(method.label = c("DL-HKSJ"),
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
    rep.res = run_method_safe(method.label = c("REML-HKSJ"),
                              method.fn = function() {
                                mod = rma( yi = .dat$yi,
                                           vi = .dat$vi,
                                           method = "REML",
                                           knha = TRUE )
                                
                                report_meta(mod, .mod.type = "rma")
                              },
                              .rep.res = rep.res )
    
    
    srr(rep.res)

    
    ### Exact
    rep.res = run_method_safe(method.label = c("Exact"),
                              method.fn = function() {
                                
                                ci = rma.exact.fast( yi = .dat$yi,
                                                     vi = .dat$vi,
                                                     plot = FALSE )
                                
                                # this method doesn't do point estimation of inference for tau
                                return( list( stats = data.frame( 
                                  MLo = ci[1],
                                  MHi = ci[2]) ) )
                                
                                
                              },
                              .rep.res = rep.res )
    
    
    ### Jeffreys1-central
    rep.res = run_method_safe(method.label = c("Jeffreys1-shortest"),
                              method.fn = function() {
                                
                                m = bayesmeta(y = .dat$yi,
                                                sigma = .dat$sei,
                                                tau.prior = "Jeffreys",
                                                interval.type = "shortest")
                                
                                # sanity check: plot posterior and prior
                                # prior is the dashed line
                                # plot(m, prior = TRUE)
                                
                                # marginal (not joint) intervals
                                tau_ci = as.numeric( m$post.interval(tau.level=0.95) ) 
                                mu_ci = as.numeric( m$post.interval(mu.level=0.95) )
                                
                                # this method doesn't do point estimation of inference for tau
                                return( list( stats = data.frame( 
                                  Mhat = m$MAP["marginal", "mu"],
                                  Shat = m$MAP["marginal", "tau"],
                                  MLo = mu_ci[1],
                                  MHi = mu_ci[2],
                                  SLo = tau_ci[1],
                                  SHi = tau_ci[2] ) ) )
                                
                                
                              },
                              .rep.res = rep.res )
    
    
    ### Jeffreys2
    rep.res = run_method_safe(method.label = c("Jeffreys2-shortest"),
                              method.fn = function() {
                                
                                m = bayesmeta(y = .dat$yi,
                                                sigma = .dat$sei,
                                                tau.prior = "overallJeffreys",
                                                interval.type = "shortest")
                                
                                # sanity check: plot posterior and prior
                                # prior is the dashed line
                                # plot(m, prior = TRUE)
                                
                                # marginal (not joint) intervals
                                tau_ci = as.numeric( m$post.interval(tau.level=0.95) ) 
                                mu_ci = as.numeric( m$post.interval(mu.level=0.95) )
                                
                                # this method doesn't do point estimation of inference for tau
                                return( list( stats = data.frame( 
                                  Mhat = m$MAP["marginal", "mu"],
                                  Shat = m$MAP["marginal", "tau"],
                                  MLo = mu_ci[1],
                                  MHi = mu_ci[2],
                                  SLo = tau_ci[1],
                                  SHi = tau_ci[2] ) ) )
                                
                                
                              },
                              .rep.res = rep.res )
    
    # this can happen if there is a failure
    if ("method.1" %in% names(rep.res)) rep.res = rep.res %>% select(-method.1)
    
    rep.res = rep.res %>% add_column(.before = 1, group = .group)
    
    
    if ( .group == unique(d$group)[1] ) rs = rep.res else rs = rbind(rs, rep.res)
    
  } # end loop over .group
  
  
  
  # POST-PROCESSING
  
  # remove unused methods
  #rsp = rs %>% filter(method %in% c( "MLE-profile", "REML", "DL", "Jeffreys") )
  
  rsp = rs %>% rename(method.pretty = method)
  
  # ratio of CI width of Jeffreys vs. the winner among the other two methods
  rsp = rsp %>% group_by(group) %>%
    mutate( MhatWidth = MHi - MLo, 
            MhatTestReject = sign(MHi) == sign(MLo),
            CI_ratio = min( MhatWidth[ method.pretty != "Jeffreys2-shortest" ], na.rm = TRUE ) / MhatWidth[ method.pretty == "Jeffreys2-shortest" ] )
  
  # extract k, which is the numeric part of the group variable
  rsp$k = as.numeric( str_extract(rsp$group, "\\d+") )
  
  View(rsp %>% select(method.pretty, Mhat, MLo, MHi, MhatTestReject, CI_ratio))
  
  
  setwd(results.dir)
  fwrite(rsp, "zito_forest_results.csv")
  
}  # end "if (rerun.analyses == TRUE)"


# sanity check for CI_ratio calculation
.group = "All-cause death (k = 3)"
temp = rsp %>% filter(group == .group)
temp %>% select(method.pretty, MhatWidth)
# narrowest/Jeffreys2:
expect_equal( round(1.04/2.21, 2), round(temp$CI_ratio[1], 2) )
  


# PRIOR AND POST PLOTS FOR ONE META-ANALYSIS  -------------------------------------------------

# pick one subset meta-analysis
.group = "All-cause death (k = 3)"
.dat = d %>% filter(group == .group)


# ~ Plot posteriors  -------------------------------------------------

m2 = bayesmeta(y = .dat$yi,
               sigma = .dat$sei,
               tau.prior = "overallJeffreys",
               interval.type = "shortest")

mu_vec = seq( -4, 4, 0.01 )
tau_vec = c( seq(0, 0.1, 0.01), seq(0.1, 0.25, 0.01), seq(0.25, 2.5, 0.05) )


dp = expand_grid( mu = mu_vec,
                  tau = tau_vec )
nrow(dp)

dp = dp %>% rowwise() %>%
  mutate( joint_post = m2$dposterior(mu = mu, tau = tau),
          mu_post = m2$dposterior(mu = mu),
          tau_post = m2$dposterior(tau = tau) )


### Posterior plot
p = ggplot(data = dp,
           aes(x = tau,
               y = mu,
               z = joint_post)) +
  
  # posterior mode for mu and 95% CI
  geom_hline(yintercept = m2$MAP["marginal", "mu"], lty = 1, color = "red") +
  # geom_vline(xintercept = m2$summary["95% lower", "mu"], lty = 2, color = "red") +
  # geom_vline(xintercept = m2$summary["95% upper", "mu"], lty = 2, color = "red") +
  
  # posterior mode for tau and 95% CI
  geom_vline(xintercept = m2$MAP["marginal", "tau"], lty = 1, color = "blue") +
  # geom_hline(yintercept = m2$summary["95% lower", "tau"], lty = 2, color = "blue") +
  # geom_hline(yintercept = m2$summary["95% upper", "tau"], lty = 2, color = "blue") +
  
  geom_contour(color = "black") +
  
  xlab( bquote( p(tau ~ "|" ~ hat(theta) ) ) ) +
  ylab( bquote( p(mu ~ "|" ~ hat(theta) ) ) ) +

  
  scale_x_continuous( limits = c(0, 0.6), breaks = seq(0, 0.6, .1) ) +
  scale_y_continuous( limits = c(-0.9, 0.3), breaks = seq(-1, 1, .1) ) +
  
  
  theme_bw(base_size = 20) +

  theme(text = element_text(face = "bold"),
        axis.title = element_text(size=20),
        legend.position = "bottom" )
p

my_ggsave(name = "zito_all_cause_death_posterior.pdf",
          .plot = p,
          .width = 10,
          .height = 8,
          .results.dir = results.dir,
          .overleaf.dir = overleaf.dir.figs)


### Marginal posterior: mu
# see Rover (2020), pg 28
p = ggplot(data = dp,
           aes(x = mu,
               y = mu_post)) +
  
  # posterior mode for mu and 95% CI
  geom_vline(xintercept = m2$MAP["marginal", "mu"], lty = 1, color = "red") +
  geom_vline(xintercept = m2$summary["95% lower", "mu"], lty = 2, color = "red") +
  geom_vline(xintercept = m2$summary["95% upper", "mu"], lty = 2, color = "red") +
  
  geom_line(linewidth = 1.1) +

  xlab( bquote(mu) ) +
  ylab( bquote( p(mu ~ "|" ~ hat(theta) ) ) ) +
  
  
  scale_x_continuous( limits = c(-3, 3), breaks = seq(-3, 3, 1) ) +
  
  theme_bw(base_size = 24) +
  
  theme(text = element_text(face = "bold"),
        axis.title = element_text(size=24),
        legend.position = "bottom",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )
p

my_ggsave(name = "zito_all_cause_death_mu_posterior.pdf",
          .plot = p,
          .width = 10,
          .height = 8,
          .results.dir = results.dir,
          .overleaf.dir = overleaf.dir.figs)



### Marginal posterior: tau
# see Rover (2020), pg 28
p = ggplot(data = dp,
           aes(x = tau,
               y = tau_post)) +
  
  # posterior mode for tau and 95% CI
  geom_vline(xintercept = m2$MAP["marginal", "tau"], lty = 1, color = "blue") +
  geom_vline(xintercept = m2$summary["95% lower", "tau"], lty = 2, color = "blue") +
  geom_vline(xintercept = m2$summary["95% upper", "tau"], lty = 2, color = "blue") +

  geom_line(linewidth = 1.1) +
  
  
  xlab( bquote(tau) ) +
  ylab( bquote( p(tau ~ "|" ~ hat(theta) ) ) ) +
  
  
  scale_x_continuous( limits = c(0, 2.5), breaks = seq(0, 2.5, .5) ) +
  
  theme_bw(base_size = 24) +
  
  theme(text = element_text(face = "bold"),
        axis.title = element_text(size=24),
        legend.position = "bottom",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )
p

my_ggsave(name = "zito_all_cause_death_tau_posterior.pdf",
          .plot = p,
          .width = 10,
          .height = 8,
          .results.dir = results.dir,
          .overleaf.dir = overleaf.dir.figs)





# ~ Plot priors: Make plotting dataframe  -------------------------------------------------

# similar to plot_prior_one_k:

# evaluate prior for various parameters
# prior is independent of .mu, so just choose one
tau_vec = c( seq(0, 0.1, 0.001), seq(0.1, 0.25, 0.01), seq(0.25, 1, 0.05) )
dp2 = expand_grid( .mu = 0,
                   .tau = tau_vec,
                   prior_name = c("Jeffreys1", "Jeffreys2" ) )

dp2 = dp2 %>% rowwise() %>%
  mutate( prior.val = exp( get_lprior(mu = .mu,
                                      tau = .tau,
                                      sei = .dat$sei,
                                      prior_name = prior_name ) ) )


# ~ Rescale the priors to have the same max height  -------------------------------------------------
dp2$prior.val.scaled = NA
dp2$prior.val.scaled[ dp2$prior_name == "Jeffreys1" ] = dp2$prior.val[ dp2$prior_name == "Jeffreys1" ] / max(dp2$prior.val[ dp2$prior_name == "Jeffreys1" ])

dp2$prior.val.scaled[ dp2$prior_name == "Jeffreys2" ] = dp2$prior.val[ dp2$prior_name == "Jeffreys2" ] / max(dp2$prior.val[ dp2$prior_name == "Jeffreys2" ])


### Points that maximize each curve
dp2 = dp2 %>% group_by(prior_name) %>%
  mutate( is.max = ifelse( prior.val.scaled == max(prior.val.scaled), TRUE, FALSE ) )
( max.points = dp2 %>% filter(is.max == TRUE) )


# ~ Make plot  -------------------------------------------------

# in same order as N.pretty
my.colors = c("#E075DB", "#F2340E")

plot = ggplot( data = dp2, 
               aes(x = .tau,
                   y = prior.val.scaled,
                   color = prior_name ) ) +
  
  geom_line(linewidth = 1.1) +
  
  geom_point( data = max.points, 
              aes(x = .tau,
                  y = prior.val.scaled,
                  color = prior_name ),
              size = 3) +
  
  xlab( bquote(tau) ) +
  ylab( bquote(p(tau)) ) +
  
  scale_x_continuous(breaks = seq( min(dp2$.tau), max(dp2$.tau), 0.1),
                     limits = c( min(dp2$.tau), max(dp2$.tau) ) ) +
  
  # scale_y_continuous(breaks = seq( min(dp$log.prior), max(dp$log.prior), 0.25),
  #                    limits = c( min(dp$.tau), max(dp$.tau) ) ) +
  
  scale_color_manual(values = my.colors, 
                     name = "") +
  
  theme_bw(base_size = 20) +
  #ggtitle( paste("k = ", .k) ) +
  
  theme(text = element_text(face = "bold"),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )

plot


my_ggsave(name = "zito_all_cause_death_priors.pdf",
          .plot = plot,
          .width = 10,
          .height = 8,
          .results.dir = results.dir,
          .overleaf.dir = overleaf.dir.figs)


# ~ bayesmeta sanity check  -------------------------------------------------

# sanity check: use bayesmeta
# prior is the dashed line
m1 = mybayesmeta(y = .dat$yi,
                 sigma = .dat$sei,
                 tau.prior = "Jeffreys",
                 interval.type = "central")

#plot(m1, prior = TRUE)


m2 = mybayesmeta(y = .dat$yi,
                 sigma = .dat$sei,
                 tau.prior = "Jeffreys2",
                 interval.type = "central")

#plot(m2, prior = TRUE)


bayesmeta_prior = function(prior_name, tau){
  if ( prior_name == "Jeffreys1" ) return( m1$tau.prior(tau) )
  if ( prior_name == "Jeffreys2" ) return( m2$tau.prior(tau) )
}


dp2 = expand_grid( .tau = tau_vec,
                   prior_name = c("Jeffreys1", "Jeffreys2" ) )

dp2 = dp2 %>% rowwise() %>%
  mutate( prior.val = bayesmeta_prior( tau = .tau,
                                       prior_name = prior_name ) )


### Rescale the priors to have the same max height
dp2$prior.val.scaled = NA
dp2$prior.val.scaled[ dp2$prior_name == "Jeffreys1" ] = dp2$prior.val[ dp2$prior_name == "Jeffreys1" ] / max(dp2$prior.val[ dp2$prior_name == "Jeffreys1" ])

dp2$prior.val.scaled[ dp2$prior_name == "Jeffreys2" ] = dp2$prior.val[ dp2$prior_name == "Jeffreys2" ] / max(dp2$prior.val[ dp2$prior_name == "Jeffreys2" ])


### Points that maximize each curve
dp2 = dp2 %>% group_by(prior_name) %>%
  mutate( is.max = ifelse( prior.val.scaled == max(prior.val.scaled), TRUE, FALSE ) )
( max.points = dp2 %>% filter(is.max == TRUE) )


### Make plot

# in same order as N.pretty
my.colors = c("#E075DB", "#F2340E")

plot = ggplot( data = dp2, 
               aes(x = .tau,
                   y = prior.val.scaled,
                   color = prior_name ) ) +
  
  geom_line(linewidth = 1.1) +
  
  geom_point( data = max.points, 
              aes(x = .tau,
                  y = prior.val.scaled,
                  color = prior_name ),
              size = 3) +
  
  xlab( bquote(tau) ) +
  ylab( bquote(p(tau)) ) +
  
  geom_vline( xintercept = 0, lty = 2 ) +
  
  scale_x_continuous(breaks = seq( min(dp2$.tau), max(dp2$.tau), 0.1),
                     limits = c( min(dp2$.tau), max(dp2$.tau) ) ) +
  
  # scale_y_continuous(breaks = seq( min(dp$log.prior), max(dp$log.prior), 0.25),
  #                    limits = c( min(dp$.tau), max(dp$.tau) ) ) +
  
  scale_color_manual(values = my.colors, 
                     name = "") +
  
  theme_bw(base_size = 16) +
  #ggtitle( paste("k = ", .k) ) +
  
  theme(text = element_text(face = "bold"),
        axis.title = element_text(size=20),
        legend.position = "bottom" )

plot

# looks identical; yay!



# FOREST PLOT ----------------------------------------------------

# retrieve existing analysis results
if ( rerun.analyses == FALSE ) {
  setwd(results.dir)
  rsp = fread("zito_forest_results.csv")
}


# set y-axis order
unique(rsp$group)
correct.order = c( "CV death and myocardial infarction (k = 2)",
                   "All-cause death (k = 3)",
                   "Myocardial infarction (k = 2)",
                   "Index ICA (k = 4)",
                   "Index revascularization (k = 4)",
                   "Downstream testing (k = 4)")

rsp$group = factor( rsp$group, levels = rev(correct.order) )
levels(rsp$group)

# reorder methods
correct.order = c("DL-HKSJ", "REML-HKSJ", "Exact", "Jeffreys1-shortest", "Jeffreys2-shortest")
rsp$method = factor(rsp$method.pretty, levels = correct.order)
levels(rsp$method)

# same colors as in analyze_sims_helper.R for prettiness

.colors = c("#246105",
                 "black",
                 "#CC9808",
                 "#E075DB",
                 "#F2340E")



# Mhat forest plot -------------------------------------------------

# find good x-axis limits
summary( exp(rsp$MLo) )
summary( exp(rsp$MHi) )
xmin = 0.05
xmax = 0.85


p = ggplot( data = rsp,
            aes( y = group,
                 x = exp(Mhat), 
                 xmin = exp(MLo), 
                 xmax = exp(MHi),
                 color = method) ) +
  
  # reference line at null
  geom_vline(xintercept = 1,
             lwd = .8,
             color = "gray") +
  
  
  geom_errorbarh( aes(xmax = exp(MHi),
                      xmin = exp(MLo)),
                  height = 0,
                  lwd = 0.8,
                  position = position_dodge(width = 0.5) ) +
  
  
  geom_point(size=3,
             position=position_dodge(width = 0.5) ) +
  
  # manually provided colors
  scale_colour_manual(values = .colors,
                      guide = guide_legend(reverse = TRUE)) +
  
  
  scale_y_discrete( name = "Outcome" ) +
  # scale_x_continuous( limits = c(xmin, xmax),
  #                     breaks = seq(xmin, xmax, 0.05) ) +
  
  # scale_x_log10(breaks = seq(0.1, 10, 1)) +
  # coord_cartesian( xlim = c(0.1, 10)) +
  
  coord_cartesian( xlim = c(0, 10)) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  
  xlab( "Pooled risk ratio with 95% CI" ) +
  
  labs(color  = "Method") +
  
  
  theme_bw(base_size = 16) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=16),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" )

p


my_ggsave(name = "zito_forest.pdf",
          .plot = p,
          .width = 13,
          .height = 10,
          .results.dir = results.dir,
          .overleaf.dir = overleaf.dir.figs)


# Shat forest plot -------------------------------------------------


# rename methods
rsp$method.pretty.tau = rsp$method.pretty
rsp$method.pretty.tau[ rsp$method.pretty.tau == "DL-HKSJ" ] = "DL-Qprofile"
rsp$method.pretty.tau[ rsp$method.pretty.tau == "REML-HKSJ" ] = "REML-Qprofile"
table(rsp$method.pretty.tau)

# reorder methods
correct.order = c("DL-Qprofile", "REML-Qprofile", "Jeffreys1-shortest", "Jeffreys2-shortest")
rsp = rsp %>% filter(method.pretty.tau %in% correct.order)
rsp$method.pretty.tau = factor(rsp$method.pretty.tau, levels = correct.order)
levels(rsp$method.pretty.tau)

# same colors as in analyze_sims_helper.R for prettiness

.colors = c("#246105",
            "black",
            #"#CC9808",
            "#E075DB",
            "#F2340E")


# find good x-axis limits
summary( rsp$SLo )
summary( rsp$SHi )
xmin = 0
xmax = 12


p = ggplot( data = rsp,
            aes( y = group,
                 x = Shat, 
                 xmin = SLo, 
                 xmax = SHi,
                 color = method.pretty.tau) ) +
  
  # reference line at null
  geom_vline(xintercept = 1,
             lwd = .8,
             color = "gray") +
  
  
  geom_errorbarh( aes(xmax = SHi,
                      xmin = SLo),
                  height = 0,
                  lwd = 0.8,
                  position = position_dodge(width = 0.5) ) +
  
  
  geom_point(size=3,
             position=position_dodge(width = 0.5) ) +
  
  # manually provided colors
  scale_colour_manual(values = .colors,
                      guide = guide_legend(reverse = TRUE)) +
  
  
  scale_y_discrete( name = "Outcome" ) +
  # scale_x_continuous( limits = c(xmin, xmax),
  #                     breaks = seq(xmin, xmax, 0.05) ) +
  
  # scale_x_log10(breaks = seq(0.1, 10, 1)) +
  # coord_cartesian( xlim = c(0.1, 10)) +
  
  coord_cartesian( xlim = c(0, 10)) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  
  xlab( bquote( bold( hat(tau) ~ "with 95% CI (log-RR scale)") ) ) +
  
  labs(color  = "Method") +
  
  
  theme_bw(base_size = 16) +
  
  theme( text = element_text(face = "bold"),
         axis.title = element_text(size=16),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = "bottom" )

p


my_ggsave(name = "zito_forest_tau.pdf",
          .plot = p,
          .width = 13,
          .height = 10,
          .results.dir = results.dir,
          .overleaf.dir = overleaf.dir.figs)


# ONE-OFF STATS FOR PAPER  -------------------------------------------------


# retrieve existing analysis results
if ( rerun.analyses == FALSE ) {
  setwd(results.dir)
  rsp = fread("zito_forest_results.csv")
}

# ~ CI width comparisons  -------------------------------------------------

update_result_csv( name = "Mean perc narrower Jeffreys2 vs winning other method",
                   value = round( 100 * ( mean(rsp$CI_ratio - 1) ) ),
                   print = TRUE )


update_result_csv( name = "Mean perc narrower Jeffreys2 vs winning other method k=2",
                   value = round( 100 * ( mean(rsp$CI_ratio[rsp$k==2] - 1 ) ) ),
                   print = TRUE )

# comparing specific CI limits
exp(rsp$MLo[ rsp$group == "CV death and myocardial infarction (k = 2)" & rsp$method.pretty == "Jeffreys2-shortest"])
exp(rsp$MHi[ rsp$group == "CV death and myocardial infarction (k = 2)" & rsp$method.pretty == "Jeffreys2-shortest"])

exp(rsp$MLo[ rsp$group == "CV death and myocardial infarction (k = 2)" & rsp$method.pretty == "REML-HKSJ"])
exp(rsp$MHi[ rsp$group == "CV death and myocardial infarction (k = 2)" & rsp$method.pretty == "REML-HKSJ"])
