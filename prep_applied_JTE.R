
# Manually enters data from Zito's CAD meta-analysis.

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


# ENTER DATA FOR EACH OUTCOME -------------------------------------------------

# each outcome is a separate forest plot; Supplementary Figs 3-6
# specifically, the CCTA vs. SPECT-MPI studies

# note that Supp Fig 5 is omitted because there was only 1 study

# these are on RR scale still
yi = c(0.66, 0.29,  # CV & MI (Fig 3)
       0.32, 0.33, 0.80, # all-cause death (Fig 4)
       0.71, 0.29,  # MI (Fig 6)
       0.74, 1.16, 1.68, 1.24,  # index ICA (Fig 7)
       1.12, 0.79, 6.85, 1.68,  # index revascularization (Fig 8)
       0.71, 1.95, 0.92, 1.08 )  # downstream testing (Fig 9)
      
hi = c(0.94, 1.41,
       3.07, 8.06, 1.15,
       1.22, 1.41, 
       0.99, 2.09, 4.06, 1.42,
       1.69, 1.96, 54.52, 2.08,
       0.94, 2.95, 1.74, 1.18 )


group = c( rep("CV death and myocardial infarction (k = 2)", 2),
           rep("All-cause death (k = 3)", 3),
           rep("Myocardial infarction (k = 2)", 2),
           rep("Index ICA (k = 4)", 4),
           rep("Index revascularization (k = 4)", 4),
           rep("Downstream testing (k = 4)", 4)
           )



d = scrape_meta(type = "RR", est = yi, hi = hi)
d$group = group


# sanity check vs. their results
for ( .g in unique(d$group) ) {
  cat("\n\n        ******* GROUP", .g)

  
  .dat = d %>% filter(group == .g)
  
  # note: CIs differ a LOT with KNHA vs. without; their analyses clearly don't use it
  m = rma( yi = .dat$yi, vi = .dat$vyi,
           method = "DL", knha = FALSE )
  
  cat("\n Without KNHA: RR = ", round( exp(m$beta), 3), ", hi = ", round( exp(m$ci.ub), 3) )
  
  m = rma( yi = .dat$yi, vi = .dat$vyi,
           method = "DL", knha = TRUE )
  
  cat("\n With KNHA: RR = ", round( exp(m$beta), 3), ", hi = ", round( exp(m$ci.ub), 3) )
  
  cat("\n\n")
}




# rename
names(d)[names(d) == "vyi"] = "vi"
d$sei = sqrt(d$vi)





# SAVE DATA  -------------------------------------------------

setwd(here())
fwrite(d, "data_zito.csv")



