
# Manually enters data from CBT & insomnia meta-analysis.

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


# # ENTER DATA - TST ONLY -------------------------------------------------
# 
# # Annals meta CBT insomnia, page 8
# # total sleep time outcome
# 
# yi = c(2.04, 10.84, -2.00, 33.70, 19.10, -4.20, 56.40, -9.67, 13.20, -12.32, 5.37, 13.50, -6.47, 10.20, 16.20, 11.78,
#   -5.23, 47.90, 62.00, 18.60,
#   65.11, 63.70, 25.80, 6.50)
# 
# hi = c(51.67, 45.63, 24.87, 72.38, 63.47, 26.60, 94.10, 13.49, 34.02, 26.43, 35.23, 58.60, 15.49, 37.77, 36.06, 47.65,
#   35.93, 85.75, 98.94, 39.65,
#   117.28, 98.98, 48.28, 49.21)
# 
# group = c( rep("Posttreatment", 16),
#              rep("Early follow-up", 4),
#              rep("Late follow-up", 4) )
# 
# d = scrape_meta(type = "raw", est = yi, hi = hi)
# d$group = group
# names(d)[names(d) == "vyi"] = "vi"
# d$sei = sqrt(d$vi)
# 
# 
# # sanity checks: reproduce their results
# rma.uni(yi = yi, vi = vi, data = d %>% filter(group=="Posttreatment"), knha = TRUE)
# rma.uni(yi = yi, vi = vi, data = d %>% filter(group=="Early follow-up"), knha = TRUE)
# rma.uni(yi = yi, vi = vi, data = d %>% filter(group=="Late follow-up"), knha = TRUE)


# ENTER DATA FOR EACH OUTCOME -------------------------------------------------

# each outcome is a separate forest plot; Figs 1-4

# ~ SOL (Fig 1) -------------------------------------------------
yi = c(-21.29, -29.00, -12.60, 2.85, -26.00, -22.30, -10.90, -25.00, -24.65, -27.70, -13.84, -9.50, -21.83, -5.30, -9.00, -29.50,
       -15.40, -6.80, -44.60, -2.00,
       -39.70, -15.50, 0.90)

hi = c(0.77, -18.67, 6.40, 24.97, 2.13, 12.29, 3.53, -10.65, -4.85, -22.37, 3.40, 10.26, -11.71, 24.09, 0.35, -3.98,
       3.43, 19.72, -10.30, 7.35,
       -4.90, 0.10, 21.05)


group = c( rep("SOL; posttreatment (k = 16)", 16),
           rep("SOL; early follow-up (k = 4)", 4),
           rep("SOL; late follow-up (k = 3)", 3) )



d1 = scrape_meta(type = "raw", est = yi, hi = hi)
d1$group = group

# sanity check vs. Figure 1
# agrees exactly :)
for ( .g in unique(d1$group) ) {
  cat("\n\n        ******* GROUP", .g)
  cat("\n\n")
  
  .dat = d1 %>% filter(group == .g)
  
  print( rma( yi = .dat$yi, vi = .dat$vyi,
              method = "DL", knha = TRUE ) )
}




# ~ TST (Fig 3)  -------------------------------------------------
yi = c(2.04, 10.84, -2.00, 33.70, 19.10, -4.20, 56.40, -9.67, 13.20, -12.32, 5.37, 13.50, -6.47, 10.20, 16.20, 11.78,
  -5.23, 47.90, 62.00, 18.60,
  65.11, 63.70, 25.80, 6.50)

hi = c(51.67, 45.63, 24.87, 72.38, 63.47, 26.60, 94.10, 13.49, 34.02, 26.43, 35.23, 58.60, 15.49, 37.77, 36.06, 47.65,
  35.93, 85.75, 98.94, 39.65,
  117.28, 98.98, 48.28, 49.21)


group = c( rep("TST; posttreatment (k = 16)", 16),
             rep("TST; early follow-up (k = 4)", 4),
             rep("TST; late follow-up (k = 4)", 4) )

d2 = scrape_meta(type = "raw", est = yi, hi = hi)
d2$group = group


# sanity check vs. Figure 3
# agrees exactly :)
for ( .g in unique(d2$group) ) {
  cat("\n\n        ******* GROUP", .g)
  cat("\n\n")
  
  .dat = d2 %>% filter(group == .g)
  
  print( rma( yi = .dat$yi, vi = .dat$vyi,
              method = "DL", knha = TRUE ) )
}



# ~ SE% (Fig 4) -------------------------------------------------
yi = c(10.95, 12.16, 6.00, 14.90, 10.90, 5.90, 15.20, 5.06, 9.90, 10.36, 13.80, 6.97, 4.80, 7.54, 7.40, 13.51, 13.05,
       10.39, 15.90, -2.70, 24.40, 8.29,
       
       18.81, 17.70, 9.70, 1.30 )

hi = c(21.78, 19.87, 11.47, 21.65, 21.54, 13.62, 24.17, 10.09, 14.58, 16.03, 15.77, 13.47, 12.28, 11.77, 12.93, 18.23, 20.12,
       
       19.10, 22.40, 8.39, 33.37, 13.31,
       30.35, 26.25, 14.93, 8.22)

group = c( rep("SE%; posttreatment (k = 17)", 17),
           rep("SE%; early follow-up (k = 5)", 5),
           rep("SE%; late follow-up (k = 4)", 4) )



d3 = scrape_meta(type = "raw", est = yi, hi = hi)
d3$group = group

# sanity check vs. Figure 4
# agrees exactly :)
for ( .g in unique(d3$group) ) {
  cat("\n\n        ******* GROUP", .g)
  cat("\n\n")
  
  .dat = d3 %>% filter(group == .g)
  
  print( rma( yi = .dat$yi, vi = .dat$vyi,
              method = "DL", knha = TRUE ) )
}



# ~ WASO (Fig 2) -------------------------------------------------
yi = c(-24.36, -16.78, -31.70, -73.60, -11.27, -23.41, -27.40, -19.07, -24.12, -9.10, -20.63, -15.00, -75.60, -20.96,
       -20.11, -70.20, -41.40,
       -29.33, -26.70, -3.30)

hi = c(9.78, 9.82, -6.72, -34.42, 10.69, -0.42, -5.97, -1.42, 11.91, 17.78, -4.95, 20.05, -48.93, -6.76,
       10.13, -31.32, -13.84,
       5.65, -0.76, 23.30)


group = c( rep("WASO; posttreatment (k = 14)", 14),
           rep("WASO; early follow-up (k = 3)", 3),
           rep("WASO; late follow-up (k = 3)", 3) )



d4 = scrape_meta(type = "raw", est = yi, hi = hi)
d4$group = group

# sanity check vs. Figure 2
# agrees exactly :)
for ( .g in unique(d4$group) ) {
  cat("\n\n        ******* GROUP", .g)
  cat("\n\n")
  
  .dat = d4 %>% filter(group == .g)
  
  print( rma( yi = .dat$yi, vi = .dat$vyi,
       method = "DL", knha = TRUE ) )
}



# combine them
d = bind_rows(d1, d2, d3, d4)
names(d)[names(d) == "vyi"] = "vi"
d$sei = sqrt(d$vi)





# SAVE DATA  -------------------------------------------------

setwd(here())
fwrite(d, "data_insomnia.csv")



