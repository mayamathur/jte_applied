
library(MetaUtility)
library(metafor)
library(meta)

# studies' RR estimates and CI upper bounds
# from Supplementary Figure 3
yi = c(0.66, 0.29)
hi = c(0.94, 1.41)

dat = scrape_meta(type = "RR", est = yi, hi = hi)


### DL
# with KNHA
m1 = rma.uni(yi = yi, 
             vi = vyi, 
             data = dat,
             method = "DL", 
             knha = TRUE)
exp(m1$b)
exp(m1$ci.ub)  # dramatically higher than reported upper limit of 0.90

# without KNHA
m2 = rma.uni(yi = yi, 
             vi = vyi, 
             data = dat,
             method = "DL", 
             knha = FALSE)
exp(m2$b)
exp(m2$ci.ub)  # agrees with the reported 0.90



### REML - results almost the same as DL
m1 = rma.uni(yi = yi, 
             vi = vyi, 
             data = dat,
             method = "REML", 
             knha = TRUE)
exp(m1$b)
exp(m1$ci.ub)  # dramatically higher than reported upper limit of 0.90

# without KNHA
m2 = rma.uni(yi = yi, 
             vi = vyi, 
             data = dat,
             method = "REML", 
             knha = FALSE)
exp(m2$b)
exp(m2$ci.ub)  # agrees with the reported 0.90
