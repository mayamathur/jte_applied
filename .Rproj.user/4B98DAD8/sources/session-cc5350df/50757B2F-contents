

# Annals meta CBT insomina, page 8
# total sleep time

yi = c(2.04, 10.84, -2.00, 33.70, 19.10, -4.20, 56.40, -9.67, 13.20, -12.32, 5.37, 13.50, -6.47, 10.20, 16.20, 11.78,
  -5.23, 47.90, 62.00, 18.60,
  65.11, 63.70, 25.80, 6.50)


lo = c(-47.59, -23.95, -28.87, -4.98, -25.27, -35.00, -18.70, -32.83, -7.62, -51.07, -24.49, -31.60, -28.43, -17.37, -3.66, -24.09,
  -46.39, 10.05, 25.06, -2.45,
  12.94, 28.42, 3.32, -26.21)

hi = c(51.67, 45.63, 24.87, 72.38, 63.47, 26.60, 94.10, 13.49, 34.02, 26.43, 35.23, 58.60, 15.49, 37.77, 36.06, 47.65,
  35.93, 85.75, 98.94, 39.65,
  117.28, 98.98, 48.28, 49.21)

group = c( rep("Posttreatment", 16),
             rep("Early follow-up", 4),
             rep("Late follow-up", 4) )

d = scrape_meta(type = "raw", est = yi, hi = hi)
d$group = group
names(d)[names(d) == "vyi"] = "vi"
d$sei = sqrt(d$vi)


# sanity checks: reproduce their results
rma.uni(yi = yi, vi = vi, data = d %>% filter(group=="Posttreatment"), knha = TRUE)
rma.uni(yi = yi, vi = vi, data = d %>% filter(group=="Early follow-up"), knha = TRUE)
rma.uni(yi = yi, vi = vi, data = d %>% filter(group=="Late follow-up"), knha = TRUE)
