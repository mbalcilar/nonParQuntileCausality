# Mehmet Balcilar, 2014-7-4

# Note: The implementation of the nonparametric causality in quantiles is based
#       only one lag. Code needs modification for higher lag orders

library(quantreg) # quantile regression
library(KernSmooth) # for bandwidth selection
library(np)
library(ggplot2)  # for graphing

#rm(list=ls(all=TRUE))

source('lrq_causality_functions_v2.R')

data <- read.csv("gold_oil.csv")

# quantiles at which th statistics is computed
#q0 <-seq(0.1, 0.90, by = 0.10)    # 0.10, 0.20, .... , 0.80, 0.90
q0 <-seq(0.05, 0.95, by = 0.05)    # 0.05, 0.10, .... , 0.90, 0.95

# causality in conditional mean
res.mean <- lrq.causality.test(x=data$Oil, y=data$Gold, type="mean", q=q0)   # note type = "mean"
print(res.mean)

# causality in conditional variance
res.variance <- lrq.causality.test(x=data$Oil, y=data$Gold, type="variance", q=q0)   # note type = "varince"
print(res.variance)

# plot the results

plot.mean <- do.causality.figure(res.mean,title="Causality in Mean")
print(plot.mean)

plot.variance <- do.causality.figure(res.variance,title="Causality in Variance")
print(plot.variance)

## save figures in png files

png(file = "figure_mean.png", height = 14, width = 22, units = "cm", 
    res = 360, bg = "transparent")
print(plot.mean)
dev.off()

png(file = "figure_variance.png", height = 14, width = 22, units = "cm", 
    res = 360, bg = "transparent")
print(plot.variance)
dev.off()


