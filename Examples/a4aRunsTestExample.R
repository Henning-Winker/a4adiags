#----------------------------------
# Residual analysis with Runs Tests
# written by Henning Winker
# henning.winker@ec.europa.eu
#-----------------------------------

#************************************
# Install library(a4adiags)
# 
# require(devtools)
# devtools::install_github("henning-winker/a4adiags")
#*******************************************************

# Load packages
library(FLCore)
library(FLa4a)
library(a4adiags)

# Load Data
data(ple4)
data(ple4.index)
stk = ple4
idx = FLIndices(ple4.index)
fit <- sca(stock=stk, indices=idx)

# Do Residual on surveys 
obs = as.data.frame(idx) # observed values
hat = as.data.frame(index(fit)) # fitted values

# For combined survey (N)
a4apar(mfrow=c(1,1))
a4aPlotRunstest(obs,hat,add=T)

# By age (1-6)
a4apar(mfrow=c(3,3),plot.cex = 0.7,labs=F)
a4aPlotRunstest(obs,hat,indexselect = 1:9,subplots = "age",add=T)
mtext("Year",side=1,outer=T,line=0.5)
mtext("Survey Residuals",side=2,outer=T,line=0.5)

# check catch at age
obs = as.data.frame(catch.n(stk))
hat = as.data.frame(catch.n(fit))
a4apar(mfrow=c(3,3),plot.cex = 0.7,labs=F)
a4aPlotRunstest(obs,hat,indexselect = 1:9,subplots = "age",add=T)
mtext("Year",side=1,outer=T,line=0.5)
mtext("Catch Residuals",side=2,outer=T,line=0.5)

#Combined
a4apar(mfrow=c(1,1))
a4aPlotRunstest(obs,hat,add=T)



