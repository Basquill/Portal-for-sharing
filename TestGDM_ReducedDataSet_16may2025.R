
# Test GDM matrix permutation script with reduced dataset

library(tidyverse)
library(gdm)

sessionInfo()

#R version 4.4.2 (2024-10-31 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 10 x64 (build 19045)

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] gdm_1.6.0-5     lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
#[7] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0

## Set up project directory ------------------------------------------------------------------------------------------------

wd <- ""
setwd(wd)

## Load data
# 2258 plots
SXY <-read.csv("./SXY_abioOnly_proxOmit_22april2025.csv",stringsAsFactors = TRUE, header = TRUE)

# Create subset of 500 plots for quick run

SXY <- sample_n(SXY, 500)

# Split site info (S); predictors (X); response variables (Y) from SXY
S <- subset(SXY, select = c(1:3)) # Plot ID plus geo-coordinates
X <- subset(SXY, select = c(9:21)) # all predictors
Y <- subset(SXY, select = c(27, 31:33, 40)) # subset of response variables

# apply cube transformation on response data -----------------------------------------------------------------------------------------

Ycube <- Y^(1/3)

# Assemble data objects needed for GDM  ----------------------------------------------------------------------------------------------

sppData <- cbind(S,Ycube)
envData <- cbind(S,X)

# assemble site pair table required for GDM

gdmTab <- formatsitepair(bioData=sppData, bioFormat =1, XColumn="EASTING", YColumn="NORTHING", abundance = TRUE, sppColumn = "species",
                         siteColumn="Plot", predData=envData)

# test fit

gdm <- gdm(data=gdmTab,geo=TRUE) # 
str(gdm) # dev explained = 31.1

##########################################################################################################
# Calculate significance of each predictor using Monte Carlo permutation with stepwise backward selection
###########################################################################################################

gc()

modTest <- gdm.varImp(gdmTab, geo=T, nPerm=50, parallel=T, cores=2, predSelect = T,
                      outFile = "modtest_AbioOnly_15may2025")

# get "Error in nullDev - varDevTab : non-numeric argument to binary operator" 
# error occurs after initiation of "Fitting GDMs to the permuted site-pair tables..."

