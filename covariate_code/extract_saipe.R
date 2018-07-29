library(data.table)
library(ggplot2)
library(haven)
library(rgdal)
library(readxl)

## Poverty/income, 2015
pov_2015 <- fread('C:/Users/ngraetz/Downloads/est15all.csv')
pov_2015[, poverty_all := as.numeric(`Poverty Percent, All Ages`)]
pov_2015[, poverty_0_17 := as.numeric(`Poverty Percent, Age 0-17`)]
pov_2015[, poverty_5_17 := as.numeric(`Poverty Percent, Age 5-17 in Families`)]
pov_2015[, median_hh_income := as.numeric(gsub(',','',`Median Household Income`))]
pov_2015[, statefips := as.character(`State FIPS Code`)]
pov_2015[nchar(statefips)==1, statefips := paste0('0',statefips)]
pov_2015[, countyfips := as.character(`County FIPS Code`)]
pov_2015[nchar(countyfips)==1, countyfips := paste0('00',countyfips)]
pov_2015[nchar(countyfips)==2, countyfips := paste0('0',countyfips)]
pov_2015[, fips := paste0(statefips, countyfips)]
pov_2015[, year := 2015]
pov_2015 <- pov_2015[, c('fips','year','poverty_all','poverty_0_17','poverty_5_17','median_hh_income')]

## Poverty/income, 2010
pov_2010 <- fread('C:/Users/ngraetz/Downloads/est10all.csv')
pov_2010[, poverty_all := as.numeric(`Poverty Percent All Ages`)]
pov_2010[, poverty_0_17 := as.numeric(`Poverty Percent Under Age 18`)]
pov_2010[, poverty_5_17 := as.numeric(`Poverty Percent Ages 5-17`)]
pov_2010[, median_hh_income := as.numeric(gsub(',','',`Median Household Income`))]
pov_2010[, statefips := as.character(`State FIPS`)]
pov_2010[, countyfips := as.character(`County FIPS`)]
pov_2010[nchar(countyfips)==1, countyfips := paste0('00',countyfips)]
pov_2010[nchar(countyfips)==2, countyfips := paste0('0',countyfips)]
pov_2010[, fips := paste0(statefips, countyfips)]
pov_2010[, year := 2010]
pov_2010 <- pov_2010[, c('fips','year','poverty_all','poverty_0_17','poverty_5_17','median_hh_income')]

## Poverty/income, 2000
pov_2000 <- read.table('https://www2.census.gov/programs-surveys/saipe/datasets/2000/2000-state-and-county/est00all.dat', fill = TRUE, header = FALSE)
## From Census doc: https://www2.census.gov/programs-surveys/saipe/technical-documentation/file-layouts/state-county/2000-estimate-layout.txt
cols <- c(1,2,6,12,18,21) 
col_names <- c('statefips','countyfips','poverty_all','poverty_0_17','poverty_5_17','median_hh_income')
setnames(pov_2000, names(pov_2000)[cols], col_names)
pov_2000 <- pov_2000[col_names]
pov_2000 <- as.data.table(pov_2000)
## Format FIPS code
pov_2000[, statefips := as.character(statefips)]
pov_2000[, countyfips := as.character(countyfips)]
pov_2000[nchar(countyfips)==1, countyfips := paste0('00',countyfips)]
pov_2000[nchar(countyfips)==2, countyfips := paste0('0',countyfips)]
pov_2000[, fips := paste0(statefips, countyfips)]
pov_2000[, year := 2000]
pov_2000 <- pov_2000[, c('fips','year','poverty_all','poverty_0_17','poverty_5_17','median_hh_income')]

## Poverty/income, 1990
pov_1990 <- read.table('https://www2.census.gov/programs-surveys/saipe/datasets/1989/1989-state-and-county/est89all.dat', fill = TRUE, header = FALSE)
## From Census doc: https://www2.census.gov/programs-surveys/saipe/technical-documentation/file-layouts/state-county/2000-estimate-layout.txt
cols <- c(1,2,6,12,18,21) 
col_names <- c('statefips','countyfips','poverty_all','poverty_0_17','poverty_5_17','median_hh_income')
setnames(pov_1990, names(pov_1990)[cols], col_names)
pov_1990 <- pov_1990[col_names]
pov_1990 <- as.data.table(pov_1990)
## Format FIPS code
pov_1990[, statefips := as.character(statefips)]
pov_1990[, countyfips := as.character(countyfips)]
pov_1990[nchar(countyfips)==1, countyfips := paste0('00',countyfips)]
pov_1990[nchar(countyfips)==2, countyfips := paste0('0',countyfips)]
pov_1990[, fips := paste0(statefips, countyfips)]
pov_1990[, year := 1989]
pov_1990 <- pov_1990[, c('fips','year','poverty_all','poverty_0_17','poverty_5_17','median_hh_income')]

## Append and save
all_saipe_pov <- rbind(pov_1990, pov_2000, pov_2010, pov_2015)
write.csv(all_saipe_pov, 'C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/saipe_pov.csv', row.names=FALSE)


