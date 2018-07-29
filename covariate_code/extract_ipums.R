library(data.table)
library(ggplot2)
library(haven)

## Load IPUMS extract from data portal.
d <- fread("C:/Users/ngraetz/Downloads/usa_00003.csv/usa_00003.csv")
d[, STATEFIP := as.character(STATEFIP)]
d[, COUNTYFIPS := as.character(COUNTYFIPS)]
d[nchar(STATEFIP)==1, STATEFIP := paste0('0',STATEFIP)]
d[nchar(COUNTYFIPS)==2, COUNTYFIPS := paste0('0',COUNTYFIPS)]
d[nchar(COUNTYFIPS)==1, COUNTYFIPS := paste0('00',COUNTYFIPS)]
d[, fips := paste0(STATEFIP, COUNTYFIPS)]

## Recode all categorical variables using codebook from IPUMS.
cb <- fread("C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/cb_ipums.csv")
all_values <- as.numeric(gsub('value_', '', grep('value_', names(cb), value=TRUE)))
all_values <- all_values[!is.na(all_values)]
for(v in unique(cb[, variable])) {
  message(paste0('Recoding ', v, '...'))
  for(val in all_values) {
    value_recode <- cb[variable==v, get(paste0('recode_', val))]
    if(!is.na(value_recode)) {
      ## Handle multiple recodes.
      if(grepl(',', cb[variable==v, get(paste0('value_', val))])) {
        value_set <- as.numeric(unlist(strsplit(cb[variable==v, get(paste0('value_', val))], split=',')))
        d[get(v) %in% value_set, (paste0(v,'_recode')) := value_recode]
      }
      ## Handle continuous recodes.
      if(grepl('-', cb[variable==v, get(paste0('value_', val))])) {
        value_set <- as.numeric(unlist(strsplit(cb[variable==v, get(paste0('value_', val))], split='-')))
        d[get(v) <= max(value_set) & get(v) >= min(value_set), (paste0(v,'_recode')) := value_recode]
      }
      ## Handle single values.
      if(!grepl('-', cb[variable==v, get(paste0('value_', val))]) & !grepl(',', cb[variable==v, get(paste0('value_', val))])) {
        value_set <- as.numeric(cb[variable==v, get(paste0('value_', val))])
        d[get(v) == value_set, (paste0(v,'_recode')) := value_recode]
      }
    }
  }
  ## Handle missing values from IPUMS.
  missing_set <- as.numeric(unlist(strsplit(cb[variable==v, value_missing], split=',')))
  d[get(v) %in% missing_set, (paste0(v,'_recode')) := NA]
  ## If value is still missing but not in official missing set, that means I just haven't codebooked that category so we will code it as other=999 for now.
  ## These people should be included in the denominators of all percents/rates at the county-level.
  d[is.na(get(v)) & !(get(v) %in% missing_set), (paste0(v,'_recode')) := 999]
}

## BPL < 150 = native-born, BPL >= 150 = foreign-born, 999 = missing.
d[BPL < 150, foreign := 0]
d[BPL >= 150, foreign := 1]
d[BPL == 999, foreign := NA]

## EDUC >= 10 = college, EDUC < 10 = no college, EDUC = 0 = missing.
d[EDUC >= 10, college := 1]
d[EDUC < 10, college := 0]
d[EDUC == 0, college := NA]

## HHINCOME = 9999999 = missing.
## Need to adjust for inflation.

## VALUEH = 0000000, 9999998, 9999999 = missing.
## Need to adjust for inflation. 
## IPUMS warns against using this longitudinally because of shifting methods of accounting.

## Wages: INCTOT = INCWAGE + INCWELFR + INCSS + others.
## Make ratio of total wages to welfare within county?
## INCWAGE = 999999, 999998 = missing.
## INCSS = 99999 = missing.
## INCWELFR = 99999 = missing.
d[INCWAGE %in% c(999999, 999998), INCWAGE := NA]
d[INCSS == 99999, INCSS := NA]
d[INCWELFR == 99999, INCSS := NA]

## POVERTY = % of poverty threshold.
d[POVERTY < 100, poverty_percent := 1]
d[POVERTY >= 100, poverty_percent := 0]

## TRANTIME = in minutes.
## TRANTIME = 0 = missing.
d[TRANTIME==0, TRANTIME := NA]

## Make percent indicators for agriculture, mining, manufacturing, construction. 
d[IND1990_recode==1, agriculture_percent := 1]
d[IND1990_recode!=1 & !is.na(IND1990_recode), agriculture_percent := 0]

d[IND1990_recode==2, mining_percent := 1]
d[IND1990_recode!=2 & !is.na(IND1990_recode), mining_percent := 0]

d[IND1990_recode==3, manufacturing_percent := 1]
d[IND1990_recode!=3 & !is.na(IND1990_recode), manufacturing_percent := 0]

d[IND1990_recode==4, construction_percent := 1]
d[IND1990_recode!=4 & !is.na(IND1990_recode), construction_percent := 0]

## Collapse weighting by PERWT. 
d_collapsed <- d[, list(foreign_percent = weighted.mean(foreign, PERWT, na.rm=TRUE),
                        college_percent = weighted.mean(college, PERWT, na.rm=TRUE),
                        poverty_percent = weighted.mean(poverty_percent, PERWT, na.rm=TRUE),
                        agriculture_percent = weighted.mean(agriculture_percent, PERWT, na.rm=TRUE),
                        mining_percent = weighted.mean(mining_percent, PERWT, na.rm=TRUE),
                        manufacturing_percent = weighted.mean(manufacturing_percent, PERWT, na.rm=TRUE),
                        construction_percent = weighted.mean(construction_percent, PERWT, na.rm=TRUE),
                        unemployed_percent = weighted.mean(EMPSTAT_recode, PERWT, na.rm=TRUE),
                        average_commute_minutes = weighted.mean(TRANTIME, PERWT, na.rm=TRUE),
                        total_income_wages = weighted.sum(INCWAGE, PERWT, na.rm=TRUE),
                        total_income_ss = weighted.sum(INCSS, PERWT, na.rm=TRUE),
                        total_income_welfare = weighted.sum(INCWELFR, PERWT, na.rm=TRUE)),
                 by = 'fips']

## Format and save.