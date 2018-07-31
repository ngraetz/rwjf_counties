library(data.table)
library(ggplot2)
repo <- 'C:/Users/ngraetz/Documents/repos/rwjf_counties/'

## Load all covariates and save as one big dataset long on county/year, wide on covariate.
input_dir <- paste0(repo, 'nmx_clean_data/')
cov_dir <- paste0(repo, 'covariate_clean_data/')

## Use all county/year/race/sex populations as template and merge everything else.
pop_white <- readRDS(paste0(input_dir, 'nhw_total_pop.RDS'))
pop_white_25_64 <- readRDS(paste0(input_dir, 'nhw_25_64_pop.RDS'))
setnames(pop_white_25_64, 'total_pop', '25_64_pop')
white_working_pop <- merge(pop_white, pop_white_25_64, by=c('fips','year','sex','race'))
white_working_pop[, perc_25_64 := `25_64_pop` / total_pop]
white_working_pop <- white_working_pop[, c('fips','year','sex','race','perc_25_64')]

pop_black <- readRDS(paste0(input_dir, 'nhb_total_pop.RDS'))
pop_black_25_64 <- readRDS(paste0(input_dir, 'nhb_25_64_pop.RDS'))
setnames(pop_black_25_64, 'total_pop', '25_64_pop')
black_working_pop <- merge(pop_black, pop_black_25_64, by=c('fips','year','sex','race'))
black_working_pop[, perc_25_64 := `25_64_pop` / total_pop]
black_working_pop <- black_working_pop[, c('fips','year','sex','race','perc_25_64')]

pop <- rbind(white_working_pop, black_working_pop)

## Merge all other covariates
for(c in c('bea_covs','bls_laus_covs','factfinder_edu','factfinder_fb','saipe_pov')) {
message(c)
cov <- fread(paste0(cov_dir,c,'.csv'))
cov[, year := as.numeric(year)]
merge_vars <- c('fips','year','sex','race')
if(!('race' %in% names(cov))) merge_vars <- merge_vars[merge_vars!='race']
if(!('sex' %in% names(cov))) merge_vars <- merge_vars[merge_vars!='sex']
## Handle Virginia FIPS in BEA
if(c=='bea_covs') {
  map <- fread(paste0(repo, 'covariate_code/bea_county_template.csv'))
  map[, fips := as.character(fips)]
  map[, bea_fips := as.character(bea_fips)]
  map[nchar(fips)==4, fips := paste0('0',fips)]
  map[is.na(bea_fips), bea_fips := fips]
  pop <- merge(pop, map, all.x=TRUE, by='fips')
  pop[, bea_fips := as.character(bea_fips)]
  pop[nchar(bea_fips)==4, bea_fips := paste0('0',bea_fips)]
  cov[, bea_fips := fips]
  cov[, fips := NULL]
  merge_vars <- c('bea_fips','year')
}
message(paste0('Merge using ', paste(merge_vars, collapse = ', ')))
pop <- merge(pop, cov, all.x=TRUE, by=merge_vars)
message(dim(pop)[1])
}

## Save.
saveRDS(pop, paste0(repo, 'covariate_clean_data/combined_covs.RDS'))



