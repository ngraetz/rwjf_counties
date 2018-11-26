## Set location of repo.
repo <- 'C:/Users/ngraetz/Documents/repos/rwjf_counties/'

## Load all libraries and functions.
library(plyr)
library(INLA)
library(data.table)
library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(Hmisc)
library(spdep)
source(paste0(repo, 'functions.R'))
source(paste0(repo, 'functions_shapley.R'))
input_dir <- paste0(repo, 'nmx_clean_data/')
cov_dir <- paste0(repo, 'covariate_clean_data/')
#out_dir <- 'C:/Users/ngraetz/Dropbox/Penn/papers/rwjf/paa_materials/'
out_dir <- paste0(repo, '/results')

## Set options for this run (data domain and covariates).
## Current datasets: 25-64 ASDR for national NHW male, national NHW female, South NHW male, South NWH female, South NHB male, South NHB female.
## Current possible covariates:
##    BEA:    "percent_transfers","percent_wage_salary_employment","income_per_capita","total_employees"
##    BLS:    "labor_force","employed","unemployed","percent_unemployment"
##    SAIPE:  "poverty_all","poverty_0_17","poverty_5_17","median_hh_income"
##    FactFinder (Census + ACS): "fb","less_12","hs","assoc","college"
##    Pops:   "perc_25_64"
##    AHRF:   "mds_pc"
#race <- 'nhw'
sex_option <- 1
domain <- 'national'
#cov_domain <- 'all'
for(cov_domain in c('ses','med','beh','all')) {
  if(cov_domain=='ses') covs <- c('college','log_hh_income','percent_transfers') ## To add: eviction_rate, perc_manufacturing
  if(cov_domain=='med') covs <- c('log_mds_pc','chr_mammography','chr_diabetes_monitoring') ## To add: insurance (SAHIE)
  if(cov_domain=='beh') covs <- c('current_smoker_prev','obesity_prev','as_heavy_drinking_prev')
  if(cov_domain=='pop') covs <- c('fb','perc_25_64') ## To add: perc_black, perc_hispanic, perc_native, net_migration
  if(cov_domain=='all') covs <- c('college','poverty_all','log_hh_income','percent_transfers',
                                  'log_mds_pc','chr_mammography','chr_diabetes_monitoring',
                                  'as_diabetes_prev','current_smoker_prev','obesity_prev','as_heavy_drinking_prev',
                                  'fb') ## To add: perc_black, perc_hispanic, perc_native, net_migration
  if(cov_domain=='custom') covs <- c('college','poverty_all','log_hh_income','percent_transfers','percent_unemployment','fb','perc_25_64')
  year_range <- c(2000,2010,2015)
  plot_trends <- FALSE

  start_year <- min(year_range)
  end_year <- max(year_range)
  sex_name <- ifelse(sex_option==1,'Male','Female')
  #output_name <- paste0(capitalize(domain), ' ', toupper(race), ' ', capitalize(sex_name))

  ## Load master shapefile.
  counties <- readOGR(paste0(repo, "/cb_2016_us_county_20m"), 'cb_2016_us_county_20m')
  counties@data <- transform(counties@data, fips = paste0(STATEFP, COUNTYFP)) # create unique county 5-digit fips
  counties <- counties[counties@data$fips != "02016", ] # Drop Aleutians West, AK - screws up plots
  counties <- counties[counties$STATEFP != '02' &
                         counties$STATEFP != '15' &
                         counties$STATEFP != '72',]
  background_counties <- counties[counties$STATEFP != '02' &
                                    counties$STATEFP != '15' &
                                    counties$STATEFP != '72',]
  background.dt <- as.data.table(fortify(background_counties, region = 'STATEFP'))
  map_colors <- c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')
  counties$state <- as.numeric(counties$STATEFP)
  states <- gUnaryUnion(counties, id = counties@data$state)
  background.dt.states <- as.data.table(fortify(states, region = 'state'))

  ## Load mortality data, subset to race/domain, merge all contextual covariates, and merge total population to create age-standardized deaths for binomial model.
  #mort <- readRDS(paste0(input_dir, 'asdr_', race, '_', domain, '.RDS'))
  mort <- readRDS(paste0(input_dir, 'age_specific_mort_cod.RDS'))
  mort[, fips := as.character(fips)]
  mort[nchar(fips)==4, fips := paste0('0',fips)]
  setnames(mort, c('total_deaths','pooled_year'), c('deaths','year'))
  pop_nhw <- readRDS(paste0(input_dir, 'nhw', '_total_pop.RDS'))
  pop_nhb <- readRDS(paste0(input_dir, 'nhb', '_total_pop.RDS'))
  pop <- rbind(pop_nhw,pop_nhb)
  #pop <- readRDS(paste0(input_dir, race, '_total_pop.RDS'))
  all_covs <- readRDS(paste0(cov_dir, 'combined_covs.RDS'))
  all_covs[, log_mds_pc := log(mds_pc+0.01)]

  total_race <- mort[, lapply(.SD, sum, na.rm=TRUE), by=c('year','sex','fips'), .SDcols=c('total_pop','black','hisp')]
  total_race[, perc_black := black / total_pop]
  total_race[, perc_hispanic := hisp / total_pop]
  all_covs <- merge(all_covs, total_race[, c('year','sex','fips','perc_black','perc_hispanic')], by=c('year','sex','fips'))

  all_covs <- as.data.table(all_covs %>% mutate_at(funs(scale(.) %>% as.vector), .vars=covs))

  ## Old method by race
  # mort <- merge(mort, pop, by=c('fips','year','sex','race'))
  # mort <- merge(mort, all_covs, by=c('fips','year','sex','race'))
  # race_option <- ifelse(race=='nhw',0,1)
  # mort <- mort[sex==sex_option & race==race_option, ]
  # mort[, deaths := round(nmx * total_pop)]
  ## Now with all-race, collapse covs using pops
  pop <- merge(pop, all_covs, by=c('fips','year','sex','race'))
  pop <- pop[year %in% year_range, lapply(.SD, weighted.mean, w=total_pop, na.rm=TRUE), by=c('fips','year','sex'), .SDcols=covs]
  mort <- merge(mort, pop, by=c('fips','year','sex'))

  mort[is.na(total_pop) | total_pop == 0, total_pop := 1]
  mort <- mort[year %in% year_range, ]
  mort <- mort[sex == sex_option, ]

  ## Drop Alaska, Hawaii, Puerto Rico
  mort <- mort[!grep(paste(c('^02','^15','^72'),collapse="|"), fips, value=TRUE), ]

  ## Merge region and metro codes.
  metro_codes <- fread(paste0(repo, 'covariate_clean_data/FIPSmetroregion.csv'))
  metro_codes[, fips := as.character(fips)]
  metro_codes[nchar(fips)==4, fips := paste0('0',fips)]
  mort <- merge(mort, metro_codes[, c('fips','metroname','regionname','metro','region')], by=c('fips'), all.x=TRUE)
  mort[metroname %in% c('Nonmetro, adjacent', 'Nonmetro, nonadjacent'), metroname := 'Nonmetro']
  mort[, metro_region := paste0(metroname, '_', regionname)]
  mort[, metro_region := factor(metro_region, levels=c('Lg central metro_Pacific',unique(mort[metro_region!='Lg central metro_Pacific', metro_region])))]

  ## subset counties in the South, according to specified domain for this model run.
  ## https://www2.census.gov/geo/pdfs/maps-data/maps/reference/us_regdiv.pdf
  ##    Texas, Oklahoma, Arkansas, Louisianna, Mississippi
  ##    Alabama, Georgia, Tennessee, Kentucky, West Virginia, Florida
  ##    South Carolina, North Carolina, Virginia, DC, Maryland, Delaware
  south_states <- c('TX','OK','AR','LA','MS','AL','TN','KY','WV','VA','MD','DE','NC','SC','GA','FL')
  south_fips <- metro_codes[state %in% south_states, fips]
  if(domain=='south') {
    message('Subsetting to Southern FIPS codes')
    mort <- mort[fips %in% south_fips, ]
  }

  ## Create a neighborhood file using Queens convention.
  message('Creating spatial weight matrix...')
  nb.FOQ <- poly2nb(counties, queen=TRUE)
  ## Create an INLA weight matrix.
  lw.FOQ <- nb2INLA("FOQ_INLA",nb.FOQ)
  an <- data.frame(1:length(counties),counties$fips)
  o <- match(mort$fips,an$counties.fips)
  ID <- an$X1.length.counties.[o]
  mort[, ID := ID]

  for(c in covs) {
    mort[is.nan(get(c)), (c) := NA]
    mort[is.na(get(c)), drop := 1]
  }
  mort <- mort[is.na(drop),]

  ## Fit INLA model, save coefficients for table.
  ## Make tables comparing coefficients
  ## 1. Metro * Regions, AR1 on year
  ## 2. Metro * Regions + Covariates, AR1 on year
  ## 3. Metro * Regions + Covariates, AR1 on year + Besag on county
  message('Fitting INLA models...')
  mort[, regionname := factor(regionname, levels = c('Pacific', unique(mort[!(regionname %in% 'Pacific'), regionname])))]
  mort[, metroname := factor(metroname, levels = c('Lg central metro', unique(mort[!(metroname %in% 'Lg central metro'), metroname])))]
  ## Create broad age groups for interacting.
  mort[agegrp <= 20, broad_age := '0_24']
  mort[agegrp %in% 25:40, broad_age := '25_44']
  mort[agegrp %in% 45:60, broad_age := '45_64']
  mort[agegrp >= 65, broad_age := '65_100']
  broad_ages <- c('0_24','25_44','45_64','65_100')
  inla_formulas <- c(paste0('deaths ~ ', paste(covs, collapse = ' + '), ' + as.factor(agegrp) + year'),
                     paste0('deaths ~ ', paste(covs, collapse = ' + '), ' + as.factor(agegrp) + as.factor(metro_region) + year'),
                     paste0('deaths ~ ', paste(covs, collapse = ' + '), ' + as.factor(agegrp) + as.factor(metro_region) + year + f(ID,model="besag",graph="FOQ_INLA")'))
  full_mort <- copy(mort)

  d <- full_mort[year %in% c(start_year,end_year), ]
  d <- dcast(d, fips + agegrp ~ year, value.var = c(covs, 'total_pop', 'metro_region', 'nmx'))
  d[, (paste0('year_', start_year)) := start_year]
  d[, (paste0('year_', end_year)) := end_year]

  ## Try running all age groups separately and binding together.
  #age_groups <- list(c(0,20), c(25,40), c(45,60), c(65,80))
  age_groups <- list(c(25,40), c(45,60))
  ## Run all permutations of covariate list
  permutations <- list()
  for(c in 1:length(covs)) permutations[[c]] <- c(0,1)
  permutations <- as.data.table(expand.grid(permutations))
  setnames(permutations, names(permutations), covs)
  permutations <- permutations[-1,]
  # Run each covariate permutation
  get_coefs <- function(f) {
    inla_formula <- inla_formulas[f]
    coefs <- rbindlist(lapply(age_groups, model_permute,
                              data=full_mort, inla_f=inla_formula,
                              coef_file=paste0(sex_option,'_',c),
                              shapley=FALSE,
                              shap_covs=covs))
  }
  all_coefs <- rbindlist(lapply(1:3, get_coefs), fill=T)
}
