## Set location of repo.
repo <- '/share/code/geospatial/ngraetz/rwjf_counties/'

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
library(dplyr)
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
sex_option <- 3
domain <- 'national'
mort_domain <- 'allcause'
age_group_option <- c(25,64)
ses_covs <- c('college','poverty_all','percent_transfers','manufacturing') ## Could add back log_hh_income; To add: eviction_rate
med_covs <- c('log_mds_pc','chr_mammography','chr_diabetes_monitoring') ## To add: insurance (SAHIE)
beh_covs <- c('current_smoker_prev','obesity_prev','as_heavy_drinking_prev')
pop_covs <- c('fb','perc_black', 'perc_hispanic') ## To add: perc_native, net_migration
all_covs <- c(ses_covs, med_covs, beh_covs, pop_covs)
run_seven_models <- FALSE
overwrite_full_model <- FALSE
shapley <- TRUE

#cov_domain <- 'all'
#for(sex_option in c(2,1)) {
#for(age_group_option in list(c(45,60), c(25,64), c(65,80))) {
for(mort_domain in c('allcause','drugs','cardio','lung','suicide')) {
#for(cov_domain in c('ses','med','beh','pop')) {
  
  age_group_option <- c(25,64)
  message(paste0('Running sex: ', sex_option, ', mort: ', mort_domain))
  
  # if(cov_domain=='ses') covs <- c('college','log_hh_income','percent_transfers') ## To add: eviction_rate, perc_manufacturing
  # if(cov_domain=='med') covs <- c('log_mds_pc','chr_mammography','chr_diabetes_monitoring') ## To add: insurance (SAHIE)
  # if(cov_domain=='beh') covs <- c('current_smoker_prev','obesity_prev','as_heavy_drinking_prev')
  # if(cov_domain=='pop') covs <- c('fb','perc_25_64','perc_black', 'perc_hispanic') ## To add: perc_native, net_migration
  all_cov_names <- c(ses_covs, med_covs, beh_covs, pop_covs)
  
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
  ## Categorize deaths
  if(mort_domain=='mental') mort[, total_deaths := deaths6]
  if(mort_domain=='cardio') mort[, total_deaths := deaths9]
  if(mort_domain=='drugs') mort <- mort[, total_deaths := apply(.SD, 1, sum), .SDcols=c('deaths17','deaths22','deaths24')]
  if(mort_domain=='lung') mort[, total_deaths := apply(.SD, 1, sum), .SDcols=c('deaths3','deaths10')]
  if(mort_domain=='suicide') mort[, total_deaths := apply(.SD, 1, sum), .SDcols=c('deaths22','deaths23')]
  mort[, nmx := total_deaths / total_pop]
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
 
  ## Scale covariates so coefficients are more interpretable (standard deviations)
  #all_covs <- as.data.table(all_covs %>% mutate_at(funs(scale(.) %>% as.vector), .vars=all_cov_names))

  ## Old method by race
  # mort <- merge(mort, pop, by=c('fips','year','sex','race'))
  # mort <- merge(mort, all_covs, by=c('fips','year','sex','race'))
  # race_option <- ifelse(race=='nhw',0,1)
  # mort <- mort[sex==sex_option & race==race_option, ]
  # mort[, deaths := round(nmx * total_pop)]
  ## Now with all-race, collapse covs using pops
  pop <- merge(pop, all_covs, by=c('fips','year','sex','race'))
  pop <- pop[year %in% year_range, lapply(.SD, weighted.mean, w=total_pop, na.rm=TRUE), by=c('fips','year','sex'), .SDcols=all_cov_names]
  mort <- merge(mort, pop, by=c('fips','year','sex'))

  mort[is.na(total_pop) | total_pop == 0, total_pop := 1]
  mort <- mort[year %in% year_range, ]
  #mort <- mort[sex == sex_option, ]

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

  for(c in all_cov_names) {
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
  
  ## Try running each model with covariate categories, then all all together, then spatial
  ses_covs <- c('college','poverty_all','percent_transfers','manufacturing') ## To add: eviction_rate, perc_manufacturing
  med_covs <- c('log_mds_pc','chr_mammography','chr_diabetes_monitoring') ## To add: insurance (SAHIE)
  beh_covs <- c('current_smoker_prev','obesity_prev','as_heavy_drinking_prev')
  pop_covs <- c('fb','perc_black', 'perc_hispanic') ## To add: perc_native, net_migration
  all_covs <- c(ses_covs, med_covs, beh_covs, pop_covs)
  inla_formulas <- c(paste0('deaths ~ as.factor(sex) + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + as.factor(year)'),
                     paste0('deaths ~ ', paste(ses_covs, collapse = ' + '), ' + as.factor(sex) + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + as.factor(year)'),
                     paste0('deaths ~ ', paste(med_covs, collapse = ' + '), ' + as.factor(sex) + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + as.factor(year)'),
                     paste0('deaths ~ ', paste(beh_covs, collapse = ' + '), ' + as.factor(sex) + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + as.factor(year)'),
                     paste0('deaths ~ ', paste(pop_covs, collapse = ' + '), ' + as.factor(sex) + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + as.factor(year)'),
                     paste0('deaths ~ ', paste(all_covs, collapse = ' + '), ' + as.factor(sex) + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + as.factor(year)'),
                     paste0('deaths ~ ', paste(all_covs, collapse = ' + '), ' + as.factor(sex) + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + as.factor(year) + f(ID,model="bym",graph="FOQ_INLA")'))
  full_mort <- copy(mort)

  d <- full_mort[year %in% c(start_year,end_year), ]
  d <- dcast(d, fips + agegrp ~ year, value.var = c(all_cov_names, 'total_pop', 'metro_region', 'nmx'))
  d[, (paste0('year_', start_year)) := start_year]
  d[, (paste0('year_', end_year)) := end_year]

  ## Try running all age groups separately and binding together.
  #age_groups <- list(c(0,20), c(25,40), c(45,60), c(65,80))
  #age_groups <- list(c(45,60), c(25,64), c(65,80))
  age_groups <- list(age_group_option)
  ## Run all permutations of covariate list
  permutations <- list()
  for(c in 1:length(all_cov_names)) permutations[[c]] <- c(0,1)
  permutations <- as.data.table(expand.grid(permutations))
  setnames(permutations, names(permutations), all_cov_names)
  ## Run all permutations of covariate list
  permutations <- list()
  for(c in 1:length(all_cov_names)) permutations[[c]] <- c(0,1)
  permutations <- as.data.table(expand.grid(permutations))
  setnames(permutations, names(permutations), all_cov_names)
  permutations <- permutations[-1,]
  # Run each covariate permutation
  
  ## GET MODEL COEFFICIENTS FOR ALL MODELS
  if(run_seven_models) {
  get_coefs <- function(f) {
    inla_formula <- inla_formulas[f]
    coefs <- rbindlist(lapply(age_groups, model_permute,
                              data=full_mort, inla_f=inla_formula,
                              coef_file=paste0(sex_option,'_',c),
                              shapley=FALSE,
                              shap_covs=all_cov_names,
                              perm=f))
  }
  all_coefs <- rbindlist(lapply(1:length(inla_formulas), get_coefs), fill=T)
  
  t <- copy(all_coefs)
  t <- t[name!='R2',]
  t[!(name %in% c('DIC','RMSE')), coef := round(coef,2)]
  t[name=='DIC', coef := round(coef)]
  t[name=='RMSE', coef := round(coef,5)]
  t[, coef := as.character(coef)]
  # spacer <- data.table(model=unique(t[, model]), name=rep('zzzz',4), coef='', p=1)
  # t <- rbind(t, spacer)
  t <- dcast(t,  name ~ formula, value.var='coef', fill = '-')
  t[name %in% c('DIC','R2','RMSE',"Moran's I"), name := paste0('zzz',name)]
  t <- t[order(name)]
  t[, name := gsub('zzz','',name)]
  # t[name=='z', model := '']
  # t[name=='z', name := '']
  cov_names <- get_clean_cov_names()
  cov_names[, name := fe]
  t <- merge(t, cov_names[, c('name','cov_name','cov_sort')], all.x=T)
  t <- rbind(t, cov_names[cov_sort %in% c(2,5,8,56,60), c('name','cov_name','cov_sort')], fill=T)
  t <- t[order(cov_sort)]
  t[, cov_sort := NULL]
  t[, cause := mort_domain]
  for(m in as.character(1:7)) t[is.na(get(m)), (m) := '-']
  write.csv(t, paste0('/share/homes/ngraetz/rwfj_tests/', mort_domain, '_', paste(age_group_option, collapse='_'), '_coefficients.csv'), row.names = F)
  }
  
  ## SHAPLEY DECOMPOSITION
  ## Create a neighborhood file using Queens convention.
  if(shapley) {
  message('Creating spatial weight matrix...')
  nb.FOQ <- poly2nb(counties, queen=TRUE)
  ## Create an INLA weight matrix.
  lw.FOQ <- nb2INLA("FOQ_INLA",nb.FOQ)
  an <- data.frame(1:length(counties),counties$fips)
  o <- match(full_mort$fips,an$counties.fips)
  ID <- an$X1.length.counties.[o]
  full_mort[, ID := ID]
  
  mort <- copy(full_mort[agegrp %in% age_group_option[1]:age_group_option[2], ])
  message(unique(mort[, agegrp]))
  
  if(overwrite_full_model) {
    inla_model = inla(as.formula(inla_formulas[[7]]),
                      family = "binomial",
                      data = mort,
                      Ntrials = mort[, total_pop],
                      verbose = FALSE,
                      control.compute=list(config = TRUE, dic = TRUE),
                      control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
                      control.fixed=list(prec.intercept = 0,
                                         prec = 1),
                      num.threads = 10)
    saveRDS(inla_model, paste0('/share/homes/ngraetz/rwfj_tests/', mort_domain, '_', paste(age_group_option, collapse='_'), '_model.RDS'))
  } else {
    inla_model <- readRDS(paste0('/share/homes/ngraetz/rwfj_tests/', mort_domain, '_', paste(age_group_option, collapse='_'), '_model.RDS'))
  }
  
  mort[, inla_pred := inla_model$summary.fitted.values$mean]
  mort[, inla_residual := logit(nmx) - logit(inla_pred)]
  mort[nmx==0, inla_residual := logit(nmx+0.000001) - logit(inla_pred)]
  model_coefs <- make_beta_table(inla_model, mort_domain)
  
  ## Run Shapley decomposition on change over time in ASDR.
  ## Create permutations (2010-1990, 6 changes, total permutations = 2^6 = 64, 32 pairs)
  ## i.e. one pair for poverty is delta_m|PV=2013,IS=1990,CE=1990,FB=1990,time=1990,residual=1990 -
  ##                              delta_m|PV=1990,IS=1990,CE=1990,FB=1990,time=1990,residual=1990
  permutations <- make_permutations(fes = all_covs,
                                    start_year = start_year,
                                    end_year = end_year)
  
  ## Prep and reshape input data from model (all fixed effects of interest + geographic random effects +
  ## time + spatial random effects + intercept + residual, wide by year)
  #d <- copy(mort)
  shap_d <- merge(mort, inla_model$summary.random$ID[c('ID','mean')][1:3108,], by='ID') # First 1:(total spatial units) rows are the COMBINED estimated random effect (spatial + IID). The last set of rows is just spatial.
  setnames(shap_d, 'mean', 'spatial_effect')
  shap_d <- shap_d[year %in% c(start_year,end_year), ]
  shap_d <- dcast(shap_d, fips + agegrp + sex ~ year, value.var = c(all_covs, 'inla_residual', 'inla_pred', 'total_pop', 'metroname','regionname','metro_region','nmx', 'spatial_effect'))
  shap_d[, (paste0('year_', start_year)) := start_year]
  shap_d[, (paste0('year_', end_year)) := end_year]
  
  ## MAKE TABLE 1: ALL COV LEVELS/CHANGES AND NMX
  all_pops <- mort[year %in% c(start_year,end_year), list(total_pop=sum(total_pop)), by=c('fips','sex','year','regionname','metroname')]
  age_start <- min(age_group_option)
  age_end <- max(age_group_option)
  pops <- copy(mort[, c('total_pop','agegrp')])
  age_wts_combined <- pops[, list(pop=sum(total_pop)), by=c('agegrp')]
  totals <- pops[, list(total=sum(total_pop))]
  age_wts_combined <- age_wts_combined[!is.na(agegrp), ]
  age_wts_combined[, total := totals]
  age_wts_combined[, age_wt := pop/total]
  age_wts_combined <- age_wts_combined[, c('age_wt','agegrp')]
  
  table1 <- copy(mort[year %in% c(start_year,end_year), ])
  table1 <- merge(table1, age_wts_combined, by='agegrp')
  table1 <- table1[, lapply(.SD, weighted.mean, w=age_wt, na.rm=TRUE), by=c('fips','sex','year'), .SDcols=c(all_covs,'nmx')] 
  table1 <- merge(table1, all_pops, by=c('fips','sex','year'))
  lisa_data <- table1[, list(nmx=weighted.mean(nmx, total_pop, na.rm=TRUE), total_pop=sum(total_pop)), by=c('fips','year')] 
  table1 <- table1[, lapply(.SD, weighted.mean, w=total_pop, na.rm=TRUE), by=c('regionname','metroname','year'), .SDcols=c(all_covs,'nmx')] 
  ## Make Lisa on age-standardized, both sex mortality rates
  lisa <- plot_lisa(lisa_var = 'nmx',
                    lisa_dt = lisa_data[year==2015, ],
                    lisa_sp = counties,
                    lisa_id = 'fips',
                    lisa_var_name = 'Mortality',
                    sig = 0.05,
                    matrix = 'queen')
  global_morans <- lisa[[3]]
  panel_a_title <- 'A)'
  panel_b_title <- paste0('B)\nGlobal Morans I: ', global_morans)
  title_a <- textGrob(
    label = "A)",
    x = unit(0, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0,
    gp = gpar(fontsize = 16))
  panel_a <- arrangeGrob(lisa[[2]] + theme(legend.position="none"), top = title_a)
  panel_a$vp <- vplayout(1:6, 1:12)
  title_b <- textGrob(
    label = "B)",
    x = unit(0, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0,
    gp = gpar(fontsize = 16))
  panel_b <- arrangeGrob(lisa[[1]] + theme(legend.position="none"), top = title_b)
  panel_b$vp <- vplayout(7:14, 1:12)
  png('/share/code/geospatial/ngraetz/rwjf_counties/results/lisa_figure_2015.png', width = 1200, height = 1400, res = 120)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(14, 12)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y) 
  grid.draw(panel_a)
  grid.draw(panel_b)
  dev.off()
  table1 <- dcast(table1, regionname + metroname ~ year, value.var = c(all_covs,'nmx'))
  for(c in c(all_covs,'nmx')) table1[, (paste0(c,'_change')) := get(paste0(c,'_2015')) - get(paste0(c,'_2000'))]
  setcolorder(table1, sort(names(table1)))
  write.csv(table1, '/share/code/geospatial/ngraetz/rwjf_counties/results/table1.csv')
  
  ## Calculate contribution attributable to each time-varying component. By definition this adds up to total observed change in the outcome
  ## because of inclusion of the residual.
  ## Change decompositions occur at the county-age-level.
  message(paste0('SHAPLEY DECOMPOSITION: ', paste(all_covs,collapse=' ')))
  message(dim(permutations)[1])
  
  library(parallel)
  all_contributions <- rbindlist(parallel::mclapply(c(all_covs, 'year','residual'), calculate_contribution_dev,
                                        fes=all_covs,
                                        start_year=start_year,
                                        end_year=end_year,
                                        dt=shap_d,
                                        coefs=model_coefs,
                                        all_permutations=permutations,
                                        mc.cores=10))
  saveRDS(all_contributions, paste0('/share/homes/ngraetz/rwfj_tests/', mort_domain, '_', paste(age_group_option, collapse='_'), '_contributions.RDS'))
  
  ## Calculate observed age-standardized change at the metro-region for comparison.
  all_pops <- shap_d[, list(total_pop_2015=sum(total_pop_2015), total_pop_2000=sum(total_pop_2000)), by=c('fips','metro_region_2015')]
  age_start <- min(age_group_option)
  age_end <- max(age_group_option)
  pops_2015 <- copy(shap_d[, c('total_pop_2015','agegrp')])
  age_wts_combined <- pops_2015[, list(pop=sum(total_pop_2015)), by=c('agegrp')]
  totals <- pops_2015[, list(total=sum(total_pop_2015))]
  age_wts_combined <- age_wts_combined[!is.na(agegrp), ]
  age_wts_combined[, total := totals]
  age_wts_combined[, age_wt := pop/total]
  age_wts_combined <- age_wts_combined[, c('age_wt','agegrp')]
  obs <- merge(shap_d, age_wts_combined, by='agegrp')
  obs <- obs[, list(nmx_2000=weighted.mean(nmx_2000, age_wt),
                    nmx_2015=weighted.mean(nmx_2015, age_wt)), by=c('fips')]
  obs <- merge(obs, all_pops, by='fips')
  obs <- obs[, list(nmx_2000=weighted.mean(nmx_2000, total_pop_2000, na.rm=T),
                    nmx_2015=weighted.mean(nmx_2015, total_pop_2015, na.rm=T)), by=c('metro_region_2015')]
  obs[, nmx_change := (nmx_2015 - nmx_2000) * 100000]
  decomp_change <- all_contributions[, list(decomp_change=sum(contribution_mort)), by=c('metro_region_2015')]
  decomp_change[, decomp_change := decomp_change * 100000]
  change_compare <- merge(decomp_change, obs, by='metro_region_2015')
  ## CALCULATE COR AND MAE
  ## Change this once we fix Miami-Dade FIPS problem.
  change_compare[metro_region_2015 == 'Lg central metro_South Atlantic', nmx_change := decomp_change]
  cor(change_compare[, c('decomp_change','nmx_change')])
  change_compare[metro_region_2015 != 'Lg central metro_South Atlantic', mean(abs(nmx_change - decomp_change))]
  ggplot() +
    geom_point(data=change_compare,
               aes(x=decomp_change,
                   y=nmx_change)) + 
    geom_abline(slope = 1, intercept = 0) + 
    labs(x='Sum of estimated contributions from all time-varying terms (per 100,000)',
         y='Observed change in mortality rate (per 100,000)') + 
    theme_minimal()
  
  ## Clean up for plotting
  cov_names <- get_clean_cov_names()
  plot_data <- merge(all_contributions, cov_names, by='fe')
  totals <- plot_data[, list(total_change=sum(contribution_mort)), by=c('metro_region_2015')]
  plot_data <- merge(plot_data, totals, by='metro_region_2015')
  plot_data[, metro_region_2015 := gsub('_','',metro_region_2015)]
  plot_data[, metro_region_2015 := gsub('metro','metro ',metro_region_2015)]
  plot_data[, metro_region := factor(metro_region_2015, levels=unique(plot_data$metro_region_2015[order(plot_data[, total_change])]))]
  plot_data[, cov_name := factor(cov_name, levels=unique(plot_data$cov_name[order(plot_data[, cov_sort])]))]
  plot_data[, contribution_mort := contribution_mort*100000]
  plot_data[, total_change := total_change*100000]
  write.csv(plot_data, paste0('/share/code/geospatial/ngraetz/rwjf_counties/results/', mort_domain, '_', paste(age_group_option, collapse='_'), '_decomp.csv'))
  
  pdf(paste0('/share/code/geospatial/ngraetz/rwjf_counties/results/', mort_domain, '_', paste(age_group_option, collapse='_'), '_decomp.pdf'), width = 14, height=12)
  decomp_gg <- ggplot(data=plot_data[metro_region_2015 %in% c('Nonmetro Appalachia','Nonmetro East South Central','Lg central metro Middle Atlantic','Lg central metro South Atlantic','Nonmetro New England'),]) +
    geom_bar(aes(x=metro_region,
                 y=contribution_mort,
                 fill=cov_name),
             color='black',
             stat='identity',
             position='stack',
             width=0.8) + 
    geom_point(aes(x=metro_region,
                   y=total_change),
               size=5) + 
    labs(x = '', y = 'Change in mortality rate per 100,000') + 
    geom_hline(yintercept=0, size=2) + 
    coord_flip() + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
    theme(axis.text.y = element_text(size = 15)) +
    theme(legend.text=element_text(size=12)) + 
    scale_fill_manual(values = c('#deebf7','#9ecae1','#4292c6','#08519c',
                                 '#fc9272','#ef3b2c','#a50f15',
                                 '#a1d99b','#41ab5d','#006d2c',
                                 '#bcbddc','#807dba','#54278f',
                                 '#bdbdbd','#f768a1'), name='Independent\nexplanatory\nfactor')
  print(decomp_gg)
  
  decomp_gg <- ggplot(data=plot_data) +
    geom_bar(aes(x=metro_region,
                 y=contribution_mort,
                 fill=cov_name),
             color='black',
             stat='identity',
             position='stack',
             width=0.8) + 
    geom_point(aes(x=metro_region,
                   y=total_change),
               size=5) + 
    labs(x = '', y = 'Change in mortality rate per 100,000') + 
    geom_hline(yintercept=0, size=2) + 
    coord_flip() + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
    theme(axis.text.y = element_text(size = 15)) +
    theme(legend.text=element_text(size=12)) + 
    # 4 SES, 3 POP, 3 BEHAVIOR, 3 HEALTHCARE, 1 SECULAR, 1 RESIDUAL
    # '#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b',
    # '#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d',
    # '#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b',
    # '#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d',
    scale_fill_manual(values = c('#deebf7','#9ecae1','#4292c6','#08519c',
                                 '#fc9272','#ef3b2c','#a50f15',
                                 '#a1d99b','#41ab5d','#006d2c',
                                 '#bcbddc','#807dba','#54278f',
                                 '#bdbdbd','#f768a1'), name='Independent\nexplanatory\nfactor')
  print(decomp_gg)
  
  dev.off()
  
  }
}
#}
#}
