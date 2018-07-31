## Set location of repo.
repo <- 'C:/Users/ngraetz/Documents/repos/rwjf_counties/'

## Load all libraries and functions.
library(INLA)
library(data.table)
library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(gridExtra)
library(grid)
library(plyr)
library(RColorBrewer)
source(paste0(repo, 'functions.R'))
source(paste0(repo, 'functions_shapley.R'))
input_dir <- paste0(repo, 'nmx_clean_data/')
cov_dir <- paste0(repo, 'covariate_clean_data/')

## Set options for this run (data domain and covariates).
## Current datasets: 25-64 ASDR for national NHW male, national NHW female, South NHW male, South NWH female, South NHB male, South NHB female.
## Current possible covariates: 
##    BEA:    "percent_transfers","percent_wage_salary_employment","income_per_capita","total_employees"  
##    BLS:    "labor_force","employed","unemployed","percent_unemployment"
##    SAIPE:  "poverty_all","poverty_0_17","poverty_5_17","median_hh_income"
##    FactFinder (Census + ACS): "fb","less_12","hs","assoc","college"
race <- 'nhw'
sex_option <- 1
domain <- 'national'
covs <- c('college','poverty_all','percent_transfers','percent_unemployment','perc_25_64','fb')
year_range <- c(2000,2010,2015)
start_year <- min(year_range)
end_year <- max(year_range)

## Load master shapefile.
counties <- readOGR("C:/Users/ngraetz/Downloads/cb_2016_us_county_20m", 'cb_2016_us_county_20m')
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
mort <- readRDS(paste0(input_dir, 'asdr_', race, '_', domain, '.RDS'))
pop <- readRDS(paste0(input_dir, race, '_total_pop.RDS'))
all_covs <- readRDS(paste0(cov_dir, 'combined_covs.RDS'))
mort <- merge(mort, pop, by=c('fips','year','sex','race'))
mort <- merge(mort, all_covs, by=c('fips','year','sex','race'))
race_option <- ifelse(race=='nhw',0,1)
mort <- mort[sex==sex_option & race==race_option, ]
mort[is.na(total_pop) | total_pop == 0, total_pop := 1]
mort[, deaths := round(nmx * total_pop)]
mort <- mort[year %in% year_range, ]

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

## Drop if missing any covariates and log. Combine nonmetro categories and create interaction variable.
message('Counties before dropping missing values in any year: ', length(unique(mort[, fips])))
drop_counties <- c()
for(fe in covs) drop_counties <- c(drop_counties, mort[is.na(get(fe)), fips])
drop_counties <- c(drop_counties, mort[is.na(metroname) | is.na(regionname), fips])
mort <- mort[!(fips %in% unique(drop_counties)), ]
message('Counties after dropping missing values in any year: ', length(unique(mort[, fips])))

## Create a neighborhood file using Queens convention.
message('Creating spatial weight matrix...')
nb.FOQ <- poly2nb(counties, queen=TRUE)
## Create an INLA weight matrix.
lw.FOQ <- nb2INLA("FOQ_INLA",nb.FOQ)
an <- data.frame(1:length(counties),counties$fips)
o <- match(mort$fips,an$counties.fips)
ID <- an$X1.length.counties.[o]
mort[, ID := ID]

## Fit INLA model, save coefficients for table.
inla_formula <- as.formula(paste0('deaths ~ ',
                           paste(covs, collapse = ' + '),
                           ' + as.factor(metro_region) + year + f(ID,model="besag",graph="FOQ_INLA")'))
message('Fitting INLA model...')
inla_model = inla(inla_formula,
                    family = "binomial",
                    data = mort,
                    Ntrials = mort[, total_pop],
                    verbose = FALSE,
                    control.compute=list(config = TRUE, dic = TRUE),
                    control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
                    control.fixed=list(prec.intercept = 0,
                                       prec = 1))

## Make full prediction, full residual, and load posterior mean for all components.
mort[, inla_pred := inla_model$summary.fitted.values$mean]
mort[, inla_residual := logit(nmx) - logit(inla_pred)]
mort[nmx==0, inla_residual := logit(nmx+0.000001) - logit(inla_pred)]
coefs <- make_beta_table(inla_model, paste0(race,' ',sex_option,' ',domain))

## Run Shapley decomposition on change over time in ASDR.
## Create permutations (2010-1990, 6 changes, total permutations = 2^6 = 64, 32 pairs)
## i.e. one pair for poverty is delta_m|PV=2013,IS=1990,CE=1990,FB=1990,time=1990,residual=1990 - 
##                              delta_m|PV=1990,IS=1990,CE=1990,FB=1990,time=1990,residual=1990
permutations <- make_permutations(fes = covs,
                                  start_year = start_year,
                                  end_year = end_year)

## Prep and reshape input data from model (all fixed effects of interest + geographic random effects + 
## time + spatial random effects + intercept + residual, wide by year)
d <- copy(mort)
d <- merge(d, inla_model$summary.random$ID[c('ID','mean')], by='ID')
setnames(d, 'mean', 'spatial_effect')
d <- d[year %in% c(start_year,end_year), ]
d <- dcast(d, fips ~ year, value.var = c(covs, 'inla_residual', 'inla_pred', 'total_pop', 'metro_region', 'nmx', 'spatial_effect'))
d[, (paste0('year_', start_year)) := start_year]
d[, (paste0('year_', end_year)) := end_year]

## Calculate contribution attributable to each time-varying component. By definition this adds up to total observed change in the outcome
## because of inclusion of the residual.
all_contributions <- rbindlist(lapply(c(covs, 'year','residual'), calculate_contribution,
                                      fes=covs,
                                      start_year=start_year,
                                      end_year=end_year))

## Collapse contributions to metro-region.
collapsed <- merge(all_contributions, unique(d[, c('fips','metro_region_2015','total_pop_2015','total_pop_2000')]), by='fips')
collapsed <- collapsed[, list(contribution_mort=wtd.mean(contribution_mort, total_pop_2015, na.rm=TRUE),
                                              total_pop_2015=sum(total_pop_2015, na.rm=TRUE)),
                                       by=c('metro_region_2015','fe')]
collapsed <- merge(collapsed, collapsed[, list(total_change=sum(contribution_mort, na.rm=TRUE)), by='metro_region_2015'], by='metro_region_2015')
collapsed[, metro_region := factor(metro_region_2015, levels=unique(collapsed$metro_region_2015[order(collapsed[, total_change])]))]

## Format for plotting
cov_names <- get_clean_cov_names()
collapsed[, fe := factor(fe, levels=c(covs,"year","residual"))]
collapsed <- merge(collapsed, cov_names, by='fe')
collapsed[, cov_name := factor(cov_name, levels=cov_names[fe %in% collapsed[, fe], cov_name])]

## Plot results of decomposition by metro-region.
sex_name <- ifelse(sex==1,'Male','Female')
pdf(paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/decomp_', paste0(capitalize(domain), ' ', toupper(race), ' ', capitalize(sex_name)), '.pdf'), width = 12, height = 8)
ggplot() +
  geom_bar(data=collapsed,
           aes(x=metro_region,
               y=contribution_mort,
               fill=cov_name),
           color='black',
           stat='identity') + 
  geom_point(data=collapsed[, list(contribution_mort=sum(contribution_mort)), by=c('metro_region')],
             aes(x=metro_region,
                 y=contribution_mort),
             size=3) + 
  labs(x = 'Metro/region', y = 'Change in ASDR', title = paste0('Change in ', capitalize(domain), ' ', toupper(race), ' ', capitalize(sex_name), ' ASDR, ',start_year,'-',end_year)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = brewer.pal(length(unique(collapsed[, cov_name])),'Spectral'), name='Component')
dev.off()

