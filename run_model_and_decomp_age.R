## Set location of repo.
repo <- '/share/code/geospatial/ngraetz/rwjf_counties/'

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
race <- 'nhw'
sex_option <- 1
domain <- 'national'
cov_domain <- 'all'
#for(cov_domain in c('med','beh','pop')) {
if(cov_domain=='ses') covs <- c('college','poverty_all','log_hh_income','percent_transfers','percent_unemployment') ## To add: eviction_rate, perc_manufacturing
if(cov_domain=='med') covs <- c('log_mds_pc','chr_mammography','chr_diabetes_monitoring') ## To add: insurance (SAHIE)
if(cov_domain=='beh') covs <- c('as_diabetes_prev','current_smoker_prev','obesity_prev','as_heavy_drinking_prev')
if(cov_domain=='pop') covs <- c('fb','perc_25_64') ## To add: perc_black, perc_hispanic, perc_native, net_migration
if(cov_domain=='all') covs <- c('college','poverty_all','log_hh_income','percent_transfers','percent_unemployment',
                                'log_mds_pc','chr_mammography','chr_diabetes_monitoring',
                                'as_diabetes_prev','current_smoker_prev','obesity_prev','as_heavy_drinking_prev',
                                'fb','perc_25_64') ## To add: perc_black, perc_hispanic, perc_native, net_migration
year_range <- c(2000,2010,2015)
plot_trends <- FALSE

start_year <- min(year_range)
end_year <- max(year_range)
sex_name <- ifelse(sex_option==1,'Male','Female')
output_name <- paste0(capitalize(domain), ' ', toupper(race), ' ', capitalize(sex_name))

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
mort <- readRDS(paste0(input_dir, 'age_specific_mort.RDS'))
setnames(mort, c('total_deaths','pooled_year'), c('deaths','year'))
pop_nhw <- readRDS(paste0(input_dir, 'nhw', '_total_pop.RDS'))
pop_nhb <- readRDS(paste0(input_dir, 'nhb', '_total_pop.RDS'))
pop <- rbind(pop_nhw,pop_nhb)
pop <- readRDS(paste0(input_dir, race, '_total_pop.RDS'))
all_covs <- readRDS(paste0(cov_dir, 'combined_covs.RDS'))
all_covs[, log_mds_pc := log(mds_pc+0.01)]

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

## Drop if missing any covariates and log. Combine nonmetro categories and create interaction variable.
# message('Counties before dropping missing values in any year: ', length(unique(mort[, fips])))
# drop_counties <- c()
# for(fe in covs) drop_counties <- c(drop_counties, mort[is.na(get(fe)), fips])
# drop_counties <- c(drop_counties, mort[is.na(metroname) | is.na(regionname), fips])
# mort <- mort[!(fips %in% unique(drop_counties)), ]
# message('Counties after dropping missing values in any year: ', length(unique(mort[, fips])))

## Create trend plots by metro (color) and region (facet)
if(plot_trends==TRUE) {
  trends <- mort[, lapply(.SD, weighted.mean, w=total_pop, na.rm=TRUE), .SDcols=covs, by=c('metroname','regionname','year')]
  cov_names <- get_clean_cov_names()
  pdf(paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/trend_plots/', domain, '_', race, '_', sex_name, '_trends.pdf'), width = 12, height = 8)
  for(c in covs) {
    gg <- ggplot() + 
      geom_line(data=trends,
                aes(x=year,
                    y=get(c),
                    color=metroname),
                size=2) +
      theme_minimal() + 
      scale_color_manual(values = brewer.pal(length(unique(collapsed[, cov_name])),'Dark2'), name='Metro') + 
      labs(x='Year',y=cov_names[fe==c, cov_name], title=cov_names[fe==c, cov_name]) + 
      facet_wrap(~regionname)
    print(gg)
  }
  dev.off()
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
inla_formulas <- c(paste0('deaths ~ as.factor(agegrp) + as.factor(metroname)'),
                   paste0('deaths ~ as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + year'),
                   paste0('deaths ~ ', paste(covs, collapse = ' + '), ' + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + year'),
                   paste0('deaths ~ ', paste(covs, collapse = ' + '), ' + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + year + f(ID,model="besag",graph="FOQ_INLA")'),
                   paste0('deaths ~ ', paste(covs, collapse = ' + '), ' + as.factor(agegrp) + as.factor(metro_region) + year + f(ID,model="besag",graph="FOQ_INLA")'),
                   paste0('deaths ~ ', paste(paste0(covs,'*as.factor(broad_age)'), collapse = ' + '), ' + as.factor(agegrp) + as.factor(metro_region) + year + f(ID,model="besag",graph="FOQ_INLA")'))
# for(f in 1:4) {
#   inla_model = inla(as.formula(inla_formulas[f]),
#                     family = "binomial",
#                     data = mort,
#                     Ntrials = mort[, total_pop],
#                     verbose = FALSE,
#                     control.compute=list(config = TRUE, dic = TRUE),
#                     control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
#                     control.fixed=list(prec.intercept = 0,
#                                        prec = 1),
#                     num.threads = 4)
#   assign(paste0('inla_model_', f), inla_model)
#   mort[, (paste0('inla_pred_', f)) := inla_model$summary.fitted.values$mean]
#   mort[, (paste0('inla_residual_', f)) := nmx - get(paste0('inla_pred_', f))]
# }
# all_tables <- lapply(1:4, make_coef_table)
# all_coefs <- rbindlist(all_tables, fill=TRUE)
# all_dic_morans <- data.table(model = rep(paste0('Model ', 1:4),3),
#                              name = c(rep("Global Moran's I",4), rep("DIC",4), rep("RMSE",4)),
#                              coef = rep(.00001,12))
# for(m in 1:4) {
#   lisa <- plot_lisa(lisa_var = paste0('inla_residual_', m),
#                     lisa_dt = mort[year==2010,],
#                     lisa_sp = counties,
#                     lisa_id = 'fips',
#                     lisa_var_name = 'Residuals',
#                     sig = 0.05)
#   all_dic_morans[model == paste0('Model ', m) & name == "Global Moran's I", coef := lisa[[3]]]
#   all_dic_morans[model == paste0('Model ', m) & name == "DIC", coef := ifelse(is.nan(get(paste0('inla_model_',m))$dic$dic) | is.infinite(get(paste0('inla_model_',m))$dic$dic), 
#                                                                               get(paste0('inla_model_',m))$dic$deviance.mean,
#                                                                               get(paste0('inla_model_',m))$dic$dic)]
#   all_dic_morans[model == paste0('Model ', m) & name == 'RMSE', coef := mort[, sqrt(weighted.mean(get(paste0('inla_residual_', m))^2, w = total_pop))]]
# }
# all_coefs <- rbind(all_coefs, all_dic_morans, fill = TRUE)
# saveRDS(all_coefs, file = paste0(out_dir, '/sex_', sex_option, '_race_', race_option, '_coef_table_just_inla.RDS'))
# 
# lisa <- plot_lisa(lisa_var = 'nmx',
#                   lisa_dt = mort[year==2015,],
#                   lisa_sp = counties,
#                   lisa_id = 'fips',
#                   lisa_var_name = 'nMx',
#                   sig = 0.05)
# global_morans <- lisa[[3]]
# panel_a_title <- 'A)'
# panel_b_title <- paste0('B)\nGlobal Morans I: ', global_morans)
# title_a <- textGrob(
#   label = "A)",
#   x = unit(0, "lines"), 
#   y = unit(0, "lines"),
#   hjust = 0, vjust = 0,
#   gp = gpar(fontsize = 16))
# panel_a <- arrangeGrob(lisa[[2]] + theme(legend.position="none"), top = title_a)
# panel_a$vp <- vplayout(1:6, 1:12)
# title_b <- textGrob(
#   label = "B)",
#   x = unit(0, "lines"), 
#   y = unit(0, "lines"),
#   hjust = 0, vjust = 0,
#   gp = gpar(fontsize = 16))
# panel_b <- arrangeGrob(lisa[[1]] + theme(legend.position="none"), top = title_b)
# panel_b$vp <- vplayout(7:14, 1:12)
# 
# png(paste0(out_dir, '/lisa_data_',m,'_sex', sex_option, '_race_', race_option, '.png'), width = 1200, height = 1400, res = 120)
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(14, 12)))
# vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y) 
# #print(lisa[[2]] + ggtitle(panel_a_title) + theme(legend.position="none"), vp = vplayout(1:6, 1:12))
# grid.draw(panel_a)
# #print(lisa[[1]] + ggtitle(panel_b_title) + theme(legend.position="none"), vp = vplayout(7:14, 1:12))
# grid.draw(panel_b)
# dev.off()

## Do everything else with the full spatial model and metro*region interaction.
## With all-age covariates             DIC = 1025807.50
## With broad-age-specific covariates  DIC =  982275.99
full_mort <- copy(mort)

# mort <- copy(full_mort)
# mort <- mort[agegrp %in% 45:60, ]
# inla_model = inla(as.formula(inla_formulas[5]),
#                     family = "binomial",
#                     data = mort,
#                     Ntrials = mort[, total_pop],
#                     verbose = FALSE,
#                     control.compute=list(config = TRUE, dic = TRUE),
#                     control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
#                     control.fixed=list(prec.intercept = 0,
#                                        prec = 1),
#                     num.threads = 4)

## Make full prediction, full residual, and load posterior mean for all components.
# mort[, inla_pred := inla_model$summary.fitted.values$mean]
# mort[, inla_residual := logit(nmx) - logit(inla_pred)]
# mort[nmx==0, inla_residual := logit(nmx+0.000001) - logit(inla_pred)]
# coefs <- make_beta_table(inla_model, paste0(race,' ',sex_option,' ',domain))
#saveRDS(coefs, file = paste0(out_dir,'coefs_', output_name, '.RDS'))

## Plot fitted age curve
# age_curve <- copy(coefs[grep('agegrp',name),])
# age_curve[, odds := exp(coef + coefs[name=='(Intercept)', coef] + coefs[name=='year', coef * 2015])]
# age_curve[, coef := odds / (1 + odds)]
# age_curve[, age := as.numeric(gsub('as.factor\\(agegrp\\)','',name))]
# age_obs <- mort[, list(deaths=sum(deaths),total_pop=sum(total_pop)), by=c('agegrp')]
# age_obs[, nmx := deaths / total_pop]
# #pdf(paste0(out_dir,'/fitted_age_pattern.pdf'))
# ggplot() + 
#   geom_line(data=age_curve,
#             aes(x=age,
#                 y=coef),
#             color='red') + 
#   geom_point(data=age_obs,
#             aes(x=agegrp,
#                 y=nmx),
#             color='black') + theme_minimal()
#dev.off()

## Run Shapley decomposition on change over time in ASDR.
## Create permutations (2010-1990, 6 changes, total permutations = 2^6 = 64, 32 pairs)
## i.e. one pair for poverty is delta_m|PV=2013,IS=1990,CE=1990,FB=1990,time=1990,residual=1990 - 
##                              delta_m|PV=1990,IS=1990,CE=1990,FB=1990,time=1990,residual=1990
# permutations <- make_permutations(fes = covs,
#                                   start_year = start_year,
#                                   end_year = end_year)

## Prep and reshape input data from model (all fixed effects of interest + geographic random effects + 
## time + spatial random effects + intercept + residual, wide by year)
# d <- copy(full_mort)
# d <- merge(d, inla_model$summary.random$ID[c('ID','mean')], by='ID')
# setnames(d, 'mean', 'spatial_effect')
# d <- d[year %in% c(start_year,end_year), ]
# d <- dcast(d, fips + agegrp ~ year, value.var = c(covs, 'inla_residual', 'inla_pred', 'total_pop', 'metro_region', 'nmx', 'spatial_effect'))
# d[, (paste0('year_', start_year)) := start_year]
# d[, (paste0('year_', end_year)) := end_year]

## Calculate contribution attributable to each time-varying component. By definition this adds up to total observed change in the outcome
## because of inclusion of the residual.
## Change decompositions occur at the county-age-level.
# all_contributions <- rbindlist(lapply(c(covs, 'year','residual'), calculate_contribution,
#                                       fes=covs,
#                                       start_year=start_year,
#                                       end_year=end_year,
#                                       d=d))

## Try running all age groups separately and binding together.
age_groups <- list(c(0,20), c(25,40), c(45,60), c(65,85))
all_contributions <- rbindlist(lapply(age_groups, shapley_ages,
                                      data=full_mort, inla_f=inla_formulas[5],
                                      coef_file=paste0(sex_option,'_',cov_domain),
                                      shapley=TRUE,
                                      shap_covs=covs))
#}

## Scatter total predicted change from decomp (sum of contributions) with observed change.
plot_metro_regions <- c('Lg central metro_Middle Atlantic',
                        'Nonmetro_East South Central')
decomp_change <- all_contributions[, list(decomp_change=sum(contribution_mort)), by=c('fips','agegrp')]
decomp_change <- merge(decomp_change, d[, c('fips','agegrp','nmx_2000','nmx_2015','total_pop_2015','metro_region_2015')], by=c('fips','agegrp'))
decomp_change[, obs_change := nmx_2015 - nmx_2000]
round(cor(decomp_change[, c('decomp_change','obs_change')], use='complete.obs')[1,2],2)
decomp_change <- decomp_change[, list(decomp_change=weighted.mean(decomp_change,total_pop_2015,na.rm=T)), by=c('agegrp','metro_region_2015')]
setnames(decomp_change, 'metro_region_2015', 'metro_region')
decomp_change <- decomp_change[metro_region %in% plot_metro_regions, ]
all_contributions_age <- merge(full_mort[year==2015, c('agegrp','fips','total_pop','metro_region')], all_contributions, by=c('agegrp','fips'))
all_contributions_age <- all_contributions_age[, list(contribution_mort=weighted.mean(contribution_mort, total_pop, na.rm=T)), by=c('fe','agegrp','metro_region')]
all_contributions_age <- all_contributions_age[metro_region %in% plot_metro_regions, ]
ggplot() +
  geom_bar(data=all_contributions_age,
           aes(x=as.factor(agegrp),
               y=contribution_mort*100000,
               fill=fe),
           color='black',
           stat='identity') + 
  geom_point(data=decomp_change,
             aes(x=as.factor(agegrp),
                 y=decomp_change*100000),
             size=3) + 
  labs(x = '', y = 'Change in mortality rate (per 100,000)', title = paste0('Change in mortality rate by age, 2000-2015.')) + 
  theme_minimal() +
  facet_wrap(~metro_region) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
  scale_fill_manual(values = brewer.pal(length(unique(all_contributions[, fe])),'Spectral'), name='Component')

## Create age-standardized contributions over some age range using 2000 Census age structure.
age_start <- 45
age_end <- 60
pops_2000 <- full_mort[year==2000 & agegrp >= age_start & agegrp <= age_end, ]
age_wts_combined <- pops_2000[, list(pop=sum(total_pop)), by=c('agegrp','sex')]
totals <- pops_2000[, list(total=sum(total_pop)), by=c('sex')]
age_wts_combined <- age_wts_combined[!is.na(agegrp), ]
age_wts_combined <- merge(age_wts_combined, totals, by='sex')
age_wts_combined[, age_wt := pop/total]
age_wts_combined <- age_wts_combined[, c('age_wt','agegrp')]
collapsed <- merge(all_contributions, age_wts_combined, by='agegrp')
collapsed <- collapsed[, list(contribution_mort=weighted.mean(contribution_mort, age_wt)), by=c('fips','fe')]

## Collapse contributions to metro-region.
all_pops <- d[agegrp >= age_start & agegrp <= age_end, list(total_pop_2015=sum(total_pop_2015), total_pop_2000=sum(total_pop_2000)), by=c('fips','metro_region_2015')]
#collapsed <- merge(all_contributions, unique(all_pops[, c('fips','metro_region_2015','total_pop_2015','total_pop_2000')]), by='fips')
collapsed <- merge(collapsed, all_pops, by='fips')
collapsed <- collapsed[, list(contribution_mort=wtd.mean(contribution_mort, total_pop_2015, na.rm=TRUE),
                              total_pop_2015=sum(total_pop_2015, na.rm=TRUE)),
                       by=c('metro_region_2015','fe')]
collapsed <- merge(collapsed, collapsed[, list(total_change=sum(contribution_mort, na.rm=TRUE)), by='metro_region_2015'], by='metro_region_2015')
collapsed[, metro_region := factor(metro_region_2015, levels=unique(collapsed$metro_region_2015[order(collapsed[, total_change])]))]

## Calculate observed age-standardized change at the metro-region for comparison.
obs <- merge(d, age_wts_combined, by='agegrp')
obs <- obs[, list(nmx_2000=weighted.mean(nmx_2000, age_wt),
                  nmx_2015=weighted.mean(nmx_2015, age_wt)), by=c('fips')]
obs <- merge(obs, all_pops, by='fips')
obs <- obs[, list(nmx_2000=weighted.mean(nmx_2000, total_pop_2000),
                  nmx_2015=weighted.mean(nmx_2015, total_pop_2015)), by=c('metro_region_2015')]
obs[, nmx_change := nmx_2015 - nmx_2000]

## Format for plotting
cov_names <- get_clean_cov_names()
collapsed[, fe := factor(fe, levels=c(covs,"year","residual"))]
collapsed <- merge(collapsed, cov_names, by='fe')
collapsed[, cov_name := factor(cov_name, levels=cov_names[fe %in% collapsed[, fe], cov_name])]

## Plot results of decomposition by metro-region.
collapsed[, metro_region_clean := gsub('_',' - ',metro_region)]
collapsed[, metro_region_clean := factor(metro_region_clean, levels=unique(collapsed$metro_region_clean[order(collapsed[, total_change])]))]
#pdf(paste0(out_dir,'/decomp_', output_name, '.pdf'), width = 11, height = 11*(2/3))
ggplot() +
  geom_bar(data=collapsed,
           aes(x=metro_region_clean,
               y=contribution_mort*100000,
               fill=cov_name),
           color='black',
           stat='identity') + 
  geom_point(data=collapsed[, list(contribution_mort=sum(contribution_mort)), by=c('metro_region_clean')],
             aes(x=metro_region_clean,
                 y=contribution_mort*100000),
             size=3) + 
  labs(x = '', y = 'Change in ASDR', title = paste0('Change in ', output_name, ' ASDR, ',start_year,'-',end_year)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
  scale_fill_manual(values = brewer.pal(length(unique(collapsed[, cov_name])),'Spectral'), name='Component')
#dev.off()

## Plot maps of each contribution by county.
pdf(paste0(out_dir,'/maps_', output_name, '.pdf'), height=6, width=9)
all_contributions[, contribution_mort_per_100000 := contribution_mort * 100000]
top_code <- quantile(all_contributions[, contribution_mort_per_100000], probs=0.90, na.rm=T)
bottom_code <- quantile(all_contributions[, contribution_mort_per_100000], probs=0.10, na.rm=T)
for(c in c(covs)) {
  m <- make_county_map(map_dt = all_contributions[fe==c, ],
                       map_sp = counties,
                       map_var = 'contribution_mort_per_100000',
                       legend_title = 'Contribution',
                       high_is_good = FALSE,
                       map_title = paste0('Contribution of ', cov_names[fe==c, cov_name], ' to change in ASDR, ', start_year,'-',end_year),
                       map_limits = c(-20,50),
                       diverge = TRUE)
  m <- m + scale_fill_gradientn(guide = guide_legend(title = 'Contribution'),
                                limits = c(-20,50),
                                breaks = c(-20,-10,0,10,20,30,40,50),
                                colours=rev(brewer.pal(10,'Spectral')),
                                values=c(-20,0,50), na.value = "#000000", rescaler = function(x, ...) x, oob = identity)
  print(m)
}
c <- 'residual'
top_code <- quantile(all_contributions[, contribution_mort], probs=0.99, na.rm=T)
bottom_code <- quantile(all_contributions[, contribution_mort], probs=0.1, na.rm=T)
m <- make_county_map(map_dt = all_contributions[fe==c, ],
                     map_sp = counties,
                     map_var = 'contribution_mort',
                     legend_title = 'Contribution',
                     high_is_good = FALSE,
                     map_title = paste0('Contribution of ', cov_names[fe==c, cov_name], ' to change in ASDR, ', start_year,'-',end_year),
                     map_limits = c(bottom_code,top_code),
                     diverge = TRUE)
print(m)
dev.off()


pdf(paste0(out_dir,'/female_transfers.pdf'), height=6, width=9)
c <- 'percent_transfers'
m <- make_county_map(map_dt = all_contributions[fe==c, ],
                     map_sp = counties,
                     map_var = 'contribution_mort',
                     legend_title = 'Contribution',
                     high_is_good = FALSE,
                     map_title = paste0('Contribution of ', cov_names[fe==c, cov_name], ' to change in ASDR, ', start_year,'-',end_year),
                     map_limits = c(bottom_code,0.0007),
                     diverge = TRUE)
print(m)
dev.off()

c <- 'fb'
m <- make_county_map(map_dt = all_contributions[fe==c, ],
                     map_sp = counties,
                     map_var = 'contribution_mort',
                     legend_title = 'Contribution',
                     high_is_good = FALSE,
                     map_title = paste0('Contribution of ', cov_names[fe==c, cov_name], ' to change in ASDR, ', start_year,'-',end_year),
                     map_limits = c(-0.0002,0.0007),
                     diverge = TRUE)
