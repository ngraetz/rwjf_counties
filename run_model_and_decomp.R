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
library(Hmisc)
library(spdep)
source(paste0(repo, 'functions.R'))
source(paste0(repo, 'functions_shapley.R'))
input_dir <- paste0(repo, 'nmx_clean_data/')
cov_dir <- paste0(repo, 'covariate_clean_data/')
out_dir <- 'C:/Users/ngraetz/Dropbox/Penn/papers/rwjf/paa_materials/'

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
covs <- c('college','poverty_all','percent_transfers','fb','percent_unemployment','perc_25_64')
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

model_fips <- 99999
for(fip in unique(mort[metro_region %in% "Nonmetro_West North Central", fips])) {
  message(fip)
  model_fips <- c(model_fips, fip)
  inla_model = inla(as.formula('deaths ~ as.factor(agegrp)'),
                    family = "binomial",
                    data = mort,
                    Ntrials = mort[, total_pop],
                    verbose = FALSE,
                    control.compute=list(config = TRUE, dic = TRUE),
                    control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
                    control.fixed=list(prec.intercept = 0,
                                       prec = 1), num.threads = 2)
}

m <- glm(cbind(round(deaths), round(total_pop-deaths)) ~ as.factor(agegrp) + year, data = mort, family = binomial)
fit_glm <- glm(cbind(round(out_migration), round(total_pop-out_migration)) ~ 
                 lag5_r_size_15_19 + lag5_out_rate + as.factor(name) + year,
               data=model_data[name %in% cor_countries_mort & !is.na(out_migration), ], family=binomial)

## Fit INLA model, save coefficients for table.
## Make tables comparing coefficients
## 1. Metro * Regions, AR1 on year
## 2. Metro * Regions + Covariates, AR1 on year
## 3. Metro * Regions + Covariates, AR1 on year + Besag on county
message('Fitting INLA models...')
mort[, regionname := factor(regionname, levels = c('Pacific', unique(mort[!(regionname %in% 'Pacific'), regionname])))]
mort[, metroname := factor(metroname, levels = c('Lg central metro', unique(mort[!(metroname %in% 'Lg central metro'), metroname])))]
inla_formulas <- c(paste0('deaths ~ as.factor(agegrp) + as.factor(metroname)'),
                   paste0('deaths ~ as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + year'),
                   paste0('deaths ~ ', paste(covs, collapse = ' + '), ' + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + year'),
                   paste0('deaths ~ ', paste(covs, collapse = ' + '), ' + as.factor(agegrp) + as.factor(metroname) + as.factor(regionname) + year + f(ID,model="besag",graph="FOQ_INLA")'),
                   paste0('deaths ~ ', paste(covs, collapse = ' + '), ' + as.factor(agegrp) + as.factor(metro_region) + year + f(ID,model="besag",graph="FOQ_INLA")'))
for(f in 1:4) {
  inla_model = inla(as.formula(inla_formulas[f]),
                    family = "binomial",
                    data = mort[fips=='01001', ],
                    Ntrials = mort[fips=='01001', total_pop],
                    verbose = FALSE,
                    control.compute=list(config = TRUE, dic = TRUE),
                    control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
                    control.fixed=list(prec.intercept = 0,
                                       prec = 1))
  assign(paste0('inla_model_', f), inla_model)
  mort[, (paste0('inla_pred_', f)) := inla_model$summary.fitted.values$mean]
  mort[, (paste0('inla_residual_', f)) := nmx - get(paste0('inla_pred_', f))]
}
all_tables <- lapply(1:4, make_coef_table)
all_coefs <- rbindlist(all_tables, fill=TRUE)
all_dic_morans <- data.table(model = rep(paste0('Model ', 1:4),3),
                             name = c(rep("Global Moran's I",4), rep("DIC",4), rep("RMSE",4)),
                             coef = rep(.00001,12))
for(m in 1:4) {
  lisa <- plot_lisa(lisa_var = paste0('inla_residual_', m),
                    lisa_dt = mort[year==2010,],
                    lisa_sp = counties,
                    lisa_id = 'fips',
                    lisa_var_name = 'Residuals',
                    sig = 0.05)
  all_dic_morans[model == paste0('Model ', m) & name == "Global Moran's I", coef := lisa[[3]]]
  all_dic_morans[model == paste0('Model ', m) & name == "DIC", coef := ifelse(is.nan(get(paste0('inla_model_',m))$dic$dic) | is.infinite(get(paste0('inla_model_',m))$dic$dic), 
                                                                              get(paste0('inla_model_',m))$dic$deviance.mean,
                                                                              get(paste0('inla_model_',m))$dic$dic)]
  all_dic_morans[model == paste0('Model ', m) & name == 'RMSE', coef := mort[, sqrt(weighted.mean(get(paste0('inla_residual_', m))^2, w = total_pop))]]
}
all_coefs <- rbind(all_coefs, all_dic_morans, fill = TRUE)
saveRDS(all_coefs, file = paste0(out_dir, '/sex_', sex_option, '_race_', race_option, '_coef_table_just_inla.RDS'))

lisa <- plot_lisa(lisa_var = 'nmx',
                  lisa_dt = mort[year==2015,],
                  lisa_sp = counties,
                  lisa_id = 'fips',
                  lisa_var_name = 'nMx',
                  sig = 0.05)
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
  
  png(paste0(out_dir, '/lisa_data_',m,'_sex', sex_option, '_race_', race_option, '.png'), width = 1200, height = 1400, res = 120)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(14, 12)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y) 
  #print(lisa[[2]] + ggtitle(panel_a_title) + theme(legend.position="none"), vp = vplayout(1:6, 1:12))
  grid.draw(panel_a)
  #print(lisa[[1]] + ggtitle(panel_b_title) + theme(legend.position="none"), vp = vplayout(7:14, 1:12))
  grid.draw(panel_b)
  dev.off()
  
## Do everything else with the full spatial model and metro*region interaction.
inla_model = inla(as.formula(inla_formulas[5]),
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
saveRDS(coefs, file = paste0(out_dir,'coefs_', output_name, '.RDS'))

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
collapsed[, metro_region_clean := gsub('_',' - ',metro_region)]
collapsed[, metro_region_clean := factor(metro_region_clean, levels=unique(collapsed$metro_region_clean[order(collapsed[, total_change])]))]
pdf(paste0(out_dir,'/decomp_', output_name, '.pdf'), width = 11, height = 11*(2/3))
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
dev.off()

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
