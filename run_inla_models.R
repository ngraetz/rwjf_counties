library(raster)
library(rgdal)
library(rgeos)
library(data.table)
library(maptools)
library(ggplot2)
library(spdep)
library(INLA)
library(spdep)
library(gridExtra)
library(grid)
library(GWmodel)
library(spgwr)
library(tables)
library(spdep)
source('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/functions.R')

for(race_option in c(0,1)) {
  for(sex_option in c(1,2)) {

## Set options
#region_subset <- c('East South Central','Appalachia')
region_subset <- NULL
year_as_factor <- TRUE
#race_option <- 1
#sex_option <- 1
race_title <- ifelse(race_option==1, 'Black', 'White')
pop_option <- 'all'
output_dir <- paste0('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/', race_title, '_', pop_option)
plot_dir <- paste0(output_dir, '/plots_fixed_region')
dir.create(output_dir)
dir.create(plot_dir)

## Load shapefile
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

## Load nmx data
if(pop_option=='all') nmx_list <- readRDS("C:/Users/ngraetz/Documents/repos/spatial_demography_2018/nmx_by_sex_nwh_black_county.rds")
if(pop_option=='25_65') nmx_list <- readRDS("C:/Users/ngraetz/Documents/repos/spatial_demography_2018/data_byrace_25_65.rds")
nmx <- nmx_list[[1]]

## Load covariates
covs <- c('eqi','exposures_by_income','exposures_by_nativity','pop_composition','typologies')
covs <- lapply(paste0('C:/Users/ngraetz/Downloads/',covs,'.csv'), fread)

## Format Brown exposure data
brown <- covs[[2]]
brown <- brown[, c('fips','twh1k_pv', 'tbh1k_pv', 'tbh1k_is', 'twh1k_is', 'twh1k_ce', 'tbh1k_ce',
                   'twh0k_pv', 'tbh0k_pv', 'tbh0k_is', 'twh0k_is', 'twh0k_ce', 'tbh0k_ce',
                   'twh9k_pv', 'tbh9k_pv', 'tbh9k_is', 'twh9k_is', 'twh9k_ce', 'tbh9k_ce')]
brown_pop <- covs[[3]]
brown_pop <- brown_pop[, c('fips', 'fb9p_kt', 'fb0p_kt', 'fb1p_kt',
                           'fw9p_kt', 'fw0p_kt', 'fw1p_kt')]
brown <- merge(brown, brown_pop, by='fips')
brown = melt(brown, id.vars = c("fips"),
             measure.vars = c('twh1k_pv', 'tbh1k_pv', 'tbh1k_is', 'twh1k_is', 'twh1k_ce', 'tbh1k_ce',
                              'twh0k_pv', 'tbh0k_pv', 'tbh0k_is', 'twh0k_is', 'twh0k_ce', 'tbh0k_ce',
                              'twh9k_pv', 'tbh9k_pv', 'tbh9k_is', 'twh9k_is', 'twh9k_ce', 'tbh9k_ce',
                              'fb9p_kt', 'fb0p_kt', 'fb1p_kt',
                              'fw9p_kt', 'fw0p_kt', 'fw1p_kt'))
brown[, variable := as.character(variable)]
brown[nchar(variable)==8, year := as.numeric(substr(variable,4,4))]
brown[nchar(variable)==7, year := as.numeric(substr(variable,3,3))]
brown[year==1, year := 2010]
brown[year==0, year := 2000]
brown[year==9, year := 1990]
brown[nchar(variable)==7, variable := paste0(substr(variable, 1, 2), substr(variable, 4, nchar(variable)))]
brown[nchar(variable)==8, variable := paste0(substr(variable, 1, 3), substr(variable, 5, nchar(variable)))]
brown <- dcast(brown, year + fips ~ variable, value.var = 'value')

## Load county health rankings data (RWJF)
rwjf <- fread('C:/Users/ngraetz/Documents/Penn/papers/county_health_rankings_2018.csv')
setnames(rwjf, 'FIPS', 'fips')
rwjf[, fips := as.character(fips)]
rwjf[nchar(fips)==4, fips := paste0('0',fips)]
rwjf <- rwjf[County!="", ]
setnames(rwjf,
         c('X..Obese',
           'X..Smokers',
           'X..Excessive.Drinking',
           'Z.Score.15', # People per PCP
           'Z.Score.17', # People per MHP
           'X..Uninsured.1',
           'X..Unemployed.1',
           'X..Severe.Housing.Problems'),
         c('Obese',
           'Smokers',
           'Excessive.Drinking',
           'PCP.Z',
           'MHP.Z',
           'Uninsured',
           'Unemployed',
           'Severe.Housing.Problems'))
rwjf <- rwjf[, c('fips','Obese','Smokers','Excessive.Drinking','PCP.Z','MHP.Z','Uninsured','Unemployed','Severe.Housing.Problems')]
#rwjf[, PCP.Ratio := as.numeric(gsub(':1','',PCP.Ratio))]
#rwjf[, MHP.Ratio := as.numeric(gsub(':1','',MHP.Ratio))]
for(v in c('Obese','Smokers','Excessive.Drinking','Uninsured','Unemployed','Severe.Housing.Problems')) {
  rwjf[, (v) := get(v) / 100]
  message(paste0(v, ' ', length(rwjf[is.na(get(v)), get(v)])))
}
for(v in c('PCP.Z','MHP.Z')) {
  message(paste0(v, ' ', length(rwjf[is.na(get(v)), get(v)])))
}
# pdf(paste0(output_dir, '/RWJF_county_ranks_data.pdf'), height=6, width=9)
# for(v in c('Obese','Smokers','Excessive.Drinking','PCP.Z','MHP.Z','Uninsured','Unemployed','Severe.Housing.Problems')) {
# m <- make_county_map(map_dt = model_data_black[year==2010,],
#                            map_sp = counties,
#                            map_var = v,
#                            legend_title = v,
#                            high_is_good = ifelse(v %in% c('PCP.Z','MHP.Z'), TRUE, FALSE),
#                            map_title = v)
# print(m)
# }
# dev.off()

## Compile model dataset
model_data <- merge(nmx, covs[[1]][, c('fips','state','air_EQI_22July2013','water_EQI_22July2013','land_EQI_22July2013','built_EQI_22July2013')])
model_data <- merge(model_data, brown, by=c('fips', 'year'))
model_data <- merge(model_data, rwjf, by='fips')

## Merge populations
if(pop_option=='all') all_pops <- fread("C:/Users/ngraetz/Downloads/all_pops_fixed.csv")
if(pop_option=='25_65') all_pops <- fread("C:/Users/ngraetz/Downloads/all_pops_25_65.csv")
model_data <- merge(model_data, all_pops[, c('year','fips','race','total_pop','total_deaths')], by=c('year','fips','race'), all.x=TRUE)

## Merge metro codes
metro_codes <- fread('C:/Users/ngraetz/Downloads/FIPSmetroregion.csv')
metro_codes[, fips := as.character(fips)]
metro_codes[nchar(fips)==4, fips := paste0('0',fips)]
model_data <- merge(model_data, metro_codes[, c('fips','metroname','regionname','metro','region')], by=c('fips'), all.x=TRUE)

## Rescale covariates 
cov_names <- c('air_EQI_22July2013','water_EQI_22July2013','land_EQI_22July2013','built_EQI_22July2013',
               'twhk_pv', 'tbhk_pv', 'tbhk_is', 'twhk_is', 'twhk_ce', 'tbhk_ce', 'fwp_kt', 'fbp_kt')
for(c in cov_names) model_data[, (c) := scales::rescale(get(c))]

## Subset to race/sex and convert rate to death counts ("age-standardized" death counts, not raw deaths)
model_data_black <- model_data[race==race_option & sex==sex_option, ]
model_data_black[is.na(total_pop) | total_pop == 0, total_pop := 1]
model_data_black[, deaths := round(nmx * total_pop)]
saveRDS(model_data_black, paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/model_data_race', race_option, '_sex', sex_option, '.RDS'))

  }
}

## Plot nmx
pdf(paste0(output_dir, '/fbp_kt.pdf'), height=6, width=9)
v <- 'fbp_kt'
m <- make_county_map(map_dt = model_data_black[year==2010,],
                     map_sp = counties,
                     map_var = v,
                     legend_title = v,
                     high_is_good = TRUE,
                     map_title = v,
                     map_limits = c(0,0.25)) + guides(fill=guide_legend(title="Foreign-\nborn, %"))
print(m)
dev.off()

# Create a neighborhood file
nb.FOQ <- poly2nb(counties, queen=TRUE)
# Create an INLA weight matrix
lw.FOQ <- nb2INLA("FOQ_INLA",nb.FOQ)
an <- data.frame(1:length(counties),counties$fips)
o <- match(model_data_black$fips,an$counties.fips)
ID <- an$X1.length.counties.[o]
model_data_black[, ID := ID]

## Series of Bayesian logistic models, starting with just metro and AR1 and adding covariates.
## 1. Metro, AR1 on year
## 2. Metro + Regions, AR1 on year
## 3. Metro + Regions + Covariates, AR1 on year
## 4. Metro + Regions + Covariates, AR1 on year + Besag on county
if(race_option == 0) {
  inla_formula1 <- deaths ~ as.factor(metroname) + f(year, model="ar1")
  inla_formula2 <- deaths ~ as.factor(metroname) + as.factor(regionname) + f(year, model="ar1")  
  inla_formula3 <- deaths ~ 
    twhk_pv + twhk_is + twhk_ce + fwp_kt + as.factor(metroname) + as.factor(regionname) + 
    Obese + Smokers + Excessive.Drinking + PCP.Z + MHP.Z + Uninsured + Unemployed + Severe.Housing.Problems +
    f(year, model="ar1")
  inla_formula4 <- deaths ~ 
    twhk_pv + twhk_is + twhk_ce + fwp_kt + as.factor(metroname) + as.factor(regionname) +
    Obese + Smokers + Excessive.Drinking + PCP.Z + MHP.Z + Uninsured + Unemployed + Severe.Housing.Problems +
    f(year, model="ar1") + f(ID,model="besag",graph="FOQ_INLA")
}
if(race_option == 1) {
  inla_formula1 <- deaths ~ as.factor(metroname) + f(year, model="ar1")
  inla_formula2 <- deaths ~ as.factor(metroname) + as.factor(regionname) + f(year, model="ar1")  
  inla_formula3 <- deaths ~ 
    tbhk_pv + tbhk_is + tbhk_ce + fbp_kt + as.factor(metroname) + as.factor(regionname) + 
    Obese + Smokers + Excessive.Drinking + PCP.Z + MHP.Z + Uninsured + Unemployed + Severe.Housing.Problems +
    f(year, model="ar1")
  inla_formula4 <- deaths ~ 
    tbhk_pv + tbhk_is + tbhk_ce + fbp_kt + as.factor(metroname) + as.factor(regionname) + 
    Obese + Smokers + Excessive.Drinking + PCP.Z + MHP.Z + Uninsured + Unemployed + Severe.Housing.Problems + 
    f(year, model="ar1") + f(ID,model="besag",graph="FOQ_INLA")
}
if(year_as_factor== TRUE) {
  if(race_option == 0) {
    inla_formula1 <- deaths ~ as.factor(metroname) + as.factor(year)
    inla_formula2 <- deaths ~ as.factor(metroname) + as.factor(regionname) + as.factor(year)  
    inla_formula3 <- deaths ~ 
      twhk_pv + twhk_is + twhk_ce + fwp_kt + as.factor(metroname) + as.factor(regionname) + 
      #Obese + Smokers + Excessive.Drinking + PCP.Z + MHP.Z + Uninsured + Unemployed + Severe.Housing.Problems +
      as.factor(year)
    inla_formula4 <- deaths ~ 
      twhk_pv + twhk_is + twhk_ce + fwp_kt + as.factor(metroname) + as.factor(regionname) + 
      #Obese + Smokers + Excessive.Drinking + PCP.Z + MHP.Z + Uninsured + Unemployed + Severe.Housing.Problems + 
      as.factor(year) + f(ID,model="besag",graph="FOQ_INLA")
  }
  if(race_option == 1) {
    inla_formula1 <- deaths ~ as.factor(metroname) + as.factor(year)
    inla_formula2 <- deaths ~ as.factor(metroname) + as.factor(regionname) + as.factor(year)
    inla_formula3 <- deaths ~ 
      tbhk_pv + tbhk_is + tbhk_ce + fbp_kt + as.factor(metroname) + as.factor(regionname) + 
      #Obese + Smokers + Excessive.Drinking + PCP.Z + MHP.Z + Uninsured + Unemployed + Severe.Housing.Problems +
      as.factor(year)
    inla_formula4 <- deaths ~ 
      tbhk_pv + tbhk_is + tbhk_ce + fbp_kt + as.factor(metroname) + as.factor(regionname) + 
      #Obese + Smokers + Excessive.Drinking + PCP.Z + MHP.Z + Uninsured + Unemployed + Severe.Housing.Problems + 
      as.factor(year) + f(ID,model="besag",graph="FOQ_INLA")
  }
}
if(!is.null(region_subset)) {
  model_data_black <- model_data_black[regionname %in% region_subset, ]
  model_data_black[, regionname := factor(regionname, levels = c(region_subset[1], unique(model_data_black[!(regionname %in% region_subset[1]), regionname])))]
}
if(is.null(region_subset)) {
  model_data_black[, regionname := factor(regionname, levels = c('Pacific', unique(model_data_black[!(regionname %in% 'Pacific'), regionname])))]
}
inla_model_1 = inla(inla_formula1,
                    family = "binomial",
                    data = model_data_black,
                    Ntrials = model_data_black[, total_pop],
                    verbose = FALSE,
                    control.compute=list(config = TRUE, dic = TRUE),
                    control.inla=list(int.strategy='eb'))
model_data_black[, inla_pred_1 := inla_model_1$summary.fitted.values$mean]
model_data_black[, inla_residual_1 := nmx - inla_pred_1]
inla_model_2 = inla(inla_formula2,
                    family = "binomial",
                    data = model_data_black,
                    Ntrials = model_data_black[, total_pop],
                    verbose = FALSE,
                    control.compute=list(config = TRUE, dic = TRUE),
                    control.inla=list(int.strategy='eb'))
model_data_black[, inla_pred_2 := inla_model_2$summary.fitted.values$mean]
model_data_black[, inla_residual_2 := nmx - inla_pred_2]
inla_model_3 = inla(inla_formula3,
                    family = "binomial",
                    data = model_data_black,
                    Ntrials = model_data_black[, total_pop],
                    verbose = FALSE,
                    control.compute=list(config = TRUE, dic = TRUE),
                    control.inla=list(int.strategy='eb'))
model_data_black[, inla_pred_3 := inla_model_3$summary.fitted.values$mean]
model_data_black[, inla_residual_3 := nmx - inla_pred_3]
inla_model_4 = inla(inla_formula4,
                    family = "binomial",
                    data = model_data_black,
                    Ntrials = model_data_black[, total_pop],
                    verbose = FALSE,
                    control.compute=list(config = TRUE, dic = TRUE),
                    control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
                    control.fixed=list(prec.intercept = 0,
                                       prec = 1))
model_data_black[, inla_pred_4 := inla_model_4$summary.fitted.values$mean]
model_data_black[, inla_residual_4 := nmx - inla_pred_4]

## Make tables comparing coefficients
all_tables <- lapply(1:4, make_coef_table)
all_coefs <- rbindlist(all_tables, fill=TRUE)
all_dic_morans <- data.table(model = rep(paste0('Model ', 1:4),3),
                             name = c(rep("Global Moran's I",4), rep("DIC",4), rep("RMSE",4)),
                             coef = rep(.00001,12))
for(m in 1:4) {
  lisa <- plot_lisa(lisa_var = paste0('inla_residual_', m),
                    lisa_dt = model_data_black[year==2010,],
                    lisa_sp = counties,
                    lisa_id = 'fips',
                    lisa_var_name = 'Residuals',
                    sig = 0.05)
  if(m==3) {
    png(paste0(plot_dir, '/inla_',m,'_sex', sex_option, '_residual_lisa.png'), width = 1200, height = 1400, res = 120)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(14, 12)))
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    print(lisa[[1]] + theme(legend.position="none"), vp = vplayout(1:8, 1:12))
    print(lisa[[2]] + theme(legend.position="none"), vp = vplayout(9:14, 1:12))
    dev.off()
  }
  
  all_dic_morans[model == paste0('Model ', m) & name == "Global Moran's I", coef := lisa[[3]]]
  all_dic_morans[model == paste0('Model ', m) & name == "DIC", coef := ifelse(is.nan(get(paste0('inla_model_',m))$dic$dic) | is.infinite(get(paste0('inla_model_',m))$dic$dic), 
                                                                              get(paste0('inla_model_',m))$dic$deviance.mean,
                                                                              get(paste0('inla_model_',m))$dic$dic)]
  all_dic_morans[model == paste0('Model ', m) & name == 'RMSE', coef := model_data_black[, sqrt(weighted.mean(get(paste0('inla_residual_', m))^2, w = total_pop))]]
}
all_coefs <- rbind(all_coefs, all_dic_morans, fill = TRUE)
saveRDS(all_coefs, file = paste0(output_dir, '/sex', sex_option, '_coef_table_just_inla.RDS'))

lisa <- plot_lisa(lisa_var = paste0('nmx'),
                  lisa_dt = model_data_black[year==2010,],
                  lisa_sp = counties,
                  lisa_id = 'fips',
                  lisa_var_name = 'nmx',
                  sig = 0.05,
                  matrix = 'nn',
                  n_neighbors = 10)
png(paste0(plot_dir, '/sex_', sex_option, 'nmx_lisa.png'), width = 1200, height = 1400, res = 120)
grid.newpage()
pushViewport(viewport(layout = grid.layout(14, 12)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(lisa[[1]] + theme(legend.position="none"), vp = vplayout(1:8, 1:12))
print(lisa[[2]] + theme(legend.position="none"), vp = vplayout(9:14, 1:12))
dev.off()

saveRDS(lisa[[4]], paste0(plot_dir, '/sex_', sex_option, '_lisa_counties.rds'))

  }
}
