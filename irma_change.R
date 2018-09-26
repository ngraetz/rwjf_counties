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
library(kableExtra)
library(Hmisc)
library(relaimpo)
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
covs <- c('percent_transfers','mds_pc','obesity','college','net_in_mig','manufacturing','poverty_all')
year_range <- c(1990,2000,2010,2015)
plot_trends <- FALSE

start_year <- min(year_range)
end_year <- max(year_range)
sex_name <- ifelse(sex_option==1,'Male','Female')
output_name <- paste0(capitalize(domain), ' ', toupper(race), ' ', capitalize(sex_name))

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

## Hack on covariates with weird year ranges.
obesity <- fread("C:/Users/ngraetz/Documents/repos/rwjf_counties/covariate_clean_data/chr_covs.csv")
obesity[year==2004, year := 1990]
obesity[year==2013, year := 2015]
mort <- merge(mort, obesity, by=c('fips','year'), all.x=TRUE)

mig <- fread("C:/Users/ngraetz/Downloads/Working Age Net In-Migration Rates by Region and Metro Category, Race, Sex, and Period.csv")
setnames(mig, c('urban.x','total in-migration rate'), c('urban','net_in_mig'))
mig <- mig[race=='White' & sex == sex_option, ]
mig[, race := '0']
mig[, race := as.numeric(race)]
mig[year==1991, year := 1990]

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

## Collapse to metro-region
mort[regionname=='Nonmetro', metro := 4]
mort[, urban := metro]
mort <- mort[, lapply(.SD, weighted.mean, w=total_pop, na.rm=TRUE), .SDcols=covs[covs!='net_in_mig'], by=c('region','urban','year','race','sex')]

## Merge net in-migration indirect estimates from Arun by region/metro.
mort <- merge(mort, mig, by=c('region','urban','year','race','sex'))

## Merge life expectancy at birth.
lt <- fread(paste0(repo, 'Annual Life Tables by Metro 4-Category and Region by Sex and Race.csv'))
lt <- lt[Age.x==0, c('year','urban','region','urban_name','region_name','race','sex','ex')]
lt[race=='White', race := '0']
lt[race=='Black', race := '1']
lt[, race := as.numeric(race)]
change <- merge(mort, lt, by=c('year','urban','region','race','sex'))

## Reshape and calculate change.
change <- change[year %in% c(min(year_range), max(year_range))]
change <- dcast(change, region_name + urban_name + race + sex ~ year, value.var = c('ex',covs))
for(c in c('ex',covs)) change[, (paste0('change_',c)) := get(paste0(c,'_',max(year_range))) - get(paste0(c,'_',min(year_range)))]
for(c in c('ex',covs)) change[, (paste0('change1_',c)) := round(get(paste0('change_',c)), 1)]
for(c in c('ex',covs)) change[, (paste0('change2_',c)) := linebreak(paste0(as.character(round(get(paste0('change_',c)), 1)), '\n(', as.character(round(get(paste0(c,'_',min(year_range))),1)),
                              '-', as.character(round(get(paste0(c,'_',max(year_range))),1)),')'))]

## Reshape wide by metro for change table.
table_1 <- dcast(change, region_name ~ urban_name, value.var = paste0('change1_',c('ex',covs)))
table_2 <- dcast(change, region_name ~ urban_name, value.var = paste0('change2_',c('ex',covs)))
write.csv(table_1, paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/results/irma_pdr/change1_', sex_name, '.csv'), row.names = FALSE)
write.csv(table_2, paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/results/irma_pdr/change2_', sex_name, '.csv'), row.names = FALSE)

## Correlation matrices
write.csv(cor(change[, c('change_ex','change_percent_transfers','change_mds_pc',"change_obesity","change_college","change_net_in_mig","change_manufacturing")]),
          paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/results/irma_pdr/cor_matrix_', sex_option,'.csv'), row.names=TRUE)

## Make scatters of change in LE and change in each covariate, 2000-2015. Facet by region, color by metro.
#   Nonmetro: #de2d26
#   Small metro: #bdbdbd
#   Suburbs: #2ca25f
#   Metro: #253494
colors <- cbind(col2rgb('#de2d26'), col2rgb('#bdbdbd'), col2rgb('#2ca25f'), col2rgb('#253494'))
colnames(colors) <- c('Nonmetro','Small/Medium Metro','Suburbs','Metro')
cov_names <- get_clean_cov_names()
change[urban_name=='Small Metro', urban_name := 'Small/Medium Metro']

pdf(paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/results/irma_pdr/change_scatters_', sex_name, '.pdf'), height=9, width=12)
change[, urban_name := factor(urban_name, levels=c('Large Central Metro','Large Metro Suburb','Small/Medium Metro','Nonmetro'))]
scatter_data <- melt(change, id.vars = c('region_name','urban_name','change_ex'), measure.vars = paste0('change_',covs), value.name = 'cov_value', variable.name = 'cov')
scatter_data[, cov := gsub('change_','',cov)]
scatter_data <- merge(scatter_data, cov_names, by.x='cov', by.y='fe')
scatter_data[cov_name=='Obesity', cov_name := 'Change in percent obese, 2004-2013']
scatter_data[cov_name=='College', cov_name := 'Change in percent college graduates']
scatter_data[cov_name=='MDs/pc', cov_name := 'Change in MDs per thousand']
scatter_data[cov_name=='Net In-mig', cov_name := 'Change in cumulative rate of net in-migration, ages 22.5-62.5']
scatter_data[cov_name=='Transfers', cov_name := 'Change in transfers as a percent share of personal income']
scatter_data[cov_name=='Manufacturing', cov_name := 'Change in percent of labor force in manufacturing']
write.csv(scatter_data, paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/results/irma_pdr/change_scatters_', sex_name, '.csv'), row.names=FALSE)
for(c in unique(scatter_data[, cov_name])) {
this_gg <- ggplot(data = scatter_data[cov_name == c]) +
  geom_point(aes(x = cov_value,
                 y = change_ex,
                 fill = urban_name),
             size = 8,
             shape = 21) +
  scale_fill_manual(values=c('#253494','#2ca25f','#bdbdbd','#de2d26')) +
  guides(fill=guide_legend(title="Metro")) +
  theme_minimal() + 
  ylab('Change in life expectancy') + 
  xlab(c) + 
  theme(panel.spacing = unit(2, "lines")) + 
  theme(panel.border = element_rect(color="black", fill=NA)) +
  theme(legend.position="bottom")
print(this_gg)
}
dev.off()

## Decompose R^2 from OLS predicting change in LE with change in each covariate.
ols_f <- paste0('change_ex ~ ', paste(c(paste0('change_',covs), 'as.factor(region_name)', 'as.factor(urban_name)'), collapse = ' + '))
ols_change <- glm(ols_f, data=change)
decomp <- calc.relimp(ols_change, type="lmg", rela=F)@lmg
decomp <- data.table(variable=names(decomp), variance_explained=decomp*100)
decomp <- rbind(decomp, data.table(variable='unexplained',variance_explained=100-sum(decomp[, variance_explained])))
write.csv(decomp, paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/results/irma_pdr/decomp_', sex_name, '.csv'), row.names = FALSE)



# library(data.table)
# library(ggplot2)
# data_dir <- 'your directory'
# out_dir <- 'your directory'
# sex_name <- 'Male'
# scatter_data <- fread(paste0(data_dir, 'change_scatters_', sex_name, '.csv'))
# pdf(paste0(out_dir, 'change_scatters_', sex_name, '.pdf'), height=9, width=12)
# for(c in unique(scatter_data[, cov_name])) {
#   this_gg <- ggplot(data = scatter_data[cov_name == c]) +
#     geom_point(aes(x = cov_value,
#                    y = change_ex,
#                    fill = urban_name),
#                size = 8,
#                shape = 21) +
#     scale_fill_manual(values=c('#253494','#2ca25f','#bdbdbd','#de2d26'), name = '') +
#     theme_minimal() + 
#     labs(y='Change in life expectancy', x=c) +  
#     theme(panel.spacing = unit(2, "lines")) + 
#     theme(panel.border = element_rect(color="black", fill=NA)) +
#     theme(legend.position="bottom")
#   print(this_gg)
# }
# dev.off()

