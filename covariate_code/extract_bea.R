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

## Load data and reshape to be by FIPS + year.
repo <- 'c:/Users/ngraetz/Documents/repos/rwjf_counties/'
d <- fread("C:/Users/ngraetz/Downloads/CA4/CA4_1969_2016__ALL_AREAS.csv", header = TRUE)
year_cols <- paste0('year',1969:2016)
d[, fips := as.character(GeoFIPS)]
d[nchar(fips)==4, fips := paste0('0',fips)]
for(v in as.character(1969:2016)) d[, (v) := as.numeric(get(v))]

per <- fread("C:/Users/ngraetz/Downloads/CA35/CA35_1969_2016__ALL_AREAS.csv", header = TRUE)
per[, fips := as.character(GeoFIPS)]
per[nchar(fips)==4, fips := paste0('0',fips)]
for(v in as.character(1969:2016)) per[, (v) := as.numeric(get(v))]

## Need to fix Virginia FIPS codes.
# template <- as.data.table(counties@data)
# write.csv(template, 'C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/county_template.csv', row.names=FALSE)
# d_codes <- unique(d[, c('fips','GeoName')])
# write.csv(va_codes, 'C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/bea_template.csv', row.names=FALSE)
# 
# counties <- readOGR("C:/Users/ngraetz/Downloads/cb_2016_us_county_20m", 'cb_2016_us_county_20m')
# counties@data <- transform(counties@data, fips = paste0(STATEFP, COUNTYFP)) # create unique county 5-digit fips
# ## Drop Alaska, Hawaii, Puerto Rico
# counties <- counties[counties$STATEFP != '02' &
#                        counties$STATEFP != '15' &
#                        counties$STATEFP != '72',]
# all <- counties$fips
# all_missing_fips <- all[!(all %in% d[, fips])]
# for(f in all_missing_fips) {
#   this_name <- template[fips==f, NAME]
#   this_state <- substr(f,1,2)
#   message(paste0(f, ' missing: ', this_name))
#   this_state <- d_codes[grep(paste0('^',this_state), fips),]
#   possible_matches <- this_state[grep(this_name, GeoName), c('fips','GeoName')]
#   for(m in possible_matches[, fips]) message(paste0('          ',m,' in BEA: ', possible_matches[fips==m, GeoName]))
# }

# merged_fips <- fread('C:/Users/ngraetz/Documents/repos/rwjf/bea_county_template.csv')
# merged_fips[, fips := as.character(fips)]
# merged_fips[nchar(fips)==4, fips := paste0('0',fips)]
# merged_fips[, merged_fips := as.character(merged_fips)]
# merged_fips[nchar(merged_fips)==4, merged_fips := paste0('0',merged_fips)]
# for(f in merged_fips[, merged_fips]) d[fips == f, fips := merged_fips[merged_fips==f, fips]]
# d <- d[, lapply(.SD, sum, na.rm=TRUE), by=c('fips','Description'), .SDcols=as.character(1969:2016)] 

d <- melt(d, id.vars = c('fips','Description'), measure.vars = as.character(1969:2016), variable.name = 'year')
d[, year := as.numeric(as.character(year))]
d[, value := as.numeric(value)]
## This one is doubled under two headings, but since we aren't using it for anything just drop.
d <- d[Description != 'Employer contributions for government social insurance', ]
d <- dcast(d, fips + year ~ Description, value.var = 'value')

## Calculate indicators of interest.
## Total personal income = (Earnings by place of work - contributions for government social insurance + adjustment for residence) + dividends/interest/rent + transfers.
d[, percent_transfers := `Plus: Personal current transfer receipts` / `Personal income (thousands of dollars)`]

d[, percent_wage_salary_employment := `Wage and salary employment` / `Total employment`]

d[, income_per_capita := `Per capita personal income (dollars) 4/`]

d[, total_employees := `Total employment`]

## Save.
d <- d[, c('fips','year','percent_transfers','percent_wage_salary_employment','income_per_capita','total_employees')]
write.csv(d, paste0(repo, 'covariate_clean_data/bea_covs.csv'), row.names=FALSE)

d <- fread('C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/bea_covs.csv')
#d <- dcast(d, fips ~ year, value.var = 'percent_transfers')
#d[, transfer_change := `2016` - `1990`]
v <- 'percent_transfers'
source("C:/Users/ngraetz/Documents/repos/spatial_demography_2018/functions.R")
counties <- readOGR("C:/Users/ngraetz/Downloads/cb_2016_us_county_20m", 'cb_2016_us_county_20m')
counties@data <- transform(counties@data, fips = paste0(STATEFP, COUNTYFP)) # create unique county 5-digit fips
counties <- counties[counties@data$fips != "02016", ] # Drop Aleutians West, AK - screws up plots 
counties <- counties[counties$STATEFP=='51', ]
background.dt <- as.data.table(fortify(counties, region = 'STATEFP'))
map_colors <- c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')
counties$state <- as.numeric(counties$STATEFP)
states <- gUnaryUnion(counties, id = counties@data$state)
background.dt.states <- as.data.table(fortify(states, region = 'state'))

mort <- fread('C:/Users/ngraetz/Documents/repos/rwjf/bea_county_template.csv')
mort[is.na(bea_fips), bea_fips := fips]
mort[, year := 1990]
mort[, bea_fips := as.character(bea_fips)]
mort[nchar(bea_fips)==4, bea_fips := paste0('0',bea_fips)]
d[, bea_fips := fips]
d[, fips := NULL]
mort <- merge(mort, d, by=c('year','bea_fips'))

pdf('C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/transfers.pdf', height=12, width=18)
m <- make_county_map(map_dt = mort[year==1990, ],
                     map_sp = counties,
                     map_var = v,
                     legend_title = v,
                     high_is_good = FALSE,
                     map_title = '1990 Transfers, %',
                     map_limits = c(0,0.4)) + guides(fill=guide_legend(title="Transfers, %"))
print(m)
m <- make_county_map(map_dt = mort[year==2015, ],
                     map_sp = counties,
                     map_var = v,
                     legend_title = v,
                     high_is_good = FALSE,
                     map_title = '2015 Transfers, %',
                     map_limits = c(0,0.4)) + guides(fill=guide_legend(title="Transfers, %"))
print(m)
m <- make_county_map(map_dt = model_data[year==1990, ],
                     map_sp = counties,
                     map_var = 'nmx',
                     legend_title = 'nmx',
                     high_is_good = FALSE,
                     map_title = '1990 Transfers, %',
                     map_limits = c(0,0.015)) + guides(fill=guide_legend(title="Transfers, %"))
print(m)
dev.off()

