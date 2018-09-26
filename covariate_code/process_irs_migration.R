## Process in- and out-migration county flows from IRS.
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

## Load pops
all_pops <- fread("C:/Users/ngraetz/Downloads/all_pops_fixed.csv")
all_pops <- all_pops[, list(total_pop=sum(total_pop)), by=c('year','fips')]

source('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/functions.R')
pdf('C:/Users/ngraetz/Documents/Penn/papers/rwjf/irs_migration/mig_rates.pdf', height=6, width=9)
pull_migration <- function(flow_type, y) {
out <- fread(paste0("C:/Users/ngraetz/Downloads/county1011/county", flow_type, "flow1011.csv"))
cols <- c('State_Code_Origin','State_Code_Dest','County_Code_Origin','County_Code_Dest')
out[ , (cols) := lapply(.SD, as.character), .SDcols = cols]
out[nchar(State_Code_Origin)==1, State_Code_Origin := paste0('0', State_Code_Origin)]
out[nchar(State_Code_Dest)==1, State_Code_Dest := paste0('0', State_Code_Dest)]
out[nchar(County_Code_Origin)==1, County_Code_Origin := paste0('00', County_Code_Origin)]
out[nchar(County_Code_Origin)==2, County_Code_Origin := paste0('0', County_Code_Origin)]
out[nchar(County_Code_Dest)==1, County_Code_Dest := paste0('00', County_Code_Dest)]
out[nchar(County_Code_Dest)==2, County_Code_Dest := paste0('0', County_Code_Dest)]
out[, fips_origin := paste0(State_Code_Origin, County_Code_Origin)]
out[, fips_dest := paste0(State_Code_Dest, County_Code_Dest)]
if(flow_type=='out') out[, fips := fips_origin]
if(flow_type=='in') out[, fips := fips_dest]
out[, year := y]
out <- out[grepl('Tot Mig-US & For|Non-Migrants', County_Name), ]
out[grep('Tot Mig-US & For', County_Name), type := 'migrant']
out[grep('Non-Migrants', County_Name), type := 'non_migrant']
out[, total := Exmpt_Num]
out <- dcast(out, fips + year ~ type, value.var = 'total')
out[, total := migrant + non_migrant]
out[, rate := migrant / total]
out[, log_migrant := log(migrant)]
return(out) 
}

in_mig <- pull_migration(flow_type='in',2010)
setnames(in_mig, c('rate','migrant'), c('in_rate','in_migrants'))
out_mig <- pull_migration(flow_type='out',2010)
setnames(out_mig, c('rate','migrant'), c('out_rate','out_migrants'))
net <- merge(in_mig[, c('year','fips','in_rate','in_migrants')], out_mig[, c('year','fips','out_rate','out_migrants')], by=c('year','fips'))
net[, net := in_rate - out_rate]

ggplot() + 
  geom_point(data=net,
             aes(x=in_migrants,
                 y=out_migrants)) +
  geom_abline(slope=1, intercept=0) + 
  labs(x='In-migration return total', y='Out-migration return total') +
  theme_minimal()

m <- make_county_map(map_dt = net,
                     map_sp = counties,
                     map_var = 'net',
                     legend_title = 'Rate',
                     high_is_good = TRUE,
                     map_title = 'Net migration rate, 2010',
                     map_limits = c(-0.02, 0.02)) + guides(fill=guide_legend(title="Rate"))
print(m)

m <- make_county_map(map_dt = out,
                     map_sp = counties,
                     map_var = 'rate',
                     legend_title = 'Rate',
                     high_is_good = TRUE,
                     map_title = ifelse(flow_type=='out', paste0('Out-migration, ', y), paste0('In-migration, ', y)),
                     map_limits = c(0,0.08)) + guides(fill=guide_legend(title="Rate"))
print(m)
m <- make_county_map(map_dt = out,
                     map_sp = counties,
                     map_var = 'log_migrant',
                     legend_title = "ln(migrants)",
                     high_is_good = TRUE,
                     map_title = ifelse(flow_type=='out', paste0('Out-migration, ', y), paste0('In-migration, ', y)),
                     map_limits = c(5,12)) + guides(fill=guide_legend(title="ln(migrants)"))
print(m)

dev.off()



mig_2010 <- read.delim('https://www2.census.gov/programs-surveys/demo/tables/geographic-mobility/2009/county-to-county-migration-2005-2009/ctyxcty_us.txt',
                header=FALSE,
                fill=TRUE,
                sep='')
mig_2010 <- as.data.table(mig_2010)



## Compile Census in, out, and net migration for 2000-2015 (migrants identified based on current residence and residence last year).
repo <- 'C:/Users/ngraetz/Documents/repos/rwjf_counties/'
## 2000
in_mig_2000 <- read.table('https://www2.census.gov/programs-surveys/demo/tables/geographic-mobility/2000/county-to-county-flows/intxt_flow.txt',fill=TRUE,header=FALSE)
in_mig_2000 <- as.data.table(in_mig_2000)
setnames(in_mig_2000, c('V1','V3'), c('fips','in_mig'))
in_mig_2000 <- in_mig_2000[, list(in_mig=sum(in_mig)), by='fips']
out_mig_2000 <- read.table('https://www2.census.gov/programs-surveys/demo/tables/geographic-mobility/2000/county-to-county-flows/outtxt_flow.txt',fill=TRUE,header=FALSE)
out_mig_2000 <- as.data.table(out_mig_2000)
setnames(out_mig_2000, c('V1','V3'), c('fips','out_mig'))
out_mig_2000 <- out_mig_2000[, list(out_mig=sum(out_mig)), by='fips']
mig_2000 <- merge(in_mig_2000, out_mig_2000, by='fips')
mig_2000[, net_mig := in_mig - out_mig]
mig_2000[, year := 2000]
mig_2000[, fips := as.character(fips)]
mig_2000[nchar(fips)==4, fips := paste0('0',fips)]
## 2010
mig_2010_fwf <- read.fwf('C:/Users/ngraetz/Downloads/2010_ctyxcty_us.txt', 388)
mig_2010 <- as.data.table(mig_2010_fwf)
mig_2010[, V1 := as.character(V1)]
mig_2010[, mig := substr(V1, 373, 380)]
mig_2010[, mig := as.numeric(gsub(' ','',mig))]
mig_2010[, dest_fips := substr(V1, 2, 6)]
mig_2010[, orig_fips := substr(V1, 8, 12)]
in_mig_2010 <- mig_2010[, list(in_mig=sum(mig)), by='dest_fips']
out_mig_2010 <- mig_2010[, list(out_mig=sum(mig)), by='orig_fips']
setnames(in_mig_2010, 'dest_fips', 'fips')
setnames(out_mig_2010, 'orig_fips', 'fips')
mig_2010 <- merge(in_mig_2010, out_mig_2010, by='fips')
mig_2010[, net_mig := in_mig - out_mig]
mig_2010[, year := 2010]
## 2015
mig_2015_fwf <- read.fwf('C:/Users/ngraetz/Downloads/2015_CtyxCty_US.txt', 388)
mig_2015 <- as.data.table(mig_2015_fwf)
mig_2015[, V1 := as.character(V1)]
mig_2015[, mig := substr(V1, 373, 380)]
mig_2015[, mig := as.numeric(gsub(' ','',mig))]
mig_2015[, dest_fips := substr(V1, 2, 6)]
mig_2015[, orig_fips := substr(V1, 8, 12)]
in_mig_2015 <- mig_2015[, list(in_mig=sum(mig)), by='dest_fips']
out_mig_2015 <- mig_2015[, list(out_mig=sum(mig)), by='orig_fips']
setnames(in_mig_2015, 'dest_fips', 'fips')
setnames(out_mig_2015, 'orig_fips', 'fips')
mig_2015 <- merge(in_mig_2015, out_mig_2015, by='fips')
mig_2015[, net_mig := in_mig - out_mig]
mig_2015[, year := 2015]
## All migration numbers by fips/year (2000 Census, 2010 pooled ACS, 2015 pooled ACS).
all_mig <- rbind(mig_2000, mig_2010, mig_2015)
write.csv(all_mig, paste0(repo, 'covariate_clean_data/census_acs_migration.csv'), row.names=FALSE)

