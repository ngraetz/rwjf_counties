library(data.table)
library(ggplot2)
library(haven)
library(rgdal)
library(readxl)
repo <- 'C:/Users/ngraetz/Documents/repos/rwjf_counties/'
input_dir <- paste0(repo, 'nmx_clean_data/')
cov_dir <- paste0(repo, 'covariate_clean_data/')

## Get educational attainment estimates from 1990 Census, 2000 Census, 2005 ACS, 2010 ACS, 2015 ACS
raw_dir <- paste0(repo, 'covariate_raw/factfinder_edu_2000_2015/')
cb <- fread(paste0(repo, 'covariate_code/factfinder_meta.csv'))
pull_edu_data <- function(r) {
meta <- cb[r,]
d <- fread(paste0(repo, 'covariate_raw/factfinder_edu_2000_2015/', meta[, file], '.csv'))
for(v in c("less_12","hs","assoc","college","total","fips")) {
  if(grepl('[+]', meta[, get(v)])) d[, (v) := get(strsplit(meta[, get(v)], '[+]')[[1]][1]) + get(strsplit(meta[, get(v)], '[+]')[[1]][2])]
  if(!grepl('[+]', meta[, get(v)])) setnames(d, meta[, get(v)], v)
}
for(v in c("less_12","hs","assoc","college")) d[, (v) := get(v) / total]
d <- d[, c("less_12","hs","assoc","college","total","fips")]
d[, year := meta[, year]]
d[, race := meta[, race]]
d[, sex := meta[, sex]]
return(d)
}
all_edu <- rbindlist(lapply(min(cb[var=='edu', row]):max(cb[var=='edu', row]), pull_edu_data))
all_edu[, total := NULL]

## Add 1990 race-specific/all-sex data from 1990 Census SF3 website. Collapse other years across sex.
pop_white <- readRDS(paste0(input_dir, 'nhw_total_pop.RDS'))
pop_black <- readRDS(paste0(input_dir, 'nhb_total_pop.RDS'))
pop <- rbind(pop_white, pop_black)

all_edu[, fips := as.character(fips)]
all_edu[nchar(fips)==4, fips := paste0('0',fips)]
all_edu <- merge(all_edu, pop, all.x=TRUE)
all_edu <- all_edu[, lapply(.SD, weighted.mean, w=total_pop, na.rm=TRUE), .SDcols=c('less_12','college'), by=c('fips','year','race')]

census_1990 <- fread(paste0(repo, 'covariate_clean_data/census_1990.csv'))
census_1990[, less_12 := less_12 / total_persons]
census_1990[, college := college / total_persons]
census_1990[, race := 0]
census_1990 <- census_1990[nchar(fips)==5, ]
census_1990 <- census_1990[, c('fips','year','race','less_12','college')]
## Some counties are included in multiple files on the Census 1990 SF3 website (it chunks out counties so no one file is too huge).
## They are virtually identical, so just take mean.
census_1990 <- census_1990[, lapply(.SD, mean, na.rm=TRUE), .SDcols=c('less_12','college'), by=c('fips','year','race')]

all_edu <- rbind(all_edu, census_1990)

## Get % foreign-born by race estimates from 2000 Census, 2005 ACS, 2010 ACS, 2015 ACS.
pull_fb_data <- function(r) {
  meta <- cb[r,]
  d <- fread(paste0(repo, 'covariate_raw/factfinder_fb_2015/', meta[, file], '.csv'))
  for(v in c("fb_total","fb","fips")) {
    if(grepl('[+]', meta[, get(v)])) d[, (v) := get(strsplit(meta[, get(v)], '[+]')[[1]][1]) + get(strsplit(meta[, get(v)], '[+]')[[1]][2])]
    if(!grepl('[+]', meta[, get(v)])) setnames(d, meta[, get(v)], v)
  }
  d <- d[, c("fb","fb_total","fips")]
  d[, year := meta[, year]]
  d[, race := meta[, race]]
  d[, sex := meta[, sex]]
  return(d)
}
all_fb <- rbindlist(lapply(min(cb[var=='fb', row]):max(cb[var=='fb', row]), pull_fb_data))
## Collapse 2015 across sex because 1990, 2000, 2010 isn't sex-specific.
all_fb <- all_fb[, list(fb=sum(fb), fb_total=sum(fb_total)), by=c('fips','year','race')]
all_fb[, fb := fb/fb_total]
all_fb <- all_fb[, c('fips','year','race','fb')]
all_fb[, fips := as.character(fips)]
all_fb[nchar(fips)==4, fips := paste0('0',fips)]
## Get the rest from the Brown database (FactFinder doesn't have 1990, 2000, 2010 by race).
brown_fb <- fread(paste0(repo, 'covariate_raw/exposures_by_nativity.csv'))
brown_fb[, nhw_fb_1990 := w9k_for/w9k_tot]
brown_fb[, nhw_fb_2000 := w0k_for/w0k_tot]
brown_fb[, nhw_fb_2010 := w1k_for/w1k_tot]
brown_fb[, nhb_fb_1990 := b9k_for/b9k_tot]
brown_fb[, nhb_fb_2000 := b0k_for/b0k_tot]
brown_fb[, nhb_fb_2010 := b1k_for/b1k_tot]
brown_fb <- melt(brown_fb, id.vars = 'fips', measure.vars = names(brown_fb)[grep('_fb_', names(brown_fb))], value.name = 'fb')
brown_fb[grep('nhw', variable), race := 0]
brown_fb[grep('nhb', variable), race := 1]
brown_fb[, year := gsub('nhw_fb_','',variable)]
brown_fb[, year := gsub('nhb_fb_','',year)]
brown_fb <- brown_fb[, c('fips','year','race','fb')]
all_fb <- rbind(brown_fb, all_fb)

## Save 
all_edu[, fips := as.character(fips)]
all_edu[nchar(fips)==4, fips := paste0('0',fips)]
write.csv(all_edu, paste0(repo, 'covariate_clean_data/factfinder_edu.csv'), row.names = FALSE)
all_fb[, fips := as.character(fips)]
all_fb[nchar(fips)==4, fips := paste0('0',fips)]
write.csv(all_fb, paste0(repo, 'covariate_clean_data/factfinder_fb.csv'), row.names = FALSE)



## Get 1990 % manufacturing from 1990 Census and 2015 from pooled FactFinder ACS.
## fips = GEO.id2
## manufacturing = HC01_EST_VC06
## total_persons = HC01_EST_VC01
manu_2015 <- fread(paste0(repo, 'covariate_raw/factfinder_manufacturing_2015/ACS_15_5YR_S2403_with_ann_no_labels.csv'))
setnames(manu_2015, c('GEO.id2','HC01_EST_VC06','HC01_EST_VC01'), c('fips','manufacturing','total_persons'))
manu_2015[, manufacturing := manufacturing / total_persons]
manu_2015 <- manu_2015[, c('fips','manufacturing')]
manu_2015[, fips := as.character(fips)]
manu_2015[nchar(fips)==4, fips := paste0('0',fips)]
manu_2015[, year := 2015]

manu_1990 <- fread(paste0(repo, 'covariate_clean_data/manufacturing_census_1990.csv'))
manu_1990[, manufacturing := manufacturing / total_persons]
manu_1990 <- manu_1990[nchar(fips)==5, ]
manu_1990 <- manu_1990[, c('fips','year','manufacturing')]
## Some counties are included in multiple files on the Census 1990 SF3 website (it chunks out counties so no one file is too huge).
## They are virtually identical, so just take mean.
manu_1990 <- manu_1990[, lapply(.SD, mean, na.rm=TRUE), .SDcols=c('manufacturing'), by=c('fips','year')]
all_manu <- rbind(manu_2015, manu_1990)
write.csv(all_manu, paste0(repo, 'covariate_clean_data/factfinder_manufacturing.csv'), row.names = FALSE)



