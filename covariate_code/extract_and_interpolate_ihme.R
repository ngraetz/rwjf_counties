library(ggplot2)
library(data.table)
library(grid)
library(gridExtra)
repo <- 'C:/Users/ngraetz/Documents/repos/rwjf_counties/'
cov_dir <- paste0(repo, 'covariate_clean_data/')

smoking <- "C:/Users/ngraetz/Documents/repos/rwjf_counties/covariate_raw/ihme/IHME_US_COUNTY_TOTAL_AND_DAILY_SMOKING_PREVALENCE_1996_2012/IHME_US_COUNTY_TOTAL_AND_DAILY_SMOKING_PREVALENCE_1996_2012.csv"
alcohol <- "C:/Users/ngraetz/Documents/repos/rwjf_counties/covariate_raw/ihme/IHME_USA_COUNTY_ALCOHOL_USE_PREVALENCE_2002_2012_NATIONAL/heavy_prop_IHME_USA_COUNTY_ALCOHOL_USE_PREVALENCE_2002_2012_NATIONAL_Y2015M04D23.csv"
diabetes <- "C:/Users/ngraetz/Documents/repos/rwjf_counties/covariate_raw/ihme/IHME_USA_COUNTY_DIABETES_PREVALENCE_1999_2012/age_standardized_prev_IHME_USA_COUNTY_DIABETES_PREVALENCE_1999_2012_NATIONAL_Y2016M08D23.csv"
pa_obesity <- "C:/Users/ngraetz/Documents/repos/rwjf_counties/covariate_raw/ihme/IHME_USA_OBESITY_PHYSICAL_ACTIVITY_2001_2011 (1).csv"

## Read and format to the same shape.
smoking <- fread(smoking)
setnames(smoking, c('total_mean','sex'), c('current_smoker_prev','sex_name'))
smoking[sex_name=='Males', sex := 1]
smoking[sex_name=='Females', sex := 2]
smoking <- smoking[!is.na(sex), c('state','county','current_smoker_prev','year','sex')]

alcohol <- fread(alcohol)
alcohol[, county := Location]
alcohol <- melt(alcohol, id.vars = c('county','State'), measure.vars = c(paste0(2005:2012, ' Females'), paste0(2005:2012, ' Males')))
alcohol[, c("year","sex_name") := tstrsplit(variable, " ", fixed=TRUE)]
alcohol[, year := as.numeric(year)]
setnames(alcohol, c('State','value'), c('state','as_heavy_drinking_prev'))
alcohol[sex_name=='Males', sex := 1]
alcohol[sex_name=='Females', sex := 2]
alcohol <- alcohol[!is.na(sex), c('state','county','as_heavy_drinking_prev','year','sex')]

diabetes <- fread(diabetes, skip = 1)
diabetes[, fips := as.character(FIPS)]
diabetes[nchar(fips)==4, fips := paste0('0',fips)]
diabetes[, county := Location]
diabetes <- melt(diabetes, id.vars = c('county','fips'), names(diabetes)[grep('Prevalence',names(diabetes))])
diabetes[, variable := as.character(variable)]
diabetes[, c("measure","year","sex_name") := tstrsplit(variable, ", ", fixed=TRUE)]
diabetes[, year := as.numeric(year)]
setnames(diabetes, 'value', 'as_diabetes_prev')
diabetes[sex_name=='Males', sex := 1]
diabetes[sex_name=='Females', sex := 2]
diabetes <- diabetes[!is.na(sex), c('county','fips','as_diabetes_prev','year','sex')]

pa_obesity <- fread(pa_obesity)
pa_obesity[, fips := as.character(fips)]
pa_obesity[nchar(fips)==4, fips := paste0('0',fips)]
pa_obesity <- melt(pa_obesity, id.vars = c('fips','Sex','Outcome','State'), measure.vars = c("Prevalence 2001 (%)","Prevalence 2002 (%)","Prevalence 2003 (%)","Prevalence 2004 (%)",
                                                                  "Prevalence 2005 (%)","Prevalence 2006 (%)","Prevalence 2007 (%)","Prevalence 2008 (%)",
                                                                  "Prevalence 2009 (%)","Prevalence 2010 (%)","Prevalence 2011 (%)"))
pa_obesity[, c("measure","year","trash") := tstrsplit(variable, " ", fixed=TRUE)]
pa_obesity[, year := as.numeric(year)]
pa_obesity <- dcast(pa_obesity, fips + Sex + year + State ~ Outcome, value.var = 'value')
pa_obesity[Sex=='Male', sex := 1]
pa_obesity[Sex=='Female', sex := 2]
setnames(pa_obesity, c('State','Any PA','Obesity','Sufficient PA'), c('state','pa_prev','obesity_prev','sufficient_pa_prev'))
pa_obesity <- pa_obesity[!is.na(sex), c('state','fips','pa_prev','obesity_prev','sufficient_pa_prev','year','sex')]

## Merge all together using fips/county from diabetes as anchor.
template <- as.data.table(expand.grid(2000:2015, 1:2, unique(diabetes[!is.na(fips), fips])))
setnames(template, c('year','sex','fips'))
template[, fips := as.character(fips)]
template[nchar(fips)==4, fips := paste0('0',fips)]
all <- merge(template, diabetes, by=c('fips','year','sex'), all.x=TRUE)
all <- merge(all, pa_obesity, by=c('fips','year','sex'), all.x=TRUE)
all <- merge(all, alcohol, by=c('state','county','year','sex'), all.x=TRUE)
all <- merge(all, smoking, by=c('state','county','year','sex'), all.x=TRUE)
all[, state := NULL]
all[, county := NULL]

## Merge CHR variables to compare
chr <- fread("C:/Users/ngraetz/Documents/repos/rwjf_counties/covariate_clean_data/chr_covs.csv")
all <- merge(all, chr, all.x=TRUE, by=c('fips','year'))

## Interpolate forward and backward assuming growth rate is constant.
cols = c('as_diabetes_prev','pa_prev','obesity_prev','as_heavy_drinking_prev','current_smoker_prev',"chr_obesity_prev","chr_mammography","chr_diabetes_monitoring")
anscols = paste("lag", cols, sep="_")
all <- all[order(year)]
all[order(year), (anscols) := data.table::shift(.SD, 1, type="lag"), .SDcols=cols, by=c('fips','sex')]
for(v in cols) all[, paste0('r_',(v)) := log(get(v)/get(paste0('lag_',v)))]
magr <- all[, lapply(.SD, mean, na.rm=TRUE), by=c('fips','sex'), .SDcols=paste0('r_',cols)]
setnames(magr, paste0('r_',cols), paste0('magr_',cols))
all <- merge(all, magr, by=c('fips','sex'))

proj_cov <- function(x, magr) {
  if(length(x[!(is.na(x))])==0) {
    return(as.numeric(rep(NA,16)))
  }
  if(length(x[!(is.na(x))])!=0) {
    temp <- data.table(year=2000:2015, cov=x, magr=magr)
    ## Get year range of this variable.
    min_y <- min(temp[!is.na(cov), year])
    max_y <- max(temp[!is.na(cov), year])
    ## Interpolate backwards.
    if(min_y != 2000) { 
      for(y in rev(2000:(min_y-1))) {
        temp[, trend := data.table::shift(cov, type='lead')]
        temp[year==y, cov := trend / exp(magr)]
      }
    }
    ## Interpolate forwards.
    if(max_y != 2015) {
      for(y in (max_y+1):2015) {
        temp[, trend := data.table::shift(cov, type='lag')]
        temp[year==y, cov := trend * exp(magr)]
      }
    }
    ## Return just the interpolated covariate trend to assign to the new column.
    return(temp[, cov])
  }
}
all_melt <- melt(all, id.vars = c('fips','year','sex',c(paste0('magr_',cols))), measure.vars = cols, value.name = 'value')
all_melt <- melt(all_melt, id.vars = c('fips','year','sex','variable','value'), measure.vars = c(paste0('magr_',cols)), value.name = 'magr')
all_melt[, variable.1 := gsub('magr_','',variable.1)]
all_melt <- all_melt[variable==variable.1, ]
all_melt[, variable.1 := NULL]
all_melt[order(year), int_value := proj_cov(value, magr), by=c('fips','sex','variable')]
all_melt[!is.na(int_value) & is.na(value), interpolated := 'int']
all_melt[is.na(interpolated), interpolated := 'raw']
all_melt[variable=='chr_obesity_prev', value := value * 100]
all_melt[variable=='chr_obesity_prev', int_value := int_value * 100]

## Save
final_int <- dcast(all_melt, fips + year + sex ~ variable, value.var = c('int_value'))
write.csv(final_int, paste0(cov_dir,'/ihme_interpolated.csv'), row.names = FALSE)

## Make plots of interpolations.
final_int <- dcast(all_melt, fips + year + sex ~ variable, value.var = c('value','int_value'))
metro_codes <- fread(paste0(repo, 'covariate_clean_data/FIPSmetroregion.csv'))
metro_codes[, fips := as.character(fips)]
metro_codes[nchar(fips)==4, fips := paste0('0',fips)]
metro_codes[metroname %in% c('Nonmetro, adjacent', 'Nonmetro, nonadjacent'), metroname := 'Nonmetro']
final_int <- merge(final_int, metro_codes, by='fips')
final_int[, metro_region := paste0(metroname, ' - ', regionname)]
all_pops <- readRDS(paste0(cov_dir,'all_county_total_pop.RDS'))
final_int <- merge(final_int, all_pops, by=c('year','fips'))

all_covs <- c('as_diabetes_prev','pa_prev','obesity_prev','as_heavy_drinking_prev','current_smoker_prev',"chr_obesity_prev","chr_mammography","chr_diabetes_monitoring")
all_covs <- c('obesity_prev','chr_obesity_prev','as_diabetes_prev')
png(paste0('C:/Users/ngraetz/Dropbox/Penn/papers/rwjf/covariates/prep_plots/all_interpolated_trends_regions.png'))
final_int[, metro_region := paste0(metroname, ' - ', regionname)]
final_int[, metro_region := regionname]
for(r in unique(final_int[, metro_region])) {
sub <- final_int[metro_region==r, ]
## Collapse over sex to fips or metro-region
sub_agg <- sub[, lapply(.SD, weighted.mean, w=total_county_pop, na.rm=TRUE), .SDcols=c(paste0('int_value_',all_covs),paste0('value_',all_covs)), by=c('metro_region','year')]
sub <- sub[, lapply(.SD, weighted.mean, w=total_county_pop, na.rm=TRUE), .SDcols=c(paste0('int_value_',all_covs),paste0('value_',all_covs)), by=c('fips','year')]
get_cov_gg <- function(c) {
if(grepl('obesity',c)) y_lims <- c(10,60)
if(grepl('diabetes',c)) y_lims <- c(0,25)
if(c=='value_obesity_prev') y_lab <- 'IHME obesity prevalence'
if(c=='value_as_diabetes_prev') y_lab <- 'IHME diabetes prevalence'
if(c=='value_chr_obesity_prev') y_lab <- 'CHR obesity prevalence (BRFSS)'
gg <- ggplot() + 
  geom_line(data=sub,
            aes(x=year,
                y=get(paste0('int_',c)),
                group=fips),
            linetype = 'dashed',
            alpha = 0.1) +
  geom_line(data=sub,
            aes(x=year,
                y=get(c),
                group=fips),
            alpha = 0.4,
            size = 1) +
  labs(x='',y=gsub('value_','',c)) +
  geom_line(data=sub_agg,
            aes(x=year,
                y=get(c)),
            size=2,
            color='red') +
  geom_line(data=sub_agg,
            aes(x=year,
                y=get(paste0('int_',c))),
            linetype = 'dashed',
            size = 2,
            color = 'red') +
  lims(y=y_lims) + 
  labs(y=y_lab) + 
  theme_minimal()
return(gg) 
}
all_cov_ggs <- lapply(paste0('value_',all_covs), get_cov_gg)
grid.arrange(
  top = textGrob(r,gp=gpar(fontsize=20,font=3)),
  grobs = all_cov_ggs,
  widths = c(1,1,1),
  layout_matrix = rbind(c(1,2,3))
)
}
dev.off()


ggplot() + 
  geom_point(data=final_int[year==2012,],
                      aes(x=int_value_obesity_prev,
                          y=int_value_chr_obesity_prev)) + 
  geom_abline(slope=1,intercept=0) + 
  theme_minimal()




