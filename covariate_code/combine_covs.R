library(data.table)
library(ggplot2)
library(RColorBrewer)
library(rgeos)
repo <- 'C:/Users/ngraetz/Documents/repos/rwjf_counties/'
source(paste0(repo, 'functions.R'))
source(paste0(repo, 'functions_shapley.R'))
make_maps <- FALSE

## Load all covariates and save as one big dataset long on county/year, wide on covariate.
input_dir <- paste0(repo, 'nmx_clean_data/')
cov_dir <- paste0(repo, 'covariate_clean_data/')

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

## Use all county/year/race/sex populations as template and merge everything else.
pop_white <- readRDS(paste0(input_dir, 'nhw_total_pop.RDS'))
pop_white_25_64 <- readRDS(paste0(input_dir, 'nhw_25_64_pop.RDS'))
setnames(pop_white_25_64, 'total_pop', '25_64_pop')
white_working_pop <- merge(pop_white, pop_white_25_64, by=c('fips','year','sex','race'))
white_working_pop[, perc_25_64 := `25_64_pop` / total_pop]
white_working_pop <- white_working_pop[, c('fips','year','sex','race','perc_25_64')]

pop_black <- readRDS(paste0(input_dir, 'nhb_total_pop.RDS'))
pop_black_25_64 <- readRDS(paste0(input_dir, 'nhb_25_64_pop.RDS'))
setnames(pop_black_25_64, 'total_pop', '25_64_pop')
black_working_pop <- merge(pop_black, pop_black_25_64, by=c('fips','year','sex','race'))
black_working_pop[, perc_25_64 := `25_64_pop` / total_pop]
black_working_pop <- black_working_pop[, c('fips','year','sex','race','perc_25_64')]

pop <- rbind(white_working_pop, black_working_pop)

## Merge all other covariates
for(c in c('bea_covs','bls_laus_covs','factfinder_edu','factfinder_fb','saipe_pov','ahrf_covs','census_acs_migration','factfinder_manufacturing','ihme_interpolated')) {
message(c)
cov <- fread(paste0(cov_dir,c,'.csv'))
cov[, year := as.numeric(year)]
merge_vars <- c('fips','year','sex','race')
if(!('race' %in% names(cov))) merge_vars <- merge_vars[merge_vars!='race']
if(!('sex' %in% names(cov))) merge_vars <- merge_vars[merge_vars!='sex']
## Handle Virginia FIPS in BEA
if(c=='bea_covs') {
  map <- fread(paste0(repo, 'covariate_code/bea_county_template.csv'))
  map[, fips := as.character(fips)]
  map[, bea_fips := as.character(bea_fips)]
  map[nchar(fips)==4, fips := paste0('0',fips)]
  map[is.na(bea_fips), bea_fips := fips]
  pop <- merge(pop, map, all.x=TRUE, by='fips')
  pop[, bea_fips := as.character(bea_fips)]
  pop[nchar(bea_fips)==4, bea_fips := paste0('0',bea_fips)]
  cov[, bea_fips := fips]
  cov[, fips := NULL]
  merge_vars <- c('bea_fips','year')
}
message(paste0('Merge using ', paste(merge_vars, collapse = ', ')))
pop <- merge(pop, cov, all.x=TRUE, by=merge_vars)
message(dim(pop)[1])
}

## Attach total all-race, all-sex, all-age county population for use as denominator on some covariates.
all_pops <- readRDS(paste0(cov_dir,'all_county_total_pop.RDS'))
pop <- merge(pop, all_pops, all.x=TRUE, by=c('fips','year'))
message(dim(pop)[1])

## Construct any variables where we need to divide by population. Also rescale things that are percents.
pop[, mds_pc := (new_total_mds / total_county_pop) * 1000]
pop[, perc_labor := (labor_force / total_county_pop) * 100]
pop[, fb := fb * 100]
pop[, manufacturing := manufacturing * 100]
pop[, college := college * 100]
pop[, less_12 := less_12 * 100]
pop[, percent_transfers := percent_transfers * 100]
pop[, perc_25_64 := perc_25_64 * 100]
pop[, net_mig_per1000 := (net_mig / total_county_pop) * 1000]
pop[, in_mig_per1000 := (in_mig / total_county_pop) * 1000]
pop[, out_mig_per1000 := (out_mig / total_county_pop) * 1000]
pop[, log_hh_income := log(median_hh_income)]

## Some labor force is over 100% in small counties. Cap for now.
pop[!is.na(perc_labor) & perc_labor >= 100, perc_labor := 100]

## Make maps of everything.
# all_covs <- c('percent_transfers','percent_wage_salary_employment','income_per_capita',
#               'percent_unemployment','less_12','college','fb','poverty_all','mds_pc','perc_labor',
#               'net_mig_per1000','in_mig_per1000','out_mig_per1000')
all_covs <- c("as_diabetes_prev","pa_prev","obesity_prev","as_heavy_drinking_prev","current_smoker_prev")
if(make_maps==TRUE) {
pdf(paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/prep_plots/all_combined_maps.pdf'), height=6, width=9)
for(v in all_covs) {
  message(paste0('Mapping ', v, '...'))
  lower <- quantile(pop[year %in% c(2000,2010,2015), get(v)], probs = 0.05, na.rm=TRUE)
  upper <- quantile(pop[year %in% c(2000,2010,2015), get(v)], probs = 0.95, na.rm=TRUE)
  for(y in c(2000,2010,2015)) {
    m <- make_county_map(map_dt = pop[year==y & race==0 & sex==1, ],
                         map_sp = counties,
                         map_var = v,
                         legend_title = v,
                         high_is_good = FALSE,
                         map_title = paste0(y, ' ', v),
                         map_limits = c(lower,upper)) + guides(fill=guide_legend(title='Value'))
    print(m)
  }
}
dev.off()
}

## Create trend plots of everything by metro (color) and region (facet)
metro_codes <- fread(paste0(repo, 'covariate_clean_data/FIPSmetroregion.csv'))
metro_codes[, fips := as.character(fips)]
metro_codes[nchar(fips)==4, fips := paste0('0',fips)]
trends <- merge(pop, metro_codes[, c('fips','metroname','regionname')], by=c('fips'), all.x=TRUE)
trends[metroname %in% c('Nonmetro, adjacent', 'Nonmetro, nonadjacent'), metroname := 'Nonmetro']
race_pop <- readRDS(paste0(input_dir, 'nhw_total_pop.RDS'))
trends <- merge(trends, race_pop, by=c('fips','year','sex','race'))
trends <- trends[race==0 & sex==1 & year %in% c(1990,2000,2010,2015), lapply(.SD, weighted.mean, w=total_pop, na.rm=TRUE), .SDcols=all_covs, by=c('metroname','regionname','year')]

  cov_names <- get_clean_cov_names()
  pdf(paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/prep_plots/all_combined_trends.pdf'), width = 12, height = 8)
  for(c in all_covs) {
    trends[is.nan(get(c)), (c) := NA]
    gg <- ggplot() + 
      geom_line(data=trends,
                aes(x=year,
                    y=get(c),
                    color=metroname),
                size=2) +
      theme_minimal() + 
      scale_color_manual(values = brewer.pal(4,'Dark2'), name='Metro') + 
      labs(x='Year',y=cov_names[fe==c, cov_name], title=cov_names[fe==c, cov_name]) + 
      facet_wrap(~regionname)
    print(gg)
  }
  dev.off()

## Save.
saveRDS(pop, paste0(repo, 'covariate_clean_data/combined_covs.RDS'))

## Make covariance matrix heat map
cov_names <- c('college','poverty_all','log_hh_income','percent_transfers',
               'percent_unemployment','mds_pc',
               'as_diabetes_prev','pa_prev','obesity_prev','as_heavy_drinking_prev','current_smoker_prev',
               'fb','perc_25_64')
cov_matrix <- cor(pop[, cov_names, with=F], use='complete.obs')
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cov_matrix)
upper_tri <- melt(upper_tri, na.rm=TRUE)
clean_cov_names <- get_clean_cov_names()
upper_tri$id <- 1:length(upper_tri$value)
upper_tri <- as.data.table(upper_tri)
setnames(upper_tri, 'Var1', 'fe')
upper_tri <- merge(upper_tri, clean_cov_names[, c('fe','cov_name')], all.x=TRUE, by='fe')
setnames(upper_tri, c('Var2','cov_name'), c('fe2','cov_name1'))
setnames(clean_cov_names, 'fe', 'fe2')
upper_tri <- merge(upper_tri, clean_cov_names[, c('fe2','cov_name')], all.x=TRUE, by='fe2')
setnames(upper_tri, 'cov_name', 'cov_name2')
upper_tri <- upper_tri[order(id)]
upper_tri[, cov_name2 := factor(cov_name2, levels = unique(upper_tri[order(id), cov_name2]))]
upper_tri[, cov_name1 := factor(cov_name1, levels = unique(upper_tri[order(id), cov_name1]))]
# upper_tri$fe <- upper_tri$Var1
# upper_tri <- merge(upper_tri, clean_cov_names[, c('fe','cov_name')], all.x=TRUE, by='fe')
# upper_tri$fe <- upper_tri$Var2
# names(upper_tri)[names(upper_tri)=='cov_name'] <- 'cov_name1'
# upper_tri <- merge(upper_tri, clean_cov_names[, c('fe','cov_name')], all.x=TRUE, by='fe')
# names(upper_tri)[names(upper_tri)=='cov_name'] <- 'cov_name2'
# upper_tri <- upper_tri[order(id),]
pdf(paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/prep_plots/all_combined_cor_matrix.pdf'), width = 12, height = 8)
ggplot(data = upper_tri, aes(cov_name2, cov_name1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value,2))) + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  labs(y='',x='') + 
  coord_fixed()
dev.off()
