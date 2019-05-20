repo <- 'C:/Users/ngraetz/Documents/repos/rwjf_counties/'
source(paste0(repo, 'functions.R'))
source(paste0(repo, 'functions_shapley.R'))
cov_names <- get_clean_cov_names()
setnames(cov_names,'fe','cov')

t <- fread('C:/Users/ngraetz/Desktop/march23_tables/table1_sd.csv')
#t[, metro_region := paste0(metroname, ' ', regionname)]
setnames(t, 'metro_region_2015', 'metro_region')
t[, V1 := NULL]
t[, total_pop := NULL]
all_vars <- unique(t[, variable])
t <- dcast(t, metro_region ~ variable + year, value.var='value')
for(v in all_vars) t[, (paste0(v, '_change')) := get(paste0(v, '_2015')) - get(paste0(v, '_2000'))]
t <- t[, c('metro_region',grep('_change',names(t),value=T)), with=F]
t <- melt(t, id.vars = c('metro_region'), variable.name = 'cov', value.name = 'change')
t <- t[!(cov %in% c('nmx_change')), ]
t[, cov := gsub('_change', '', cov)]
t <- merge(t, cov_names, by='cov')
t_capped <- copy(t)
t_capped[change >= 1, change := 1]
t_capped[change <= -1, change := -1]
ordered_covs <- unique(t_capped[, c('cov_name','cov_sort')])
ordered_covs <- ordered_covs[order(cov_sort)]
ordered_covs <- ordered_covs[, cov_name]
t_capped[, cov_name := factor(cov_name, levels = ordered_covs)]

t <- fread('C:/Users/ngraetz/Desktop/march23_tables/table1.csv')
#t[, metro_region := paste0(metroname, ' ', regionname)]
setnames(t, 'metro_region_2015', 'metro_region')
t[, V1 := NULL]
t[, total_pop := NULL]
all_vars <- unique(t[, variable])
t[variable=='log_mds_pc', value := round(exp(value),1)]
t <- dcast(t, metro_region ~ variable + year, value.var='value')
for(v in all_vars) t[, (paste0(v, '_change')) := get(paste0(v, '_2015')) - get(paste0(v, '_2000'))]
t <- t[, c('metro_region',grep('_change',names(t),value=T)), with=F]
t <- melt(t, id.vars = c('metro_region'), variable.name = 'cov', value.name = 'change')
t <- t[!(cov %in% c('nmx_change')), ]
t[, cov := gsub('_change', '', cov)]
vars_to_scale <- c('chr_diabetes_monitoring','chr_mammography','perc_black','perc_hispanic','percent_transfers')
for(v in vars_to_scale) t[cov==v, change := change * 100]
t <- merge(t, cov_names, by='cov')
t[, cov_name := factor(cov_name, levels = ordered_covs)]

library(RColorBrewer)
display.brewer.pal(11, "RdBu")
t_capped[, metro_region := gsub('_',' ',metro_region)]
t[, metro_region := gsub('_',' ',metro_region)]
pdf('C:/Users/ngraetz/Documents/repos/rwjf_counties/covariate_changes_new_transfers.pdf', width=12, height=12)
ggplot() +
  geom_tile(data=t_capped,
            aes(x=cov_name,
                y=metro_region,
                fill=change),
            color='black') +
  geom_text(data=t,
            aes(x=cov_name,
                y=metro_region,
                label=round(change,2))) +
  #scale_fill_distiller(palette = 'RdBu', type='div', limits=c(-1,1), name='Change\n(standard\ndeviations)') +
  #scale_fill_manual(values=c('#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3')) + 
  scale_fill_gradient2(high = '#d6604d', low = '#4393c3', name='Change\n(standard\ndeviations)') +
  labs(x='',y='') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## DECOMPOSITION RESULTS
pdf('C:/Users/ngraetz/Documents/repos/rwjf_counties/decomp_results.pdf', width=12, height=12)
for(cause in c('allcause','drugs','cardio','lung','suicide')) {
dc <- fread(paste0("C:/Users/ngraetz/Desktop/march23_tables/",cause,"_25_64_decomp.csv"))
dc[, V1 := NULL]
dc[fe=='residual', cov_name := 'Residual']
dc[fe=='year', cov_name := 'Year']
dc <- dc[, c('metro_region','contribution_mort','cov_name')]
#dc <- merge(dc, cov_names)
#ordered_covs_dc <- c(ordered_covs, 'Year', 'Residual')
ordered_covs_dc <- cov_names[order(cov_sort)]
ordered_covs_dc <- unique(ordered_covs_dc$cov_name)
dc[, cov_name := factor(cov_name, levels = c(ordered_covs_dc,'Year'))]
dc_capped <- copy(dc)
ifelse(cause=='allcause', lim <- 30, lim <- 10)
dc_capped[contribution_mort >= lim, contribution_mort := lim]
dc_capped[contribution_mort <= -lim, contribution_mort := -lim]
decomp_gg <- ggplot() + 
  geom_tile(data=dc_capped,
            aes(x=cov_name,
                y=metro_region,
                fill=contribution_mort),
            color='black') + 
  geom_text(data=dc,
            aes(x=cov_name,
                y=metro_region,
                label=round(contribution_mort,1))) + 
  #scale_fill_distiller(palette = 'RdBu', type='div', limits=c(-lim,lim), name='Change in ASDR') +
  scale_fill_gradient2(high = '#d6604d', low = '#4393c3', name='Change in ASDR') +
  labs(x='',y='',title=paste0('Decomposition contributions for observed ', cause, ' mortality from each time-varying explanatory variable.')) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(decomp_gg)
}
dev.off()
