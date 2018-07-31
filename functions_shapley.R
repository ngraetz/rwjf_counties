make_permutations <- function(fes, start_year, end_year) {
  fes <- c(fes, 'year','residual')
  permutations <- list()
  for(fe in 1:length(fes)) permutations[[fe]] <- c(start_year,end_year)
  permutations <- as.data.table(expand.grid(permutations))
  setnames(permutations, names(permutations), fes)
  return(permutations)
}

## Predict on all permutations.
calculate_contribution <- function(fe, fes, start_year, end_year) {
  ## Grab all permutations where this fixed effect is changing (2010) and calculate difference vs. difference if it had not changed (1990).
  ## The difference of these two differences is the "contribution" of change in that fixed effect WITHIN this permutation (though this
  ## difference seems to be identical across permutations).
  message(paste0('Calculating contribution from ', fe, '...'))
  fe_permutations <- permutations[get(fe)==end_year, ]
  other_fes <- c(fes, 'year','residual')
  other_fes <- other_fes[other_fes!=fe]
  inv_logit <- function(x) {exp(x)/(1+exp(x))}
  calculate_permutation <- function(p, dt) {
    
    ## Select permutation.
    message(p)
    this_dt <- copy(dt)
    p_option <- fe_permutations[p,]
    
    ## Assign values to be added from all effects besides the target effect.
    for(other_fe in other_fes) {
      if(!(other_fe %in% c('residual'))) {
        ## Assign relevant fixed effect values.
        if(p_option[, get(other_fe)]==end_year) this_dt[, (paste0('p_',other_fe)) := get(paste0(other_fe,'_', end_year)) * coefs[name==other_fe, coef]]
        if(p_option[, get(other_fe)]==start_year) this_dt[, (paste0('p_',other_fe)) := get(paste0(other_fe,'_', start_year)) * coefs[name==other_fe, coef]]
      }
      if(other_fe=='residual') {
        ## Assign relevant value for residual.
        if(p_option[, residual]==end_year) this_dt[, p_residual := get(paste0('inla_residual_',end_year))] ## MAKE SURE RESIDUAL IS IN LOGIT SPACE - it comes in normal space from INLA object.
        if(p_option[, residual]==start_year) this_dt[, p_residual := get(paste0('inla_residual_',start_year))]
      }
    }
    
    ## Generate full prediction for this permutation based on whether FE value stayed the same or changed over the period.
    ## Assign target FE value based on change (2010 value) or no change (1990 value), and then add in all other effects.
    if(!(fe %in% c('residual'))) {
      this_dt[, (paste0('p_with_change_',p)) := (get(paste0(fe,'_',end_year)) * coefs[name==fe, coef])]
      this_dt[, (paste0('p_without_change_',p)) := (get(paste0(fe,'_',start_year)) * coefs[name==fe, coef])]
    }
    if(fe=='residual') {
      this_dt[, (paste0('p_with_change_',p)) := get(paste0('inla_residual_',end_year))]
      this_dt[, (paste0('p_without_change_',p)) := get(paste0('inla_residual_',start_year))]
    }
    
    ## Add intercept and random effects.
    for(r in unique(this_dt[, get(paste0('metro_region_',end_year))])) {
      if(r=='Lg central metro_Pacific') this_dt[get(paste0('metro_region_',end_year))==r, metro_region_int := 0]
      if(r!='Lg central metro_Pacific') this_dt[get(paste0('metro_region_',end_year))==r, metro_region_int := coefs[grep(r,name), coef]]
    }
    this_dt[, (paste0('p_with_change_',p)) := get((paste0('p_with_change_',p))) + coefs[name=='(Intercept)', coef] + metro_region_int + get(paste0('spatial_effect_',start_year))]
    this_dt[, (paste0('p_without_change_',p)) := get((paste0('p_without_change_',p))) + coefs[name=='(Intercept)', coef] + metro_region_int + get(paste0('spatial_effect_',start_year))]
    
    ## This change (2010 prediction based just on FE of interest - 1990 prediction) is the same across all permutations... 
    ## But then we are just adding a constant to both sides derived from this specific permutation of all other effects...?
    #message(this_dt[, get((paste0('p_with_change_',p)))][1])
    for(other_fe in other_fes) {
      this_dt[, (paste0('p_with_change_',p)) := get((paste0('p_with_change_',p))) + get(paste0('p_',other_fe))]
      this_dt[, (paste0('p_without_change_',p)) := get((paste0('p_without_change_',p))) + get(paste0('p_',other_fe))]
    }
    
    ## Generate difference in full prediction for this permutation attributable to change in this FE value over the period.
    ## The difference attributable to this effect in this permutation needs to be calculated in normal space. This is how we handle non-linearities.
    ## If we decompose life expectancy, here is where you would want to convert both scenarios in this permutation to normal space, 
    ## calculate life expectancies, and return the difference.
    this_dt[, diff := inv_logit(get(paste0('p_with_change_',p))) - inv_logit(get(paste0('p_without_change_',p)))]
    this_dt[, p_with_change := get(paste0('p_with_change_',p))]
    this_dt[, p_without_change := get(paste0('p_without_change_',p))]
    this_dt <- this_dt[, c('fips', 'p_with_change', 'p_without_change', 'diff')]
    this_dt[, p := p]
    return(this_dt)
    
  }
  
  message('Calculating all permutations...')
  all_diffs <- as.data.table(rbind.fill(lapply(1:dim(fe_permutations)[1], calculate_permutation, dt=d)))
  ## As this is a Shapley decomposition, here is where we "average over" potential path dependencies (i.e. all the different permutations).
  ## A more complex generalized decomposition could be used, such as g-computation (actually estimate all the path dependencies, decompose
  ## direct/indirect change attributable via bootstrap and stochastic simulation through periods.
  all_diffs <- all_diffs[, list(contribution_mort=mean(diff)), by=c('fips')]
  all_diffs[, fe := fe]
  return(all_diffs)
  
}


get_clean_cov_names <- function() {
  metro_names <- c('Lg central metro', 'Lg fringe metro','Md/Sm metro','Nonmetro')
  region_names <- c('Pacific', 'Appalachia', "East South Central","Mountain","West South Central","New England","South Atlantic","East North Central","West North Central","Middle Atlantic")
  metro_region_names <- apply(expand.grid(metro_names, region_names), 1, paste, collapse="_")
  new_covs <- c("college","poverty_all","percent_transfers","percent_unemployment","perc_25_64","fb")
  new_cov_names <- c('College','Poverty','Transfers','Unemployment','Percent 25-64','Foreign-born')
  cov_names <- data.table(fe = c(paste0('as.factor(year)',c(1990,2000,2010)),
                                 paste0('as.factor(metro_region)',metro_region_names),
                                 'air_EQI_22July2013','water_EQI_22July2013','land_EQI_22July2013',
                                 new_covs, 'tbhk_pv','tbhk_is','tbhk_ce','fbp_kt','twhk_pv','twhk_is','twhk_ce','fwp_kt','metro',"Global Moran's I","DIC",'RMSE','year','residual'),
                          cov_name = c('1990','2000','2010', metro_region_names, 'Air quality','Water quality','Land quality',
                                       new_cov_names, 'Poverty','Isolation','College education',
                                       'Foreign-born','Poverty','Isolation','College education','Foreign-born', 'Metro level',"Global Moran's I","DIC",'RMSE','Secular trend','Residual'),
                          cov_sort = c(1:66))
}
