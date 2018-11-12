make_permutations <- function(fes, start_year, end_year) {
  fes <- c(fes, 'year','residual')
  permutations <- list()
  for(fe in 1:length(fes)) permutations[[fe]] <- c(start_year,end_year)
  permutations <- as.data.table(expand.grid(permutations))
  setnames(permutations, names(permutations), fes)
  return(permutations)
}

## Predict on all permutations.
calculate_contribution <- function(fe, fes, start_year, end_year, dt, coefs, all_permutations) {
  ## Grab all permutations where this fixed effect is changing (2010) and calculate difference vs. difference if it had not changed (1990).
  ## The difference of these two differences is the "contribution" of change in that fixed effect WITHIN this permutation (though this
  ## difference seems to be identical across permutations).
  message(paste0('Calculating contribution from ', fe, '...'))
  fe_permutations <- all_permutations[get(fe)==end_year, ]
  other_fes <- c(fes, 'year','residual')
  other_fes <- other_fes[other_fes!=fe]
  inv_logit <- function(x) {exp(x)/(1+exp(x))}
  calculate_permutation <- function(p, perm_dt, max_perm) {

    ## Message progress
    if(nchar(as.character(p/100))==1) message(paste0(as.character(p),' / ', dim(fe_permutations)[1]))

    ## Select permutation.
    this_dt <- copy(perm_dt)
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
    if(p==1) message('Age groups: ', paste(unique(this_dt[, agegrp]), sep=' '))
    for(a in unique(this_dt[, agegrp])) {
      if(a==min(unique(this_dt[, agegrp]))) this_dt[agegrp==a, age_int := 0]
      if(a!=min(unique(this_dt[, agegrp]))) this_dt[agegrp==a, age_int := coefs[name==paste0('as.factor(agegrp)',a), coef]]
    }
    this_dt[, (paste0('p_with_change_',p)) := get((paste0('p_with_change_',p))) + coefs[name=='(Intercept)', coef] + age_int + metro_region_int + get(paste0('spatial_effect_',start_year))]
    this_dt[, (paste0('p_without_change_',p)) := get((paste0('p_without_change_',p))) + coefs[name=='(Intercept)', coef] + age_int + metro_region_int + get(paste0('spatial_effect_',start_year))]

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
    this_dt <- this_dt[, c('fips', 'agegrp', 'p_with_change', 'p_without_change', 'diff')]
    this_dt[, p := p]
    return(this_dt)

  }

  message('Calculating all permutations...')
  all_diffs <- as.data.table(rbind.fill(lapply(1:dim(fe_permutations)[1], calculate_permutation, perm_dt=dt, max_perm=dim(fe_permutations)[1])))
  ## As this is a Shapley decomposition, here is where we "average over" potential path dependencies (i.e. all the different permutations).
  ## A more complex generalized decomposition could be used, such as g-computation (actually estimate all the path dependencies, decompose
  ## direct/indirect change attributable via bootstrap and stochastic simulation through periods.
  all_diffs <- all_diffs[, list(contribution_mort=mean(diff)), by=c('fips','agegrp')]
  all_diffs[, fe := fe]
  return(all_diffs)

}


get_clean_cov_names <- function() {
  metro_names <- c('Lg central metro', 'Lg fringe metro','Md/Sm metro','Nonmetro')
  region_names <- c('Pacific', 'Appalachia', "East South Central","Mountain","West South Central","New England","South Atlantic","East North Central","West North Central","Middle Atlantic")
  metro_region_names <- apply(expand.grid(metro_names, region_names), 1, paste, collapse="_")
  new_covs <- c("percent_wage_salary_employment","income_per_capita","less_12","college","poverty_all",
                "percent_transfers","percent_unemployment","perc_25_64","fb","perc_labor","mds_pc",
                'net_mig_per1000','in_mig_per1000','out_mig_per1000','obesity','net_in_mig','manufacturing',
                "as_diabetes_prev","pa_prev","obesity_prev","as_heavy_drinking_prev","current_smoker_prev",'log_hh_income',
                'log_mds_pc','chr_mammography','chr_diabetes_monitoring')
  new_cov_names <- c("Perc wage vs salary","Income/pc","Less HS",'College','Poverty','Transfers','Unemployment',
                     'Percent 25-64','Foreign-born','Labor force','MDs/pc','Net-mig/1000','In-mig/1000','Out-mig/1000',
                     'Obesity','Net In-mig','Manufacturing',
                     "Diabetes","Physical activity","Obesity","Heavy drinking","Smoking","HH income",
                     'MDs/pc','Mammography','Diabetes monitoring')
  cov_names <- data.table(fe = c(paste0('as.factor(year)',c(1990,2000,2010)),
                                 paste0('as.factor(metro_region)',metro_region_names),
                                 'air_EQI_22July2013','water_EQI_22July2013','land_EQI_22July2013',
                                 new_covs, 'tbhk_pv','tbhk_is','tbhk_ce','fbp_kt','twhk_pv','twhk_is','twhk_ce','fwp_kt','metro',"Global Moran's I","DIC",'RMSE','year','residual'),
                          cov_name = c('1990','2000','2010', metro_region_names, 'Air quality','Water quality','Land quality',
                                       new_cov_names, 'Poverty','Isolation','College education',
                                       'Foreign-born','Poverty','Isolation','College education','Foreign-born', 'Metro level',"Global Moran's I","DIC",'RMSE','Secular trend','Residual'),
                          cov_sort = c(1:86))
}


logit <- function(x) {
  return(log(x / (1-x)))
}
inv_logit <- function(x) {
  return(exp(x)/(1+exp(x)))
}


shapley_ages <- function(ages,data,inla_f,coef_file,shapley,shap_covs) {

  mort <- copy(data[agegrp %in% ages[1]:ages[2], ])
  message(unique(mort[, agegrp]))
  inla_model = inla(as.formula(inla_f),
                    family = "binomial",
                    data = mort,
                    Ntrials = mort[, total_pop],
                    verbose = FALSE,
                    control.compute=list(config = TRUE, dic = TRUE),
                    control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
                    control.fixed=list(prec.intercept = 0,
                                       prec = 1),
                    num.threads = 4)

  ## Save coefs table
  coefs <- make_beta_table(inla_model, paste0('age_',  ages[1], '_',  ages[2], '_', sex_option))
  saveRDS(coefs, file = paste0(repo,'/coefs/age_',  ages[1], '_',  ages[2], '_', coef_file,'.RDS'))
  message(paste0('Saving coefs: ',repo,'/coefs/age_',  ages[1], '_',  ages[2], '_',coef_file,'.RDS'))

  if(shapley==FALSE) {
    return(NULL)
  }
  if(shapley==TRUE) {
  ## Make full prediction, full residual, and load posterior mean for all components.
  mort[, inla_pred := inla_model$summary.fitted.values$mean]
  mort[, inla_residual := logit(nmx) - logit(inla_pred)]
  mort[nmx==0, inla_residual := logit(nmx+0.000001) - logit(inla_pred)]
  model_coefs <- make_beta_table(inla_model, paste0(race,' ',sex_option,' ',domain))
  #saveRDS(coefs, file = paste0(out_dir,'coefs_', output_name, '.RDS'))

  ## Run Shapley decomposition on change over time in ASDR.
  ## Create permutations (2010-1990, 6 changes, total permutations = 2^6 = 64, 32 pairs)
  ## i.e. one pair for poverty is delta_m|PV=2013,IS=1990,CE=1990,FB=1990,time=1990,residual=1990 -
  ##                              delta_m|PV=1990,IS=1990,CE=1990,FB=1990,time=1990,residual=1990
  permutations <- make_permutations(fes = shap_covs,
                                    start_year = start_year,
                                    end_year = end_year)
  message(dim(permutations))

  ## Prep and reshape input data from model (all fixed effects of interest + geographic random effects +
  ## time + spatial random effects + intercept + residual, wide by year)
  #d <- copy(mort)
  shap_d <- merge(mort, inla_model$summary.random$ID[c('ID','mean')], by='ID')
  setnames(shap_d, 'mean', 'spatial_effect')
  shap_d <- shap_d[year %in% c(start_year,end_year), ]
  shap_d <- dcast(shap_d, fips + agegrp ~ year, value.var = c(shap_covs, 'inla_residual', 'inla_pred', 'total_pop', 'metro_region', 'nmx', 'spatial_effect'))
  shap_d[, (paste0('year_', start_year)) := start_year]
  shap_d[, (paste0('year_', end_year)) := end_year]

  ## Calculate contribution attributable to each time-varying component. By definition this adds up to total observed change in the outcome
  ## because of inclusion of the residual.
  ## Change decompositions occur at the county-age-level.
  message(paste(shap_covs,collapse=' '))
  all_contributions <- rbindlist(lapply(c(shap_covs, 'year','residual'), calculate_contribution,
                                        fes=shap_covs,
                                        start_year=start_year,
                                        end_year=end_year,
                                        dt=shap_d,
                                        coefs=model_coefs,
                                        all_permutations=permutations))
  return(all_contributions)
  }
}


model_permute <- function(ages,data,inla_f,coef_file,shapley,shap_covs,perm) {

  mort <- copy(data[agegrp %in% ages[1]:ages[2], ])
  inla_model = inla(as.formula(inla_f),
                    family = "binomial",
                    data = mort,
                    Ntrials = mort[, total_pop],
                    verbose = FALSE,
                    control.compute=list(config = TRUE, dic = TRUE),
                    control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
                    control.fixed=list(prec.intercept = 0,
                                       prec = 1),
                    num.threads = 4)

  ## Save coefs table
  coefs <- make_beta_table(inla_model, paste0('age_',  ages[1], '_',  ages[2], '_', sex_option))
  coefs[, coef := exp(coef)]
  mort[, inla_pred := inla_model$summary.fitted.values$mean]
  mort[, inla_residual := logit(nmx) - logit(inla_pred)]
  mort[nmx==0, inla_residual := logit(nmx+0.000001) - logit(inla_pred)]
  all_fit <- data.table(model = rep(paste0('age_',  ages[1], '_',  ages[2], '_', sex_option), 3),
                        name = c('DIC','RMSE','R2'),
                        coef = rep(.00001,3))
  all_fit[name == "DIC", coef := ifelse(is.nan(inla_model$dic$dic) | is.infinite(inla_model$dic$dic), inla_model$dic$deviance.mean, inla_model$dic$dic)]
  all_fit[name == 'RMSE', coef := mort[, sqrt(weighted.mean(get(paste0('inla_residual'))^2, w = total_pop))]]
  all_fit[name == 'R2', coef := 1 - (sum((mort$nmx-inla_model$summary.fitted.values$mean)^2)/sum((mort$nmx-mean(mort$nmx))^2))]
  coefs <- rbind(coefs, all_fit, fill=T)
  coefs <- coefs[name %in% c(covs,'DIC','RMSE','R2'), c('model','name','coef')]
  coefs[, p := perm]
  return(coefs)
}



