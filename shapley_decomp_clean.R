library(INLA)
library(data.table)
library(ggplot2)

run_shapley_decomp <- function(race_option, sex_option) {
  
message(paste0('Decomposing ', race_option, ' ', sex_option))
  
  start_year <- 1990
  end_year <- 2010
  if(race_option==0) fes <- c('twhk_pv','twhk_is','twhk_ce','fwp_kt')
  if(race_option==1) fes <- c('tbhk_pv','tbhk_is','tbhk_ce','fbp_kt')
    
## Load cleaned dataset for this set of arguments.
  model_data_black <- readRDS(paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/model_data_race', race_option, '_sex', sex_option, '.RDS'))
  message('Counties in data: ', length(unique(model_data_black[, fips])))

## Drop if missing any covariates and log.
  drop_counties <- c()
  for(fe in fes) drop_counties <- c(drop_counties, model_data_black[is.na(get(fe)), fips])
  drop_counties <- c(drop_counties, model_data_black[is.na(metroname) | is.na(regionname), fips])
  model_data_black <- model_data_black[!(fips %in% unique(drop_counties)), ]
  message('Counties after dropping missing values in any year: ', length(unique(model_data_black[, fips])))
  model_data_black[metroname %in% c('Nonmetro, adjacent','Nonmetro, nonadjacent'), metroname := 'Nonmetro']
  model_data_black[, regionname := factor(regionname, levels=c('Pacific',unique(model_data_black[regionname!='Pacific', regionname])))]
  model_data_black[, metroname := factor(metroname, levels=c('Lg central metro',unique(model_data_black[metroname!='Lg central metro', metroname])))]

message('Creating spatial weight matrix...')

## Create a neighborhood file using Queens convention.
  nb.FOQ <- poly2nb(counties, queen=TRUE)
## Create an INLA weight matrix.
  lw.FOQ <- nb2INLA("FOQ_INLA",nb.FOQ)
  an <- data.frame(1:length(counties),counties$fips)
  o <- match(model_data_black$fips,an$counties.fips)
  ID <- an$X1.length.counties.[o]
  model_data_black[, ID := ID]
  
## Make formulas
  if(race_option == 0) {
    inla_formula4 <- deaths ~ 
      twhk_pv + twhk_is + twhk_ce + fwp_kt + as.factor(metroname) + as.factor(regionname) + 
      year + f(ID,model="besag",graph="FOQ_INLA")
  }
  if(race_option == 1) {
    inla_formula4 <- deaths ~ 
      tbhk_pv + tbhk_is + tbhk_ce + fbp_kt + as.factor(metroname) + as.factor(regionname) + 
      year + f(ID,model="besag",graph="FOQ_INLA")
  }
message('Fitting INLA model...')
  inla_model_4 = inla(inla_formula4,
                      family = "binomial",
                      data = model_data_black,
                      Ntrials = model_data_black[, total_pop],
                      verbose = FALSE,
                      control.compute=list(config = TRUE, dic = TRUE),
                      control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
                      control.fixed=list(prec.intercept = 0,
                                         prec = 1))
## Make full prediction, full residual, and load posterior mean for all components.
  model_data_black[, inla_pred_4 := inla_model_4$summary.fitted.values$mean]
  model_data_black[, inla_residual_4 := logit(nmx) - logit(inla_pred_4)]
  model_data_black[nmx==0, inla_residual_4 := logit(nmx+0.000001) - logit(inla_pred_4)]
  coefs <- make_coef_table(4)

## Create permutations (2010-1990, 6 changes, total permutations = 2^6 = 64, 32 pairs)
## i.e. one pair for poverty is delta_m|PV=2013,IS=1990,CE=1990,FB=1990,time=1990,residual=1990 - 
##                              delta_m|PV=1990,IS=1990,CE=1990,FB=1990,time=1990,residual=1990
  make_permutations <- function(fes, start_year, end_year) {
    fes <- c(fes, 'year','residual')
    permutations <- list()
    for(fe in 1:length(fes)) permutations[[fe]] <- c(start_year,end_year)
    permutations <- as.data.table(expand.grid(permutations))
    setnames(permutations, names(permutations), fes)
    return(permutations)
  }
  permutations <- make_permutations(fes = fes,
                                    start_year = start_year,
                                    end_year = end_year)

## Prep and reshape input data from model (all fixed effects of interest + geographic random effects + time + spatial random effects + intercept + residual, wide by year)
  d <- copy(model_data_black)
  d <- merge(d, inla_model_4$summary.random$ID[c('ID','mean')], by='ID')
  setnames(d, 'mean', 'spatial_effect')
  d <- d[year %in% c(start_year,end_year), ]
  setnames(d, c('inla_residual_4','inla_pred_4'), c('inla_residual','inla_pred'))
  d <- dcast(d, fips ~ year, value.var = c(fes, 'inla_residual', 'inla_pred', 'total_pop','metroname','regionname','nmx','spatial_effect'))
  d[, (paste0('year_', start_year)) := start_year]
  d[, (paste0('year_', end_year)) := end_year]

## Predict on all permutations.
  calculate_contribution <- function(fe) {
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
          if(p_option[, residual]==end_year) this_dt[, p_residual := inla_residual_2010] ## MAKE SURE RESIDUAL IS IN LOGIT SPACE - it comes in normal space from INLA object.
          if(p_option[, residual]==start_year) this_dt[, p_residual := inla_residual_1990]
        }
      }
      
      ## Generate full prediction for this permutation based on whether FE value stayed the same or changed over the period.
      ## Assign target FE value based on change (2010 value) or no change (1990 value), and then add in all other effects.
      if(!(fe %in% c('residual'))) {
        this_dt[, (paste0('p_with_change_',p)) := (get(paste0(fe,'_2010')) * coefs[name==fe, coef])]
        this_dt[, (paste0('p_without_change_',p)) := (get(paste0(fe,'_1990')) * coefs[name==fe, coef])]
      }
      if(fe=='residual') {
        this_dt[, (paste0('p_with_change_',p)) := inla_residual_2010]
        this_dt[, (paste0('p_without_change_',p)) := inla_residual_1990]
      }
      
      ## Add intercept and random effects.
      for(r in unique(this_dt[, regionname_2010])) {
        if(r=='Pacific') this_dt[regionname_2010==r, region_int := 0]
        if(r!='Pacific') this_dt[regionname_2010==r, region_int := coefs[grep(r,name), coef]]
      }
      for(m in unique(this_dt[, metroname_2010])) {
        if(m=="Lg central metro") this_dt[metroname_2010==m, metro_int := 0]
        if(m!="Lg central metro") this_dt[metroname_2010==m, metro_int := coefs[grep(m,name), coef]]
      }
      this_dt[, (paste0('p_with_change_',p)) := get((paste0('p_with_change_',p))) + coefs[name=='(Intercept)', coef] + metro_int + region_int + spatial_effect_1990]
      this_dt[, (paste0('p_without_change_',p)) := get((paste0('p_without_change_',p))) + coefs[name=='(Intercept)', coef] +  + metro_int + region_int + spatial_effect_1990]
      
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
  
  ## Calculate contribution attributable to each time-varying component. By definition this adds up to total observed change in the outcome
  ## because of inclusion of the residual.
  all_contributions <- rbindlist(lapply(c(fes, 'year','residual'), calculate_contribution))
  
  ## Collapse contributions to metro-region.
  all_contributions <- merge(all_contributions, unique(d[, c('fips','metroname_2010','regionname_2010','total_pop_2010','total_pop_1990')]), by='fips')
  all_contributions <- all_contributions[, list(contribution_mort=wtd.mean(contribution_mort, total_pop_2010, na.rm=TRUE),
                                                total_pop_2010=sum(total_pop_2010, na.rm=TRUE)),
                                         by=c('metroname_2010','regionname_2010','fe')]
  all_contributions[, metro_region := paste0(metroname_2010, ' ', regionname_2010)]
  all_contributions <- merge(all_contributions, all_contributions[, list(total_change=sum(contribution_mort, na.rm=TRUE)), by='metro_region'], by='metro_region')
  all_contributions[, metro_region := factor(metro_region, levels=unique(all_contributions$metro_region[order(all_contributions[, total_change])]))]
  
  ## Format for plotting
  all_contributions[, fe := factor(fe, levels=c(fes,"year","residual"))]
  all_contributions <- merge(all_contributions, cov_names, by='fe')
  all_contributions[, cov_name := factor(cov_name, levels=cov_names[fe %in% all_contributions[, fe], cov_name])]
  
  return(all_contributions) 

}

metro_names <- c('Lg central metro', 'Lg fringe metro','Md/Sm metro','Nonmetro, adjacent','Nonmetro, nonadjacent')
region_names <- c('Pacific', 'Appalachia', "East South Central","Mountain","West South Central","New England","South Atlantic","East North Central","West North Central","Middle Atlantic")
cov_names <- data.table(fe = c(paste0('as.factor(year)',c(1990,2000,2010)),
                               paste0('as.factor(metroname)',metro_names),
                               paste0('as.factor(regionname)',region_names),
                               'air_EQI_22July2013','water_EQI_22July2013','land_EQI_22July2013',
                               'tbhk_pv','tbhk_is','tbhk_ce','fbp_kt','twhk_pv','twhk_is','twhk_ce','fwp_kt','metro',"Global Moran's I","DIC",'RMSE','year','residual'),
                        cov_name = c('1990','2000','2010', metro_names, region_names, 'Air quality','Water quality','Land quality',
                                     'Poverty','Isolation','College education',
                                     'Foreign-born','Poverty','Isolation','College education','Foreign-born', 'Metro level',"Global Moran's I","DIC",'RMSE','Secular trend','Residual'),
                        cov_sort = c(1:35))

all_contributions_white_male <- run_shapley_decomp(race_option = 0, sex_option = 1)
all_contributions_white_female <- run_shapley_decomp(race_option = 0, sex_option = 2)
all_contributions_black_male <- run_shapley_decomp(race_option = 1, sex_option = 1)
all_contributions_black_female <- run_shapley_decomp(race_option = 1, sex_option = 2)

pdf(paste0('C:/Users/ngraetz/Documents/Penn/papers/rwjf/decomp2.pdf'), width = 12, height = 8)
ggplot() +
  geom_bar(data=all_contributions_white_male,
           aes(x=metro_region,
               y=contribution_mort,
               fill=cov_name),
           color='black',
           stat='identity') + 
  geom_point(data=all_contributions_white_male[, list(contribution_mort=sum(contribution_mort)), by=c('metro_region')],
             aes(x=metro_region,
                 y=contribution_mort),
             size=3) + 
  labs(x = 'Metro/region', y = 'Change in ASDR, 1990-2010', title = 'Change in white male age-standardized mortality, 1990-2010') + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f'), name='Component')
ggplot() +
  geom_bar(data=all_contributions_white_female,
           aes(x=metro_region,
               y=contribution_mort,
               fill=cov_name),
           color='black',
           stat='identity') + 
  geom_point(data=all_contributions_white_female[, list(contribution_mort=sum(contribution_mort)), by=c('metro_region')],
             aes(x=metro_region,
                 y=contribution_mort),
             size=3) + 
  labs(x = 'Metro/region', y = 'Change in ASDR, 1990-2010', title = 'Change in white female age-standardized mortality, 1990-2010') + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f'), name='Component')
ggplot() +
  geom_bar(data=all_contributions_black_male,
           aes(x=metro_region,
               y=contribution_mort,
               fill=cov_name),
           color='black',
           stat='identity') + 
  geom_point(data=all_contributions_black_male[, list(contribution_mort=sum(contribution_mort)), by=c('metro_region')],
             aes(x=metro_region,
                 y=contribution_mort),
             size=3) + 
  labs(x = 'Metro/region', y = 'Change in deaths, 1990-2010', title = 'Change in black male age-standardized mortality, 1990-2010') + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f'), name='Component')
ggplot() +
  geom_bar(data=all_contributions_black_female,
           aes(x=metro_region,
               y=contribution_mort,
               fill=cov_name),
           color='black',
           stat='identity') + 
  geom_point(data=all_contributions_black_female[, list(contribution_mort=sum(contribution_mort)), by=c('metro_region')],
             aes(x=metro_region,
                 y=contribution_mort),
             size=3) + 
  labs(x = 'Metro/region', y = 'Change in deaths, 1990-2010', title = 'Change in black female age-standardized mortality, 1990-2010') + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f'), name='Component')
dev.off()









