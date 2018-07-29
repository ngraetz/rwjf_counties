library(INLA)
library(data.table)
inv_logit <- function(x) {
  exp(x)/(1+exp(x))
}
logit <- function(x) {
  log(x/(1 - x))
}

race_option <- 0
sex_option <- 1

start_year <- 1990
end_year <- 2010

if(race_option==0) fes <- c('twhk_pv','twhk_is','twhk_ce','fwp_kt')
if(race_option==1) fes <- c('tbhk_pv','tbhk_is','tbhk_ce','fbp_kt')

model_data_black <- readRDS('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/model_data_black.RDS')

# Create a neighborhood file
nb.FOQ <- poly2nb(counties, queen=TRUE)
# Create an INLA weight matrix
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
## Fit model and extract parts that vary over time.
inla_model_4 = inla(inla_formula4,
                    family = "binomial",
                    data = model_data_black,
                    Ntrials = model_data_black[, total_pop],
                    verbose = FALSE,
                    control.compute=list(config = TRUE, dic = TRUE),
                    control.inla=list(int.strategy='eb', h = 1e-3, tolerance = 1e-6),
                    control.fixed=list(prec.intercept = 0,
                                       prec = 1))
model_data_black[, inla_pred_4 := inla_model_4$summary.fitted.values$mean]
model_data_black[, inla_residual_4 := logit(nmx) - logit(inla_pred_4)]
coefs <- make_coef_table(4)

## Calculate total change and decompose to contribution of each time-varying component.
d <- copy(model_data_black)
d <- merge(d, inla_model_4$summary.random$ID[c('ID','mean')], by='ID')
setnames(d, 'mean', 'spatial_effect')
d <- d[year %in% c(start_year,end_year), ]
setnames(d, c('inla_residual_4','inla_pred_4'), c('inla_residual','inla_pred'))
d <- dcast(d, fips ~ year, value.var = c(fes, 'inla_residual', 'inla_pred', 'total_pop','metroname','regionname','nmx','spatial_effect'))
d[, (paste0('year_', start_year)) := start_year]
d[, (paste0('year_', end_year)) := end_year]
d <- d[!is.na(metroname_1990) & !is.na(regionname_1990)]

## Make sure preds line up 
for(r in unique(d[, regionname_2010])) {
  if(r=='Pacific') d[regionname_2010==r, region_int := 0]
  if(r!='Pacific') d[regionname_2010==r, region_int := coefs[grep(r,name), coef]]
}
for(m in unique(d[, metroname_2010])) {
  if(m=="Lg central metro") d[metroname_2010==m, metro_int := 0]
  if(m!="Lg central metro") d[metroname_2010==m, metro_int := coefs[grep(m,name), coef]]
}
d[, pred_1990 := coefs[name=='(Intercept)', coef] + region_int + metro_int + spatial_effect_1990]
for(fe in c(fes,'year')) d[, pred_1990 := pred_1990 + get(paste0(fe, '_1990')) * coefs[name==fe, coef]]
d[, pred_1990 := pred_1990 + inla_residual_1990]
d[, pred_1990 := inv_logit(pred_1990)]

## Calculate changes in normal space: logit(2010) - logit(1990)


## Calculate change attributable to each time-varying linear predictor.
d[, c_residual := get(paste0('inla_residual_', end_year)) - get(paste0('inla_residual_', start_year))]
for(fe in c(fes,'year')) {
  ## Add each time-varying linear predictor for the first and last year of change calculated.
  for(y in c(start_year, end_year)) d[, (paste0('pred_', fe, '_', y)) := get(paste0(fe, '_', y)) * coefs[name==fe, coef]]
  ## Calculate change based on just that predictor * coefficient.
  d[, (paste0('c_', fe)) := get(paste0('pred_', fe, '_', end_year)) - get(paste0('pred_', fe, '_', start_year))]
}


d[, observed_change := logit(get(paste0('nmx_', end_year))) - logit(get(paste0('nmx_', start_year)))]
d[, inla_pred_change := logit(inla_pred_2010) - logit(inla_pred_1990)]

d[, pred_change := c_residual + c_year]
for(fe in c(fes,'year')) d[, pred_change := pred_change + get(paste0('c_', fe))]


head(inla_model_4$summary.random$ID[c('ID','mean')])


