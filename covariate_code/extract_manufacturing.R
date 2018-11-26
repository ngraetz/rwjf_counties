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

manu_2010 <- fread(paste0(repo, 'covariate_raw/factfinder_manufacturing_2010/ACS_10_5YR_S2403_with_ann_nolabels.csv'))
setnames(manu_2010, c('GEO.id2','HC01_EST_VC06','HC01_EST_VC01'), c('fips','manufacturing','total_persons'))
manu_2010[, manufacturing := manufacturing / total_persons]
manu_2010 <- manu_2010[, c('fips','manufacturing')]
manu_2010[, fips := as.character(fips)]
manu_2010[nchar(fips)==4, fips := paste0('0',fips)]
manu_2010[, year := 2010]

## (Male + female) / total = (VD07 + VD34) / VD01
manu_2000 <- fread(paste0(repo, 'covariate_raw/factfinder_manufacturing_2000/DEC_00_SF3_P049_with_ann_nolabels.csv'))
manu_2000[, fips := GEO.id2]
manu_2000[, manufacturing := (VD07 + VD34) / VD01]
manu_2000 <- manu_2000[, c('fips','manufacturing')]
manu_2000[, fips := as.character(fips)]
manu_2000[nchar(fips)==4, fips := paste0('0',fips)]
manu_2000[, year := 2000]

manu <- rbind(manu_2000, manu_2010, manu_2015)
write.csv(manu, paste0(repo, 'covariate_clean_data/factfinder_manufacturing.csv'), row.names = FALSE)
