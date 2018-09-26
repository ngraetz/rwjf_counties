library(data.table)

chr <- fread("C:/Users/ngraetz/Downloads/CHR_TRENDS_CSV_2017.csv")
chr[, statecode := as.character(statecode)]
chr[, countycode := as.character(countycode)]
chr[nchar(statecode)==1, statecode := paste0('0',statecode)]
chr[nchar(countycode)==1, countycode := paste0('00',countycode)]
chr[nchar(countycode)==2, countycode := paste0('0',countycode)]
chr[, fips := paste0(statecode, countycode)]
chr <- chr[nchar(fips) == 5, ]
chr <- dcast(chr, yearspan + fips ~ measurename, value.var='rawvalue')
#chr <- chr[yearspan %in% c('2012-2014','2003-2005'), ]
# chr[yearspan == '2012-2014', year := 2013]
# chr[yearspan == '2003-2005', year := 2004]
setnames(chr, c('Adult obesity','Mammography screening',"Diabetes monitoring"), c('chr_obesity_prev','chr_mammography','chr_diabetes_monitoring'))
chr <- chr[, c('fips','yearspan','chr_obesity_prev','chr_mammography','chr_diabetes_monitoring')]
chr[, chr_obesity_prev := as.numeric(chr_obesity_prev)]
chr[, chr_mammography := as.numeric(chr_mammography)]
chr[, chr_diabetes_monitoring := as.numeric(chr_diabetes_monitoring)]
## Clean up years by variable
chr_obesity <- chr[!is.na(chr_obesity_prev), ]
chr_obesity[, year := as.numeric(substr(yearspan, 1, 4)) + 1]
chr_obesity <- chr_obesity[, c('fips','year','chr_obesity_prev')]
chr_mammography <- chr[!is.na(chr_mammography), ]
chr_mammography[yearspan == '2006-2007', yearspan := '2006']
chr_mammography[, year := as.numeric(yearspan)]
add <- chr_mammography[year==2006, ]
add[, year := 2007]
chr_mammography <- rbind(chr_mammography, add)
chr_mammography <- chr_mammography[, c('fips','year','chr_mammography')]
chr_diabetes_monitoring <- chr[!is.na(chr_diabetes_monitoring), ]
chr_diabetes_monitoring[yearspan == '2006-2007', yearspan := '2006']
chr_diabetes_monitoring[, year := as.numeric(yearspan)]
add <- chr_diabetes_monitoring[year==2006, ]
add[, year := 2007]
chr_diabetes_monitoring <- rbind(chr_diabetes_monitoring, add)
chr_diabetes_monitoring <- chr_diabetes_monitoring[, c('fips','year','chr_diabetes_monitoring')]
chr <- merge(chr_obesity, chr_mammography, all.x=TRUE, all.y=TRUE, by=c('fips','year'))
chr <- merge(chr, chr_diabetes_monitoring, all.x=TRUE, all.y=TRUE, by=c('fips','year'))
write.csv(chr, "C:/Users/ngraetz/Documents/repos/rwjf_counties/covariate_clean_data/chr_covs.csv", row.names=FALSE)

