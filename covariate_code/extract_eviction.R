
e <- fread('https://eviction-lab-data-downloads.s3.amazonaws.com/CA/counties.csv')
e <- e[, c('GEOID','year','eviction-rate')]
