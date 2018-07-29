library(readxl)
library(data.table)
library(ggplot2)

## Pull directly from BLS website and format (could save these as flat files somewhere).
pull_from_bls_website <- function(y) {
  tmp = tempfile(fileext = ".xlsx")
  download.file(url = paste0("https://www.bls.gov/lau/laucnty", y, ".xlsx"), destfile = tmp, mode="wb")
  d <- as.data.table(read_excel(tmp, skip = 4, col_names = TRUE))
  setnames(d, names(d), c('laus_code','state_fips','county_fips','county_name','year','skip','labor_force','employed','unemployed','percent_unemployment'))
  d[, fips := paste0(state_fips, county_fips)]
  d <- d[, c('year','fips','labor_force','employed','unemployed','percent_unemployment')]
  return(d)
}

file_years <- c(as.character(90:99), paste0('0', as.character(0:9)), as.character(10:17))
d <- rbindlist(lapply(file_years , pull_from_bls_website))

## Save.
write.csv(d, 'C:/Users/ngraetz/Documents/Penn/papers/rwjf/covariates/bls_laus_covs.csv', row.names=FALSE)


