library(data.table)
library(foreign)

## Looks like this is only by race, not sex.
meta <- as.data.table(read.dbf("C:/Users/ngraetz/Downloads/all/tables.dbf"))
meta[grep('EDU',TEXT)]
meta[TEXT=='RACE BY EDUCATIONAL ATTAINMENT', which=TRUE]
meta[1322:1350,]

meta[grep('INDUSTRY',TEXT)]
meta[TEXT=='INDUSTRY', which=TRUE]
meta[1978:2000,]


# 12:  P058       NA                                      White:      NA
# 13:  P058 P0580001                         Less than 9th grade  STF310
# 14:  P058 P0580002               9th to 12th grade, no diploma  STF310
# 15:  P058 P0580003 High school graduate (includes equivalency)  STF310
# 16:  P058 P0580004                     Some college, no degree  STF310
# 17:  P058 P0580005                            Associate degree  STF310
# 18:  P058 P0580006                           Bachelor's degree  STF310
# 19:  P058 P0580007             Graduate or professional degree  STF310

## Function to find correct table (STF310) path on website.
find_table_file <- function(site) {
  thepage = readLines(site)
  thepage <- thepage[grep('stf310', thepage)]
  thepage <- strsplit(thepage, split = 'href=\"')
  thepage <- strsplit(thepage[[1]][2], split = '\">')
  file <- thepage[[1]][1]
  return(file)
}

## Function to read correct table and aggregate education variables for White.
pull_edu_table <- function(site) {
  site <- paste0('https://www2.census.gov/census_1990/CD90_3A_', site, '/')
  file <- find_table_file(site)
  temp <- tempfile()
  download.file(paste0(site, file),temp)
  target_table <- as.data.table(read.dbf(temp))
  unlink(temp)
  target_table[, fips := paste0(as.character(STATEFP),as.character(CNTY))]
  ## Construct total college, total less than high school.
  target_table[, college := P0580006 + P0580007]
  target_table[, less_12 := P0580001 + P0580002]
  target_table[, list(less_12=sum(less_12, na.rm=TRUE), college=sum(college, na.rm=TRUE)), by='fips']
  return(target_table)
}

## There is a separate numbered webpage per chunk of counties (01-61).
all_sites <- as.character(c(1:61))
all_sites[nchar(all_sites)==1] <- paste0('0',all_sites[nchar(all_sites)==1])
all_sites[all_sites=='09'] <- '9' ## For some reason every site has a leading 0 except 9...

## Pull for each chunk of counties and append.
all_edu <- rbindlist(lapply(all_sites, pull_edu_table))




