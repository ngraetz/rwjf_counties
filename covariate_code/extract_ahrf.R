library(data.table)
library(ggplot2)

## Load data (for some absolutely insane reason, these data are stored as a single block of text, "V1" here).
d <- as.data.table(read.delim("C:/Users/ngraetz/Downloads/AHRF_2016-2017/DATA/ahrf2017.asc", fill=TRUE, header = FALSE))
d[, V1 := as.character(V1)]

## Load map from technical documentation to parse columns correctly.
col_map <- fread("C:/Users/ngraetz/Documents/repos/rwjf/AHRF_column_parsing_map.csv")
col_map <- col_map[col_col != 'COL-COL' & !is.na(col_col) & col_col != '', ]

## Parse columns according to "COL-COL" mapping of start- and end-characters in the technical documentation.
for(c in unique(col_map[, col_col])) {
  ## Grab variable name for this column from map.
  v <- col_map[col_col == c, variable_name]
  message(paste0('Parsing ', v, '...'))
  ## Grab year of this indicator from map, if it exists.
  y <- col_map[col_col == c, year]
  ## Parse characters according to start and end character from "COL-COL" variable in map.
  c_split <- as.numeric(unlist(strsplit(c, '-')))
  d[, (v) := substr(V1, c_split[1], c_split[2])]
  ## Add "year" to the variable name if it is year-specific (data are stored wide on year).
  #if(y != '') d[, (v) := paste0(get(v), '_', y)]
  if(y != '') setnames(d, v, paste0(v, '_', y)) 
}

## Format and save only variables we want for now. 
d[, fips := paste0(`FIPS State Code`, `FIPS County Code`)]
d_long <- melt(d, id.vars = c('fips','County Name'), measure.vars = names(d)[grep('Total Active M.D.s Non-Federal',names(d))], value.name = 'total_mds')
d_long[, variable := as.character(variable)]
d_long[, year := as.numeric(gsub('Total Active M.D.s Non-Federal_', '', variable))]
d_long[, new_total_mds := as.numeric(total_mds)]

# all <- counties$fips
# all_missing_fips <- all[!(all %in% d_long[, fips])]
# for(f in all_missing_fips) {
#   this_name <- template[fips==f, NAME]
#   this_state <- substr(f,1,2)
#   message(paste0(f, ' missing: ', this_name))
#   this_state <- d_codes[grep(paste0('^',this_state), fips),]
#   possible_matches <- this_state[grep(this_name, GeoName), c('fips','GeoName')]
#   for(m in possible_matches[, fips]) message(paste0('          ',m,' in BEA: ', possible_matches[fips==m, GeoName]))
# }

m <- make_county_map(map_dt = d_long[year==1990, ],
                     map_sp = counties,
                     map_var = 'new_total_mds',
                     legend_title = 'new_total_mds',
                     high_is_good = TRUE,
                     map_title = '1990 Transfers, %',
                     map_limits = c(0,100)) + guides(fill=guide_legend(title="Transfers, %"))
print(m)

d_long <- d_long[, c('fips','year','new_total_mds')]
write.csv(d_long, paste0(repo, 'covariate_clean_data/ahrf_covs.csv'), row.names=FALSE)




