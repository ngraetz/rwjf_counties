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
  if(y != '') d[, (v) := paste0(get(v), '_', y)]
}

