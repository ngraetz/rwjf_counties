library(raster)
library(rgdal)
library(rgeos)
library(data.table)
library(maptools)
library(ggplot2)
library(spdep)
library(INLA)
library(spdep)
library(gridExtra)
library(grid)
library(GWmodel)
library(spgwr)
library(tables)
source('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/functions.R')

## Load shapefile
counties <- readOGR("C:/Users/ngraetz/Downloads/cb_2016_us_county_20m", 'cb_2016_us_county_20m')
counties@data <- transform(counties@data, fips = paste0(STATEFP, COUNTYFP)) # create unique county 5-digit fips
counties <- counties[counties@data$fips != "02016", ] # Drop Aleutians West, AK - screws up plots 
counties <- counties[counties$STATEFP != '02' &
                       counties$STATEFP != '15' &
                       counties$STATEFP != '72',]
background_counties <- counties[counties$STATEFP != '02' &
                                  counties$STATEFP != '15' &
                                  counties$STATEFP != '72',]
background.dt <- as.data.table(fortify(background_counties, region = 'STATEFP'))
map_colors <- c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')
counties$state <- as.numeric(counties$STATEFP)
states <- gUnaryUnion(counties, id = counties@data$state)
background.dt.states <- as.data.table(fortify(states, region = 'state'))

## Load covariates
covs <- c('eqi','exposures_by_income','exposures_by_nativity','pop_composition','typologies')
covs <- lapply(paste0('C:/Users/ngraetz/Downloads/',covs,'.csv'), fread)

## Format Brown exposure data
brown <- covs[[2]]
brown <- brown[, c('fips','twh1k_pv', 'tbh1k_pv', 'tbh1k_is', 'twh1k_is', 'twh1k_ce', 'tbh1k_ce',
                   'twh0k_pv', 'tbh0k_pv', 'tbh0k_is', 'twh0k_is', 'twh0k_ce', 'tbh0k_ce',
                   'twh9k_pv', 'tbh9k_pv', 'tbh9k_is', 'twh9k_is', 'twh9k_ce', 'tbh9k_ce')]
brown_pop <- covs[[3]]
brown_pop <- brown_pop[, c('fips', 'fb9p_kt', 'fb0p_kt', 'fb1p_kt',
                           'fw9p_kt', 'fw0p_kt', 'fw1p_kt')]
brown <- merge(brown, brown_pop, by='fips')
brown = melt(brown, id.vars = c("fips"),
             measure.vars = c('twh1k_pv', 'tbh1k_pv', 'tbh1k_is', 'twh1k_is', 'twh1k_ce', 'tbh1k_ce',
                              'twh0k_pv', 'tbh0k_pv', 'tbh0k_is', 'twh0k_is', 'twh0k_ce', 'tbh0k_ce',
                              'twh9k_pv', 'tbh9k_pv', 'tbh9k_is', 'twh9k_is', 'twh9k_ce', 'tbh9k_ce',
                              'fb9p_kt', 'fb0p_kt', 'fb1p_kt',
                              'fw9p_kt', 'fw0p_kt', 'fw1p_kt'))
brown[, variable := as.character(variable)]
brown[nchar(variable)==8, year := as.numeric(substr(variable,4,4))]
brown[nchar(variable)==7, year := as.numeric(substr(variable,3,3))]
brown[year==1, year := 2010]
brown[year==0, year := 2000]
brown[year==9, year := 1990]
brown[nchar(variable)==7, variable := paste0(substr(variable, 1, 2), substr(variable, 4, nchar(variable)))]
brown[nchar(variable)==8, variable := paste0(substr(variable, 1, 3), substr(variable, 5, nchar(variable)))]
brown <- dcast(brown, year + fips ~ variable, value.var = 'value')
cov_names <- c('twhk_pv', 'tbhk_pv', 'tbhk_is', 'twhk_is', 'twhk_ce', 'tbhk_ce', 'fwp_kt', 'fbp_kt')
for(c in cov_names) brown[, (c) := scales::rescale(get(c))]

titles <- c('A)','B)','C)','D)')
t <- 1
for(sex_option in c(1,2)) {
  t <- 1
  for(race_option in c(0,1)) {
    
    ## Set options
    year_as_factor <- TRUE
    #race_option <- 1
    #sex_option <- 1
    race_title <- ifelse(race_option==1, 'Black', 'White')
    pop_option <- 'all'
    
    ## Load nmx data
    if(pop_option=='all') nmx_list <- readRDS("C:/Users/ngraetz/Documents/repos/spatial_demography_2018/nmx_by_sex_nwh_black_county.rds")
    if(pop_option=='25_65') nmx_list <- readRDS("C:/Users/ngraetz/Documents/repos/spatial_demography_2018/data_byrace_25_65.rds")
    nmx <- nmx_list[[1]]
    
    if(sex_option==1) lims <- c(0,0.015)
    if(sex_option==2) lims <- c(0,0.012)
    
    nmx_map <- make_county_map(map_dt = nmx[year==2010 & race==race_option & sex==sex_option & fips %in% brown[!is.na(fbp_kt) & year==2010, fips],],
                               map_sp = counties,
                               map_var = 'nmx',
                               legend_title = 'Mortality\nrate',
                               high_is_good = FALSE,
                               map_title = titles[t],
                               map_limits = lims)
    assign(paste0('map_', race_option, '_', sex_option), nmx_map)
    t <- t + 1
      
  }
}

gLegend<-function(a.plot){
  
  if ("ggplot" %in% class(a.plot)) {
    tmp <- ggplot_gtable(ggplot_build(a.plot))
  } else if ("grob" %in% class(a.plot)) {
    tmp <- .gplot
  }
  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#png(paste0('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/nmx_2010.png'), width = 2400, height = 1400, res = 120)
pdf(paste0('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/male_nmx_2010.pdf'), width = 13, height = 5)
# Initialize plot with master title
p.legend <- gLegend(map_0_1)
p.legend$vp <- viewport(layout.pos.row = 1:6, layout.pos.col = 11)
grid.newpage()
pushViewport(viewport(layout = grid.layout(6, 11)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y, clip = "off")
# Plot all data coverage maps
#print(tbl, as.table=TRUE)
print(map_0_1 + theme(legend.position="none"), vp = vplayout(1:6, 1:5))
print(map_1_1 + theme(legend.position="none"), vp = vplayout(1:6, 6:10))
# print(map_1_1 + theme(legend.position="none"), vp = vplayout(7:12, 1:5))
# print(map_1_2 + theme(legend.position="none"), vp = vplayout(7:12, 6:10))
grid.draw(p.legend)
dev.off()
#png(paste0('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/nmx_2010.png'), width = 2400, height = 1400, res = 120)
pdf(paste0('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/female_nmx_2010.pdf'), width = 13, height = 5)
# Initialize plot with master title
p.legend <- gLegend(map_0_2)
p.legend$vp <- viewport(layout.pos.row = 1:6, layout.pos.col = 11)
grid.newpage()
pushViewport(viewport(layout = grid.layout(6, 11)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y, clip = "off")
# Plot all data coverage maps
#print(tbl, as.table=TRUE)
print(map_0_2 + theme(legend.position="none"), vp = vplayout(1:6, 1:5))
print(map_1_2 + theme(legend.position="none"), vp = vplayout(1:6, 6:10))
# print(map_1_1 + theme(legend.position="none"), vp = vplayout(7:12, 1:5))
# print(map_1_2 + theme(legend.position="none"), vp = vplayout(7:12, 6:10))
grid.draw(p.legend)
dev.off()




make_county_map <- function(map_dt, map_sp, map_var, high_is_good, map_title, map_limits=NULL, diverge=FALSE, diverge_point=0, legend_title, manual_colors = NULL) {
  map_colors <- c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')
  sp_map <- merge(map_sp, map_dt, by='fips')
  sp_fort <- fortify(sp_map, region='fips')
  sp_fort <- as.data.table(merge(sp_fort, sp_map@data, by.x = "id", by.y = "fips"))
  if(high_is_good==FALSE) map_colors <- rev(map_colors)
  # Cap data at limits
  if(!is.null(map_limits)) {
    sp_fort[get(map_var) < min(map_limits), (map_var) := min(map_limits)]
    sp_fort[get(map_var) > max(map_limits), (map_var) := max(map_limits)]
  }
  this_gg <- ggplot() +
    # geom_polygon(data = background.dt,
    #              aes(x = long,
    #                  y = lat,
    #                  group = group),
    #              fill = "grey",
    #              alpha = 0.5) + 
    geom_polygon(data = sp_fort[!is.na(get(map_var)),],
                 aes(x = long,
                     y = lat,
                     group = group,
                     fill = get(map_var)),
                 color = 'black',
                 size = 0.1) +
    geom_path(data = background.dt.states,
              aes(x = long,
                  y = lat,
                  group = group),
              color = 'black',
              size = 0.2) +
    theme_classic() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    ggtitle(map_title)
  #message(map_colors)
  if(is.null(map_limits) & is.null(manual_colors)) this_gg <- this_gg + scale_fill_gradientn(colors = map_colors, guide = guide_legend(title = legend_title))
  if(!is.null(map_limits) & is.null(manual_colors)) this_gg <- this_gg + scale_fill_gradientn(colors = map_colors, limits = map_limits, na.value = map_colors[length(map_colors)], guide = guide_legend(title = legend_title))
  if(diverge) {
    values <- c(min(map_limits),diverge_point,max(map_limits))
    redblue <- c('#08519c','#ffffff','#67000d')
    #this_gg <- this_gg + scale_fill_gradientn(colours=redblue, values=values, na.value = "#000000", rescaler = function(x, ...) x, oob = identity)
    this_gg <- this_gg + scale_fill_gradientn(guide = guide_legend(title = legend_title), limits = map_limits, colours=redblue, values=values, na.value = "#000000", rescaler = function(x, ...) x, oob = identity)
  }
  if(!is.null(manual_colors)) this_gg <- this_gg + scale_fill_manual(values = manual_colors, guide = guide_legend(title = legend_title))
  return(this_gg)
}

titles <- c('A)','B)','C)','D)','E)','F)','G)','H)')
t <- 1
for(mv in c('pov','edu','iso','fb')) {
  t <- 1
for(race_option in c(0,1)) {

    if(mv=='fb') {
      black_var <- 'fbp_kt'
      white_var <- 'fwp_kt'
      ml <- c(0,0.2)
      good <- TRUE
    }
    if(mv=='pov') {
      black_var <- 'tbhk_pv'
      white_var <- 'twhk_pv'
      ml <- c(0,0.4)
      good <- FALSE
    }
    if(mv=='edu') {
      black_var <- 'tbhk_ce'
      white_var <- 'twhk_ce'
      ml <- c(0,0.4)
      good <- TRUE
    }
    if(mv=='iso') {
      black_var <- 'tbhk_is'
      white_var <- 'twhk_is'
      ml <- c(0,0.8)
      good <- FALSE
    }
    nmx_map <- make_county_map(map_dt = brown[year==2010,],
                               map_sp = counties,
                               map_var = ifelse(race_option==1, black_var, white_var),
                               legend_title = 'Proportion',
                               high_is_good = good,
                               map_title = titles[t],
                               map_limits = ml)
    assign(paste0('map_', race_option, '_', mv), nmx_map)
    t <- t + 1
}
}
pdf(paste0('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/fb_2010.pdf'), width = 13, height = 5)
# Initialize plot with master title
pov.legend <- gLegend(map_1_fb)
pov.legend$vp <- viewport(layout.pos.row = 1:6, layout.pos.col = 11)
# edu.legend <- gLegend(map_1_edu)
# edu.legend$vp <- viewport(layout.pos.row = 7:12, layout.pos.col = 11)
# iso.legend <- gLegend(map_1_iso)
# iso.legend$vp <- viewport(layout.pos.row = 13:18, layout.pos.col = 11)
# fb.legend <- gLegend(map_1_fb)
# fb.legend$vp <- viewport(layout.pos.row = 19:24, layout.pos.col = 11)
grid.newpage()
pushViewport(viewport(layout = grid.layout(6, 11)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y, clip = "off")
print(map_0_fb + theme(legend.position="none"), vp = vplayout(1:6, 1:5))
print(map_1_fb + theme(legend.position="none"), vp = vplayout(1:6, 6:10))
# print(map_0_edu + theme(legend.position="none"), vp = vplayout(7:12, 1:5))
# print(map_1_edu + theme(legend.position="none"), vp = vplayout(7:12, 6:10))
# print(map_0_iso + theme(legend.position="none"), vp = vplayout(13:18, 1:5))
# print(map_1_iso + theme(legend.position="none"), vp = vplayout(13:18, 6:10))
# print(map_0_fb + theme(legend.position="none"), vp = vplayout(19:24, 1:5))
# print(map_1_fb + theme(legend.position="none"), vp = vplayout(19:24, 6:10))
grid.draw(pov.legend)
# grid.draw(edu.legend)
# grid.draw(iso.legend)
# grid.draw(fb.legend)
dev.off()

pdf(paste0('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/edu_2010.pdf'), width = 13, height = 5)
# Initialize plot with master title
pov.legend <- gLegend(map_1_edu)
pov.legend$vp <- viewport(layout.pos.row = 1:6, layout.pos.col = 11)
grid.newpage()
pushViewport(viewport(layout = grid.layout(6, 11)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y, clip = "off")
print(map_0_edu + theme(legend.position="none"), vp = vplayout(1:6, 1:5))
print(map_1_edu + theme(legend.position="none"), vp = vplayout(1:6, 6:10))
grid.draw(pov.legend)
dev.off()

pdf(paste0('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/iso_2010.pdf'), width = 13, height = 5)
# Initialize plot with master title
pov.legend <- gLegend(map_1_iso)
pov.legend$vp <- viewport(layout.pos.row = 1:6, layout.pos.col = 11)
grid.newpage()
pushViewport(viewport(layout = grid.layout(6, 11)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y, clip = "off")
print(map_0_iso + theme(legend.position="none"), vp = vplayout(1:6, 1:5))
print(map_1_iso + theme(legend.position="none"), vp = vplayout(1:6, 6:10))
grid.draw(pov.legend)
dev.off()




########## LISA COUNTIES
black_male <- readRDS("C:/Users/ngraetz/Documents/repos/spatial_demography_2018/Black_all/plots_fixed_region/sex_1_lisa_counties.rds")
black_female <- readRDS("C:/Users/ngraetz/Documents/repos/spatial_demography_2018/Black_all/plots_fixed_region/sex_2_lisa_counties.rds")
white_male <- readRDS("C:/Users/ngraetz/Documents/repos/spatial_demography_2018/White_all/plots_fixed_region/sex_1_lisa_counties.rds")
white_female <- readRDS("C:/Users/ngraetz/Documents/repos/spatial_demography_2018/White_all/plots_fixed_region/sex_2_lisa_counties.rds")
black_male_counties <- black_male[morans_col_cat=='High-High', fips]
black_female_counties <- black_female[morans_col_cat=='High-High', fips]
white_male_counties <- white_male[morans_col_cat=='High-High', fips]
white_female_counties <- white_female[morans_col_cat=='High-High', fips]
iso_fips <- black_male[tbhk_is >= 0.75 & year==2010, fips]
sum(black_male[morans_col_cat=='High-High', total_pop]) / sum(black_male[, total_pop])
make_county_map(map_dt = black_male[fips %in% iso_fips, ],
                map_sp = counties,
                map_var = 'fbp_kt',
                legend_title = 'Proportion',
                high_is_good = FALSE,
                map_title = 'prop')
saveRDS(iso_fips, 'C:/Users/ngraetz/Documents/repos/spatial_demography_2018/Black_all/plots_fixed_region/black_iso_fips.rds')

lisa_lts <- fread('C:/Users/ngraetz/Documents/repos/spatial_demography_2018/lisa_regions.csv')
lisa_lts <- lisa_lts[Age.x == 0, ]
lisa_lts[, region := paste0(region,'_',sex)]
ggplot() +
  geom_line(data = lisa_lts,
            aes(x = year,
                y = ex,
                color = region),
            size = 2) +
  theme_minimal()