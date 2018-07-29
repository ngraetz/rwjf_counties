
plot_sig_clusters <- function(d,map_var,year_var,map_year,plot_dir,title_add) {

## Descriptive spatial 
shape.shp <- merge(counties, d[get(year_var)==map_year, c('fips', map_var), with=FALSE], by='fips')
shape.shp$varofint <- shape.shp@data[[map_var]]
#attach(shape.shp@data)

message(paste0('calculating LISA for ', map_var))

# Create neighborhood weight matrices
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)
summary(nb.FOQ)
nb.FOQ.cor <- nb.FOQ

# Create a nearest neighbor weigh matrix
# Coordinates of UTM-projected shapefile
id <- as.numeric(paste(shape.shp$fips))

coords <- coordinates(shape.shp)
# Coordinates from reprojected Shapefile with LotLang-projection 
coords.ll <- coordinates(shape.shp)
l.5NN <- knearneigh(coords.ll, k=5, longlat=T)
nb.5NN <- knn2nb(l.5NN, row.names=id)
#plot.nb(nb.5NN, coords)

# Create a distance based matrix (distance between centroids)
# Distance of 20 km
d20km <- dnearneigh(coords.ll, 0, 20, row.names = id, longlat=T)
# Distance of 50 km
d50km <- dnearneigh(coords.ll, 0, 50, row.names = id,longlat=T)

# Which variables are available in the shape.file attribute table
names(shape.shp)

# Define variable of interest
varofint <- shape.shp$varofint

# Define name of the variable of interest 
varofint.name <- map_var

# Define neighborhood matrix (type nb) 
# (choice in this case nb.FOQ.cor, d50km, nb.5NN)
nb <- nb.5NN

# Define weight style (W=row-standardized)
ws <- c("W")

# Define significance level for the cluster maps 
# (0.0001,0.001,0.01,0.05)
significance <- 0.05

# p-adjustment method (can be "none", "bonferroni", "holm",
# "hochberg","hommel","fdr")
p.ad.meth <- c("none")

# Should the cluster map show only regions belonging to significent clusters, 
# or all regions
plot.only.significant <- TRUE

# Transfer nb-object in listwise object
listw <- nb2listw(nb, style=ws, zero.policy=TRUE)

# Create lagged values of variable of interest
varofint[is.na(varofint)] <- 0
varlag <- lag.listw(listw, varofint)

# Map of Variable and Lagged Variable
png(paste0(plot_dir, 'spatially_lagged.png'), width = 1200, height = 800, res = 120)
par(mfrow=c(1,2))                                      
var.data <- data.frame(varofint, varlag)
var.names <-c(varofint.name, paste("Spatially Lagged\n", title_add))
m <- length(var.data)
n <- length(id)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
  varofint1 <- sort(var.data[,i])
  mean <- mean(varofint1)
  sd <- sd(varofint1)
  min <- min(varofint1)          
  max <- max(varofint1)
  #bins <- c(min,mean-sd, mean, mean+sd, max)
  bins <- quantile(varofint, probs = c(0,.25,.5,.75,1))
  lb <- c(length(bins)-1)
  colpal <- brewer.pal(length(bins-1), "PiYG")
  colors <- colpal[findInterval(varofint1, bins, rightmost.closed=T)]
  plot(shape.shp, col=colors)
  title(var.names[[i]], cex=1)
  legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
                                                "-", round(bins[-1],2)),cex=0.8)
}
dev.off()

# Calculate Lisa Test
lisa.FOQ <- localmoran(varofint,listw, alternative="two.sided",
                       p.adjust.method=p.ad.meth)


# Get significance level
n <- length(id)
vec <- c(1:n)
vec <- ifelse(lisa.FOQ[,5] < significance, 1,0)

m.varofint <- mean(varofint)
m.varlag <- mean(varlag, na.rm=TRUE)

# Derive sector
varlag[is.na(varlag)] <- 0
sec <- c(1:n)
for (i in 1:n) {
  if (varofint[[i]]>=m.varofint & varlag[[i]]>=m.varlag) sec[i] <- 1
  if (varofint[[i]]<m.varofint & varlag[[i]]<m.varlag) sec[i] <- 2
  if (varofint[[i]]<m.varofint & varlag[[i]]>=m.varlag) sec[i] <- 3
  if (varofint[[i]]>=m.varofint & varlag[[i]]<m.varlag) sec[i] <- 4
}

# Define colors for sectors
sec.all <- sec
colors1 <- c(1:n)
for (i in 1:n) {
  if (sec.all[i]==1) colors1[i] <- "brown2"
  if (sec.all[i]==2) colors1[i] <- "royalblue3"
  if (sec.all[i]==3) colors1[i] <- "lightblue"
  if (sec.all[i]==4) colors1[i] <- "pink"
  if (sec.all[i]==0) colors1[i] <- "white"
}

# Mark all non-significant regions white
loc.m.data <- sec*vec
colors2 <- colors1
loc.m.data[is.na(loc.m.data)] <- 0
for (i in 1:n) {
  if (loc.m.data[i]==0) colors2[i] <- "white"
}

# Cluster map
par(mfrow=c(1,1),mar=c(5, 4, 4, 2) + 0.1)
if (plot.only.significant==TRUE) 
{plot(shape.shp, col=colors2, border="grey25",lwd=0.7)} else 
{plot(shape.shp, col=colors1, border="grey25",lwd=0.7)}
legend("bottomright",fill=c("brown2","royalblue3","lightblue",
                            "pink","white"),
       legend=c("High-High","Low-Low","Low-High","High-Low"),
       border = "grey25", cex=0.8, bg="white", 
       box.col="white")
title(paste("Significant Clusters\n",title_add))

png(paste0(plot_dir, 'significance_map.png'), width = 1200, height = 800, res = 120)
# Significance map
sig.data <- data.frame(lisa.FOQ[,5])
m <- length(sig.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
  signific <- sig.data[,i]
  lb <- 5
  bins <- c(0,0.0001,0.001,0.01,0.05,1)
  colpal <- c(rev(brewer.pal(lb+1, "YlGn")[-c(1:2)]),"white")
  colors[,i] <- colpal[findInterval(signific, bins, rightmost.closed=T)]
  colors[vec==0,i] <- c("white")
  plot(shape.shp, col=colors[,i], border="grey25", lwd=0.7)
  title(paste0("Significance Level: ", title_add))    
  colpalad <- colpal
  colpalad[c(which(bins==significance):length(colpalad))] <- c("white")
  binsp <- c("0","0.0001","0.001","0.01","0.05",1)
  legend("bottomright",fill=colpalad,
         legend=paste(binsp[-length(bins)],"-",binsp[-1]), border="grey25",
         cex=0.8, bg="white", box.col="white")
}
dev.off()

}

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
    geom_polygon(data = background.dt,
                 aes(x = long,
                     y = lat,
                     group = group),
                 fill = "grey",
                 alpha = 0.5) + 
    geom_polygon(data = sp_fort[!is.na(get(map_var)),],
                 aes(x = long,
                     y = lat,
                     group = group,
                     fill = get(map_var)),
                 color = 'black',
                 size = 0.2) +
    geom_path(data = background.dt.states,
              aes(x = long,
                  y = lat,
                  group = group),
              color = 'black',
              size = 1) +
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

plot_inla_betas <- function(res_fit) {
  model_data <- as.data.table(res_fit$summary.fixed)
  names(model_data)[names(model_data)=='0.025quant'] <- 'lower'
  names(model_data)[names(model_data)=='0.975quant'] <- 'upper'
  model_data[, cov := row.names(res_fit$summary.fixed)]
  model_data <- model_data[cov!='(Intercept)']
  model_data <- as.data.table(model_data)
  cov_names <- data.table(cov = c('air_EQI_22July2013','water_EQI_22July2013','land_EQI_22July2013',
                                  'tbhk_pv','tbhk_is','tbhk_ce','fbp_kt','twhk_pv','twhk_is','twhk_ce','fwp_kt','metro','year'),
                          cov_name = c('Air quality','Water quality','Land quality',
                                       'Poverty','Isolation','College education','Foreign-born','Poverty','Isolation','College education','Foreign-born',
                                       'Metro level','Year'),
                          cov_sort = c(1:13))
  model_data <- merge(model_data, cov_names, by='cov')
  model_data[, name := factor(cov_name, levels = cov_name[order(cov_sort)])]
  selected_gg <- ggplot(model_data, aes(x=name, y=mean)) +
    geom_point(size=3) +
    geom_hline(yintercept=0, color='red') +
    geom_errorbar(aes(ymax = upper, ymin = lower), width=0.25) +
    theme(axis.text.x = element_text(angle = 75, hjust = 1, size=15)) +
    theme_minimal() +
    theme(strip.text.x = element_text(size = 15),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    coord_flip() 
  return(selected_gg)
}

plot_lisa <- function(lisa_var,lisa_var_name,lisa_dt,lisa_sp,lisa_id,sig=1,n_neighbors=5,matrix='nn') {
  # Set up neighbors
  test.sp <- merge(lisa_sp, lisa_dt, by=lisa_id)
  test.sp <- test.sp[!is.na(test.sp@data[[lisa_var]]), ]
  id <- as.numeric(paste(test.sp[[lisa_id]]))
  coords.ll <- coordinates(test.sp)
  if(matrix == 'queen') {
    nb.FOQ <- poly2nb(test.sp, queen=TRUE)
    nb.lw.w <- nb2listw(nb.FOQ, zero.policy = TRUE)
  }
  if(matrix == 'nn') {
    l.5NN <- knearneigh(coords.ll, k=n_neighbors, longlat=T)
    nb.5NN <- knn2nb(l.5NN, row.names=id)
    nb.lw.w <- nb2listw(nb.5NN, zero.policy = TRUE)
  }
  # Calculate Lisa Test (local Moran's I)
  lisa.FOQ <- localmoran(test.sp[[lisa_var]],
                         nb.lw.w, alternative="two.sided",
                         p.adjust.method='none')
  # Calculate global Moran's I
  mt.res <- moran.test(test.sp[[lisa_var]], nb.lw.w, zero.policy = TRUE)
  # Map local Moran's I
  test.sp$morans <- lisa.FOQ[,1]
  test.sp$morans_sig <- lisa.FOQ[,5]
  test.dt <- as.data.table(test.sp@data)
  
  ## Create lagged var for categorizing 
  varlag <- lag.listw(nb.lw.w, test.sp[[lisa_var]])
  m.varofint <- mean(test.sp[[lisa_var]])
  m.varlag <- mean(varlag)
  test.dt[, varlag := varlag]
  #quantNorm = function(x){qnorm(rank(x,ties.method = "average")/(length(x)+1))}
  #test.dt[, morans_sig_normal := quantNorm(1-morans_sig)]
  test.dt[morans_sig <= 1, morans_sig_cat := 1]
  test.dt[morans_sig <= 0.40, morans_sig_cat := 2]
  test.dt[morans_sig <= 0.20, morans_sig_cat := 3]
  test.dt[morans_sig <= 0.10, morans_sig_cat := 4]
  test.dt[morans_sig <= 0.05, morans_sig_cat := 5]
  test.dt[morans_sig <= 0.01, morans_sig_cat := 6]
  test.dt[morans_sig <= 0.001, morans_sig_cat := 7]
  test.dt[, morans_sig_cat := as.character(morans_sig_cat)]
  test.dt[get(lisa_var) <= m.varofint & varlag <= m.varlag & morans_sig <= 0.05, morans_col_cat := 'Low-Low']
  test.dt[get(lisa_var) >= m.varofint & varlag >= m.varlag & morans_sig <= 0.05, morans_col_cat := 'High-High']
  test.dt[is.na(morans_col_cat), morans_col_cat := 'Negative']
  morans_color_map <- c("Low-Low"='blue',"High-High"='red',"Negative"='grey')
  
  map <- make_county_map(map_dt = test.dt[morans_sig <= sig, ],
                         map_var = 'morans_col_cat',
                         legend_title = 'Morans I',
                         high_is_good = FALSE,
                         map_title = paste0('Local Morans I Test (LISA)\nGlobal Morans I: ', round(mt.res$estimate[1], 2), ', p-value = ', round(mt.res$p.value, 4)),
                         map_sp = lisa_sp,
                         manual_colors = morans_color_map)
  
  pop_breaks = quantile(test.dt$total_pop, probs=c(.2,.3,.4,.5,.7,.9))
  i <- 6
  for(b in rev(pop_breaks)) {
    test.dt[total_pop <= b, pop_cat := as.character(i)]
    i <- i - 1
  }
  test.dt[is.na(pop_cat), pop_cat := '7']
  lisa_gg <- ggplot() + 
    geom_point(data = test.dt,
               aes(x = get(lisa_var),
                   y = varlag,
                   alpha = morans_sig_cat,
                   size = pop_cat,
                   fill = morans_col_cat),
               color = 'black',
               shape = 21,
               stroke = 0.1) +
    geom_hline(yintercept=m.varlag, color='black', size=1) +
    geom_vline(xintercept=m.varofint, color='black', size=1) +
    scale_alpha_manual(values = c("1"=.3,"2"=.3,"3"=.5,"4"=.5,"5"=1,"6"=1,"7"=1), guide = FALSE) + 
    scale_size_manual(values = c("1"=3,"2"=4,"3"=5,"4"=6,"5"=7,"6"=13,"7"=16), guide = FALSE) + 
    scale_fill_manual(values = morans_color_map, guide = FALSE) + 
    #scale_size_continuous(range = c(1,10), name = 'Population', breaks = quantile(test.dt$total_pop, probs=c(.2,.5,.8)), values = c(1,3,7)) + 
    #scale_size(range = c(3, 10), breaks = quantile(test.dt$total_pop, probs=c(.2,.5,.8))) + 
    theme_minimal() +
    xlab(lisa_var_name) + 
    ylab(paste0('Lagged ', lisa_var_name)) 
    #theme(legend.position="none")
  
  lisa_list <- list(map,lisa_gg,round(mt.res$estimate[1], 2),test.dt)
  return(lisa_list)
}

make_upper_lower <- function(d) {
  for(r in 1:length(d[, coef])) {
    draws <- exp(rnorm(1000, mean = d[r, coef], sd = d[r, se]))
    d[r, upper := quantile(draws, probs = .975)]
    d[r, lower := quantile(draws, probs = .025)]
  }
  return(d)
}
make_coef_table <- function(x) {
  m <- get(paste0('inla_model_',x))
  coefs <- data.table(model=paste0('Model ', x),
                      name=rownames(m$summary.fixed),
                      coef=m$summary.fixed$mean,
                      se=m$summary.fixed$sd)
  coefs <- make_upper_lower(coefs)
  return(coefs)
}


