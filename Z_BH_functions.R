# A neat function that given a vector of values, will return colors codes.
library(plotrix)

# A neat function that given a vector of values, will return colors codes.
val2col <- function(x, min=NULL, max=NULL, na.rm=F, col) {
  if (is.null(min)) {min=min(x, na.rm=na.rm)}
  if (is.null(max)) {max=max(x, na.rm=na.rm)}
  
  x <- round( ( (x-min) / (max-min) ) * (length(col)-1)  + 1 )
  x[x<1] = 1
  x[x>length(col)] = length(col)
  return(col[x])
}

# decide the extent for plot
decideExtent <- function(turretinfo, scene.extent){
  plane.extent <- extent(min(turretinfo$longs.gmb), max(turretinfo$longs.gmb), min(turretinfo$lats.gmb),  max(turretinfo$lats.gmb))
  # scene_extent out of plane extent or not
  extent.logi <- (scene.extent[1] < plane.extent[1]) | (scene.extent[2] > plane.extent[2]) | 
    (scene.extent[3] < plane.extent[3]) | (scene.extent[4] > plane.extent[4])
  if(extent.logi){
    ext <- extent(min(scene.extent[1], plane.extent[1]), max(scene.extent[2], plane.extent[2]),
                  min(scene.extent[3], plane.extent[3]),  max(scene.extent[4], plane.extent[4]))
  }else{ext <- plane.extent}
  return(ext)
}



######################################################################################################################
######################################################################################################################
################################################# Read and plot the data #############################################
######################################################################################################################
######################################################################################################################
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# wrap into functions
# Find all pixels within the buffer of the target from SpectralGeo datatables
pixelsReader <- function(target.sp, files, sides, ratio=0.76, threshold=4e-6){
  search.r <- gBuffer(target.sp, width=threshold, byid=TRUE)
  pixels <- data.table()
  for(i in 1:length(files)){
    tmp <-  loadRData(files[i]); tmp$side <- sides[i]
    idx <- which(!is.na(over(SpatialPoints(tmp[,c(2,3)], proj4string=crs(search.r)), search.r)))
    if(length(idx) < 1) print(sprintf("Increase the threshold, no pixel found within the buffer for %s",  basename(files[i])))
    pixels <- rbind(pixels, tmp[idx,])
  }
  # assign which side the scene is at
  setcolorder(pixels, c(1, ncol(pixels), 2:(ncol(pixels)-1)))
  colnames(pixels)[10:265] <- paste("wv", 1:256, sep="")
  return(pixels)
}

remove.outliers <- function(dtable, usecols){
  stats.points <- Reduce('*', dtable[, ..usecols])
  outliers.id <- which(stats.points %in% boxplot.stats(stats.points)$out)
  return(dtable[-outliers.id, ])
}

# generate an image of selected scenes and one scene to extract pixels
plot_extractedPixels <- function(target.sp, specfiles, idscene, filename=NULL, 
                                 ratio=0.76, threshold=4e-6, textpos=4, textoffset=0.5){
  search.r <- gBuffer(target.sp, width=threshold, byid=TRUE)
  ext1 <- extent(search.r)
  # plot one scene
  im <-  loadRData(specfiles[idscene])
  # which target point is selected from all target points, could be more than 1, just use first one
  tarid <- which(!is.na(over(target.sp, footprints.spdf[unique(im$sid),])[,1]))[1]
  # points within search buffer of selected one target
  idx <- which(!is.na(over(SpatialPoints(im[,c(2,3)], proj4string=crs(search.r)), search.r[tarid])))
  
  # start to plot
  if(!is.null(filename)) jpeg(filename, width=3000, height=1200, res=300)
  # 1. plot all extracted scenes
  par(mfrow=c(1,2), mar=c(3,3,1,0.5), mgp = c(2, 0.7, 0))
  ext <- extent(gBuffer(target.sp, width=.005))
  plot(ext,col="NA", xlab="Longitude", ylab="Latitude", xaxt='n', yaxt='n')
  xlab <- pretty(range(ext[1], ext[2]), 4)
  ylab <- pretty(range(ext[3], ext[4]), 4)
  axis(1,at=xlab,labels=xlab)
  axis(2,at=ylab,labels=ylab)
  plot(footprints.spdf,add=T,col="grey")
  plot(footprints.spdf[sids,],col="blue",add=T)
  plot(target.sp, add=T,col="red",lwd=2)
  # 2. show extracted pixels in one image given the target
  ext.center <- extent(c(target[tarid,]-1e-5, target[tarid,]+2e-5)[c(1,3,2,4)])
  plot(ext.center, col="white", xlab="Longitude", ylab="Latitude", xaxt='n', yaxt='n')
  xlab <- pretty(range(ext.center[1], ext.center[2]), 4)
  ylab <- pretty(range(ext.center[3], ext.center[4]), 4)
  axis(1,at=xlab,labels=xlab)
  axis(2,at=ylab,labels=ylab)
  plot(SpatialPointsDataFrame(im[,c(2,3)], im, proj4string = crs.ll), col="grey40", add=T)
  plot(target.sp[tarid], col="red", lwd=2, add=T)
  plot(search.r, lwd=1, add=T)
  # This is due to 1 degree of longitude is shorter than 1 degree of latitude
  points(im[idx, c(2,3)], col=rainbow(length(idx)), cex=1.5)
  mapply(function(i, col) lines(c(target[tarid,1], im[i, 2]),c(target[tarid,2], im[i, 3]), col=col), idx, rainbow(length(idx)))
  dist_m <- format(distm(im[idx, c(2,3)], target[tarid,], fun=distGeo), digits=3)
  text(unlist(im[idx, c(2,3)][,1]), unlist(im[idx, c(2,3)][,2]),  paste(dist_m, "m"), 
      adj=textoffset, cex=0.8)
  
  if(!is.null(filename)) dev.off()
}


######################################################################################################################
######################################################################################################################
############################################## Spectral Processing and Plot ##########################################
######################################################################################################################
######################################################################################################################
scenepoints2shapefile <- function(scenes, folderid, res.dir, crs.ll=CRS("+proj=longlat")){
  tmp <- NULL; tmp.center <- NULL;
  for(i in 1:length(scenes)){
    scene <- scenes[[i]]
    # all ground points and center point
    center.grnd <- scene[round(nrow(scene)/2), c("sid", "sid_2","folderid", "longs", "lats", "alt")]
    center.grnd$type <- "ground"
    center.gmb <- scene[round(nrow(scene)/2), c("sid", "sid_2", "folderid","longs.gmb", "lats.gmb", "alt.gmb")]
    colnames(center.gmb) <- c("sid", "sid_2", "folderid", "longs", "lats", "alt")
    center.gmb$type <- "plane"
    tmp.center <- rbind(tmp.center, center.grnd, center.gmb)
    # all gmb points and center point
    grnd <- scene[, c("sid", "sid_2","folderid", "longs", "lats", "alt")]
    grnd$type <- "ground"
    gmb <- scene[, c("sid", "sid_2", "folderid","longs.gmb", "lats.gmb", "alt.gmb")]
    colnames(gmb) <- c("sid", "sid_2", "folderid", "longs", "lats", "alt")
    gmb$type <- "plane"
    tmp <- rbind(tmp, grnd, gmb)
  }
  rownames(tmp.center) <- 1:nrow(tmp.center)
  center.spdf <- SpatialPointsDataFrame(tmp.center[, c("longs", "lats")], data=tmp.center, proj4string = crs.ll)
  points.spdf <- SpatialPointsDataFrame(tmp[, c("longs", "lats")], data=tmp, proj4string = crs.ll)
  writeOGR(center.spdf, res.dir, paste(folderid, "_center", sep=""), driver="ESRI Shapefile", overwrite_layer = T)
  writeOGR(points.spdf, res.dir, paste(folderid, "_points", sep=""), driver="ESRI Shapefile", overwrite_layer = T)
}

# Read the geolocation ifnormation
# dimensions: 367, 258, 94686, 3  (nrow, ncol, ncell, nlayers), here the nrow, ncol, ncell may vary for different images
# nlayers: represent for lon, lat, dem or spectral wavelength
specinfo_read <- function(file, spec=NULL){
  ras <- brick(file)      # nrow, ncol, ncell, 256 (256 bands)/3 (lat, lon, dem)
  if(!is.null(spec)){
    ras <- ras[[spec]]
  }
  df <- data.table(matrix(as.vector(ras), ncol=nlayers(ras)))
  if(nlayers(ras)==3) colnames(df)<- c("lon", "lat", "dem")
  else colnames(df) <- names(ras)
  return(df)
}

# read geo and hyperspectral info of one scene
dem_spec_scene<- function(scene, dir.name, spec=NULL){
  hdr.fname <- id_scene[unlist(scene[1,c(1,2)]),]$filename
  geo.fname <- gsub("hdr$","geo",hdr.fname)
  img.fname <- gsub("hdr$","img",hdr.fname)
  # location and dem file
  geo.full.fname <- paste(dir.name,"/BH_",substring(geo.fname,1,15),"/nav/",geo.fname,sep="")
  # hyperspectral image of 256 bands
  img.full.fname <- paste(dir.name,"/BH_",substring(geo.fname,1,15),"/hsi_lwir_",c(1,2),"/cal/",img.fname,sep="")
  # whether geo and spectral files exist
  idx <- which(!file.exists(geo.full.fname) | !file.exists(img.full.fname))
  if (length(idx) > 0) stop(paste("geo and spectral files not exist for", hdr.fname[idx]))
  if (!is.null(spec)){
    # spectral file y read specific wavelengths
    f <- function(x, y, sid){return(cbind(sid = sid, specinfo_read(x), specinfo_read(y, spec)))}
  }
  else{f <- function(x, y, sid){return(cbind(sid = sid, specinfo_read(x), specinfo_read(y)))}}
  # list of two matrix for each side scene
  dem_spec_info <- mapply(f, geo.full.fname, img.full.fname, scene[1,c(1,2)], SIMPLIFY=F)
  return(dem_spec_info)
}

dem_spec_scene2center <- function(dem_spec_info, plane.center){
  # plane.center has three values: c("longs.gmb", "lats.gmb", "alt.gmb")
  dmat <- distm(dem_spec_info[,c("lon", "lat")], plane.center[c(1,2)], fun=distGeo)
  dem_spec_info$center.dist <- dmat
  dem_spec_info$center.alt <- plane.center[3]- dem_spec_info$dem
  dem_spec_info$center.range <- sqrt(dem_spec_info$center.dist^2 + dem_spec_info$center.alt^2)
  dem_spec_info$center.angle <- 180*asin(dem_spec_info$center.alt/dem_spec_info$center.range)/pi
  setcolorder(dem_spec_info, c(1:4, (ncol(dem_spec_info)-3):ncol(dem_spec_info), 5:(ncol(dem_spec_info)-4)))
  return(dem_spec_info)
}

##########################
plot_scenes <- function(scenes, ext, cex.point = 0.1, cex.center = 0.5, foot.col="grey", gmb.col="grey", gmb.mean.col="blue",
                        target.col="grey", target.mean.col="red", plottime=NULL, sleeptime = 1, centerline = F){
  plot(ext, xlab="Longitude", ylab="Latitude", col="white", xaxt="n",yaxt="n")
  # beautify the x y lables
  x1 = floor(ext[1]*100)/100; x2 = ceiling(ext[2]*100)/100
  y1 = floor(ext[3]*100)/100; y2 = ceiling(ext[4]*100)/100
  xlab = pretty(range(x1, x2), 4)
  ylab = pretty(range(y1, y2), 4)
  axis(1,at=xlab,labels=xlab)
  axis(2,at=ylab,labels=ylab)
  for(i in 1:length(scenes)){
    scene <- scenes[[i]]
    sids <- c(scene$sid[1],  scene$sid_2[1])
    plot(footprints.spdf[sids,], col=foot.col, add=T)
    points(scene[, c("longs", "lats")], col= target.col, cex = cex.point)  # targets on the ground
    points(scene[, c("longs.gmb", "lats.gmb")], col= gmb.col, cex = cex.point)  # plane
    points(scene[round(nrow(scene)/2), c("longs", "lats")], col= target.mean.col, cex = cex.center, pch=19)
    points(scene[round(nrow(scene)/2), c("longs.gmb", "lats.gmb")], col= gmb.mean.col, cex = cex.center, pch=19)  # plane
    if(centerline) lines(scene[round(nrow(scene)/2), c("longs", "longs.gmb")], scene[round(nrow(scene)/2), c("lats", "lats.gmb")], col="gray")
    if(!is.null(plottime) & i %in% plottime[,1]){
      idx <- which(plottime[,1]==i)
      text(scene[round(nrow(scene)/2), "longs.gmb"], scene[round(nrow(scene)/2), "lats.gmb"],
           as.ITime(scene[round(nrow(scene)/2), "time"]), 
           pos=plottime[idx,2], offset=plottime[idx,3], srt=plottime[idx,4])
    }
    Sys.sleep(sleeptime)
  }
}

# plot spectral
plot_dem_spec_twosides<-function(dem_spec, id_var, col = rev(brewer.pal(11,"Spectral")), minmax = NULL, 
                                 ext=NULL, axi_lable=T, legend=F, add=F){
  side1.spdf <- SpatialPointsDataFrame(dem_spec[[1]][,c("lon", "lat")], data=dem_spec[[1]], proj4string = crs.ll)
  side2.spdf <- SpatialPointsDataFrame(dem_spec[[2]][,c("lon", "lat")], data=dem_spec[[2]], proj4string = crs.ll)
  if (!add){
    if(is.null(ext)) ext <- extent(rbind(side1.spdf@coords,side2.spdf@coords))
    if(!axi_lable) plot(ext, xlab="", ylab="", col="white", xaxt="n",yaxt="n")
    else{
      plot(ext, xlab="Longitude", ylab="Latitude", col="white", xaxt="n",yaxt="n")
      xlab = pretty(range(ext[1], ext[2]), 4)
      ylab = pretty(range(ext[3], ext[4]), 3)
      axis(1,at=xlab,labels=xlab)
      axis(2,at=ylab,labels=ylab)
    }
    if(legend){
      lgvalues <- seq(minmax[1],minmax[2],10)
      color.legend(ext[2]-0.0022,ext[3]+(ext[4]-ext[3])/4,ext[2]-0.0002,ext[4]-(ext[4]-ext[3])/4, 
                   c(lgvalues[1], lgvalues[length(lgvalues)%/%2], lgvalues[length(lgvalues)]), 
                   colorRampPalette(col)(length(lgvalues)), align="lt",gradient="y", cex=0.7)
      mtext(bquote(10^-6~"Wcm"^-2*"sr"^-1*mu*"m"^-1),side=4, cex=0.7, line=-0.9, srt=360)
    }
  }
  if(!is.null(minmax)){
    plot(side1.spdf, col=val2col(side1.spdf@data[, id_var], min=minmax[1], max=minmax[2], col=col), pch=20,cex=0.5, add=T)
    plot(side2.spdf, col=val2col(side2.spdf@data[, id_var], min=minmax[1], max=minmax[2], col=col), pch=20,cex=0.5, add=T)
  }else{
    plot(side1.spdf, col=val2col(side1.spdf@data[, id_var], col=col), pch=20,cex=0.5, add=T)
    plot(side2.spdf, col=val2col(side2.spdf@data[, id_var], col=col), pch=20,cex=0.5, add=T)
  }
}
