#
#
#  "`-''-/").___..--''"`-._
# (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
# (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#   _ ..`--'_..-_/  /--'_.' ,'
# (il),-''  (li),'  ((!.-'
#
#
# Author: Guido Cervone (cervone@psu.edu) and Fangcao Xu (xfangcao@psu.edu)
#         Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#         Department of Geography and Institute for CyberScience
#         The Pennsylvania State University
#

library(reshape2)
library(doParallel)
library(foreach)
source("Z_global_variables.R")

### Read different types of files
MODTRAN.readfiles.digit <- function(filepaths){
  files <- lapply(filepaths, function(x) head(read.csv(x, skip=6, header = FALSE, stringsAsFactors = FALSE),-1))
  files <- lapply(files,setNames, c("freq", "path_trans", "path_emiss","path_thermal_scat","surface_emiss", "path_multi_scat",
                                    "sing_scat","grnd_rflt","drct_rflt","total_rad","irrad_ref","irrad_obs","-natlog_path_trans","direct_emiss",
                                    "ToA_irrad", "bbody_temp"))
  for(i in 1: length(files)){
    for(j in 1:ncol(files[[i]])){
      files[[i]][,j] <- as.numeric(files[[i]][,j])
    }
    files[[i]]$wavenumber <- files[[i]]$freq
    files[[i]] <- cbind(files[[i]][,17], files[[i]][,-c(1,17)])
    colnames(files[[i]])[1] <- "wavenumber"
    rownames(files[[i]]) <- 1: nrow(files[[i]])
  }
  return(files)
}

MODTRAN.readfiles.scan <- function(filepaths){
  # remove the first 6 and last rows
  files <- lapply(filepaths, function(x) head(read.csv(x, skip=6, header = FALSE, stringsAsFactors = FALSE),-1))
  files <- lapply(files, setNames, c("wavelength", "path_trans", "path_emiss","path_thermal_scat","surface_emiss", "path_multi_scat",
                                   "sing_scat","grnd_rflt","drct_rflt","total_rad","irrad_ref","irrad_obs","-natlog_path_trans","direct_emiss",
                                   "ToA_irrad", "bbody_temp"))
  files <- lapply(files, function(x) do.call(cbind, lapply(x, as.numeric)))
}


MODTRAN.combineGeo <- function(filelist, geoData){
  ### warning() here for duplicately appending one set of geometric information to each row of the associated result
  if(is.null(dim(geoData))){ #geoData is a vector
    # "only a vector of elevation angles is provided for the MODTRAN.combineGeo function"
    dataplot <- lapply(1:length(filelist), function(i) cbind(filelist[[i]], geoData[i], row.names = NULL))
    dataplot <- do.call(rbind, dataplot)
    colnames(dataplot)[ncol(dataplot)] <- "theta"
    
  }
  else{ #geoData is a matrix/data frame
    # "a matrix/dataframe of geo info is provided for the MODTRAN.combineGeo function"
    dataplot <- lapply(1:length(filelist), function(i) cbind(filelist[[i]], geoData[i,], row.names = NULL))
    dataplot <- do.call(rbind, dataplot)
    if(ncol(geoData)== 1){ # geoData[i,] is one value vector without original column name
      colnames(dataplot)[ncol(dataplot)] <- colnames(geoData)
    }
  }
  return(as.data.frame(dataplot))
}

MODTRAN.calculate.downwelling <- function(data.table){
  for(i in 1:nrow(data.table)){
    if(data.table$direct_emiss[i]==1){
      data.table$downwelling[i] <- NA
    }else{
      if(data.table$path_trans[i]==0){
        data.table$downwelling[i] <- 0
      }else{
        data.table$downwelling[i] <- (data.table$grnd_rflt[i]/data.table$path_trans[i])/(1-data.table$direct_emiss[i])
      }
    }
  }
  return(data.table)
}

# rewrite the table into a matrix with columns = elevation angles; rows axis= wavelength, value = radiance
# We assume that if wv is not passed, then table must have a column called wavelength
# We also assume that fi theta is not passed, then table must have a column called theta
MODTRAN.mat <- function(table, var, wv=NULL, theta=NULL, range=NULL, bh = FALSE) {
  # wavelength
  if (is.null(wv)) { wv <- table$wavelength} else {
    if (class(wv)=="character") {
      valid <- which(names(table) == wv )
      if( length(valid) > 0) { wv <- table[,valid] } else {return( NULL)}}}
  # theta
  if (is.null(theta)) {theta <- table$theta } else {
    if (class(theta)=="character") {
      valid <- which(names(table) == theta )
      if( length(valid) > 0) {theta <- table[,valid]} else {return( NULL)}}}
  # theta and range are all varying
  if (bh){
    if (is.null(range)) {range <- table$range } else {
      if (class(range)=="character") {
        valid <- which(names(table) == range )
        if( length(valid) > 0) {range <- table[,valid]} else {return( NULL)}}}
    anrange <- data.frame(unique(cbind(theta, range)))
    anrange <- do.call(paste, c(anrange, sep='_'))
  }else{
    theta <- rev(unique(theta))
  }
  wv <- unique(wv)
  # table here is ordered by wavelength from 7.5 to 11.9975 um and then ordered by the elvation angle
  id <- names(table) == var   # names() or colnames(); id: [FALSE, FALSE, TRUE...]
  if (any(id)) {
    mat = matrix(table[,id],nrow=length(wv))  # put the value into a matrix, byrow = FALSE
    rownames(mat) <- wv
    if(bh){
      colnames(mat) <- anrange
    }
    else{
      mat <- mat[,ncol(mat):1] # flip so the angles are sorted smallest to largest
      colnames(mat) <- theta
    }
    return(mat)
  } else {return( NULL)}
}


# change how many cores to use 
# elevation.agngles: elevation angle
# downNA whether allow NA to write into downwelling files or not
# bh: blue heron simulated data which not only has varying angle but also has varying range
writeMODTRAN2Matrix <- function(filesfolder, outputfolder, filename.pattern=c("Geometry", "Day","Time","Reflectivity"), folderpattern=folder.patterns, 
                                elevation.angles=eles, ncores=3, downNA = FALSE, writeinsame=FALSE, solar=FALSE, bh = FALSE){
  temp.scan <- list.files(filesfolder, pattern="*scan.csv",full.names=TRUE)
  temp.scan <- mixedsort(temp.scan)
  # whether to write the output into different folders or the same folder
  if(writeinsame == FALSE){
    outputfolder <- paste(outputfolder,folderpattern, sep="")
    names(outputfolder) <- c("trans","grnd","down","up1","up2", "surf", "solar1", "solar2","total")
    sapply(outputfolder, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))
  }else{
    outputfolder <- rep(paste(outputfolder,"/", sep=""), 9)
    names(outputfolder) <- c("trans","grnd","down","up1","up2", "surf", "solar1", "solar2", "total")
  }
  # make sure the geometry, day, and time info are provide
  geo.index <- which(filename.pattern=='Geometry')
  if(length(geo.index)==0 || length(which(filename.pattern=='Day'))==0 || length(which(filename.pattern=='Time'))==0) 
    stop("filename must contain Geometry Day and Time info")

  # separate different geos by the same day, time, reflectivity/temperature....
  vars <- return.index(temp.scan, range=1:length(filename.pattern), names = filename.pattern)
  vars.bygeo <-split(vars, f= lapply(filename.pattern[-geo.index], function(x) vars[,x])) 
  index <- data.frame(matrix(as.numeric(unlist(strsplit(names(vars.bygeo), "\\."))), ncol=length(filename.pattern[-geo.index]), byrow=T)) 
  names(index) <- filename.pattern[-geo.index]
  index$Time <- ifelse(index$Time-4 <0, index$Time+20, index$Time-4)  # change GMT time to local time
  # check whether some geo are missing in simulation, this function is using the zenith angles
  missing.set <- missing.data.simulation(vars.bygeo, if(is.null(dim(elevation.angles))) 90-elevation.angles else 90-elevation.angles[,1])
  if(length(missing.set)!=0){
    missing.index <- as.numeric(names(missing.set))
    for (i in missing.index){
      text <- paste("Missing elevation angles found in", do.call(paste, c(as.list(rep(" {} {}", length(filename.pattern[-geo.index]))), sep=",")), sep="")
      for( j in 1:length(filename.pattern[-geo.index])){
        text <- sub("\\{\\}",filename.pattern[-geo.index][j],  text)
        text <- sub("\\{\\}",index[i, j],  text)
      }
      print(text)}}
  # read the data into dataframe
  temp.scan.bygeo <-  lapply(vars.bygeo, function(x) temp.scan[as.numeric(rownames(x))]) # get original file paths for files in the same group
  table.scan.bygeo <- lapply(temp.scan.bygeo, function(x) tryCatch({MODTRAN.calculate.downwelling(MODTRAN.combineGeo(MODTRAN.readfiles.scan(x), elevation.angles))}, 
                                                                         error=function(e){stop("0 entry in filelist for some specific day, time ...")}))
  # do parallel processing
  registerDoParallel(cores=ncores)  
  getDoParWorkers()
  foreach(i=1:length(table.scan.bygeo)) %dopar%{
    table.scan <- table.scan.bygeo[[i]]
    idxname <- paste(index[i,], collapse ='_')
    # write downwelling with NA values into csv when downNA = T or no NA value
    if (downNA || !all(is.na(table.scan$downwelling))){
      mat.down <- MODTRAN.mat(table.scan,"downwelling", bh = bh)
      write.csv(mat.down, file = paste(outputfolder["down"],"Radiance_", idxname,"_down.csv",sep = ""))
    }
    if(solar){
      mat.solar1 <- MODTRAN.mat(table.scan,"sing_scat", bh = bh)
      mat.solar2 <- MODTRAN.mat(table.scan,"path_multi_scat", bh = bh)
      write.csv(mat.solar1, file = paste(outputfolder["solar1"],"Radiance_",idxname,"_solar1.csv",sep = ""))
      write.csv(mat.solar2, file = paste(outputfolder["solar2"],"Radiance_",idxname,"_solar2.csv",sep = ""))
    }
    mat.trans <- MODTRAN.mat(table.scan,"path_trans", bh = bh)
    mat.total_rad <- MODTRAN.mat(table.scan,"total_rad", bh = bh)
    mat.up1 <- MODTRAN.mat(table.scan,"path_emiss", bh = bh)
    mat.up2 <- MODTRAN.mat(table.scan,"path_thermal_scat", bh = bh)
    mat.grnd <- MODTRAN.mat(table.scan,"grnd_rflt", bh = bh)
    mat.surf <- MODTRAN.mat(table.scan,"surface_emiss", bh = bh)
    write.csv(mat.trans, file = paste(outputfolder["trans"],"Radiance_",idxname, "_trans.csv",sep = ""))
    write.csv(mat.total_rad, file = paste(outputfolder["total"],"Radiance_",idxname, "_total.csv",sep = ""))
    write.csv(mat.up1, file = paste(outputfolder["up1"],"Radiance_", idxname,"_up1.csv",sep = ""))
    write.csv(mat.up2, file = paste(outputfolder["up2"],"Radiance_", idxname, "_up2.csv",sep = ""))
    write.csv(mat.grnd, file = paste(outputfolder["grnd"],"Radiance_",idxname,"_grnd.csv",sep = ""))
    write.csv(mat.surf, file = paste(outputfolder["surf"],"Radiance_",idxname,"_surf.csv",sep = ""))
    return(TRUE)
  }
}


# Only write the MODTRAN simulated files for a specific azimuth, day, time or reflectivity
writeMODTRAN2Matrix.singleset <- function(filesfolder, outputfolder, folderpattern=folder.patterns, day=107, time=14,
                                          elevation.angles=eles, reflec=NULL, azimuth=NULL, downAllNA = FALSE, writeinsame=FALSE, solar=FALSE){
  temp.scan <- list.files(filesfolder, pattern="*scan.csv",full.names=TRUE)
  temp.scan <- mixedsort(temp.scan)
  # read the data into dataframe
  files.scan <- MODTRAN.readfiles.scan(temp.scan)
  table.scan.geo <- MODTRAN.combineGeo(files.scan, elevation.angles)  # one data frame for each set of day, time and reflectivity
  table.scan<- MODTRAN.calculate.downwelling(table.scan.geo)
  # whether to write the output into different folders or same folder
  if(writeinsame == FALSE){
    outputfolder <- paste(outputfolder,folderpattern, sep="")
    names(outputfolder) <- c("trans","grnd","down","up1","up2", "surf", "solar1", "solar2","total")
    sapply(outputfolder, function(x) if(!dir.exists(x)) dir.create(x))
  }else{
    outputfolder <- rep(paste(outputfolder,"/", sep=""), 9)
    names(outputfolder) <- c("trans","grnd","down","up1","up2", "surf", "solar1", "solar2", "total")
  }
  # whether has a specific azimuth or not
  if(is.null(azimuth)) azimuth<-"" else azimuth<-paste(azimuth, "_", sep="")
  # whether has a specific reflectivity or not
  if(is.null(reflec)) reflec<-"" else reflec<-paste(reflec, "_", sep="")
  
  # write downwelling with NA values into csv
  if(downAllNA && !all(is.na(table.scan$downwelling))){ 
    # if not all NA for downwelling radiance
    mat.down <- MODTRAN.mat(table.scan,"downwelling")
    write.csv(mat.down, file = paste(outputfolder["down"],"Radiance_",azimuth, day,"_",time,"_",reflec,"down.csv",sep = ""))
  }else{
    # if not any NA for downwelling radiance
    if(!any(is.na(table.scan$downwelling))){
      mat.down <- MODTRAN.mat(table.scan,"downwelling")
      write.csv(mat.down, file = paste(outputfolder["down"],"Radiance_",azimuth, day,"_",time,"_",reflec,"down.csv",sep = ""))
    }
  }
  
  if(solar){
    mat.solar1 <- MODTRAN.mat(table.scan,"sing_scat")
    mat.solar2 <- MODTRAN.mat(table.scan,"path_multi_scat")
    write.csv(mat.solar1, file = paste(outputfolder["solar1"],"Radiance_",azimuth, day,"_",time,"_",reflec,"solar1.csv",sep = ""))
    write.csv(mat.solar2, file = paste(outputfolder["solar2"],"Radiance_",azimuth, day,"_",time,"_",reflec,"solar2.csv",sep = ""))
  }
  mat.trans <- MODTRAN.mat(table.scan,"path_trans")
  mat.total_rad <- MODTRAN.mat(table.scan,"total_rad")
  mat.up1 <- MODTRAN.mat(table.scan,"path_emiss")
  mat.up2 <- MODTRAN.mat(table.scan,"path_thermal_scat")
  mat.grnd <- MODTRAN.mat(table.scan,"grnd_rflt")
  mat.surf <- MODTRAN.mat(table.scan,"surface_emiss")
  write.csv(mat.trans, file = paste(outputfolder["trans"],"Radiance_",azimuth, day,"_",time,"_",reflec, "trans.csv",sep = ""))
  write.csv(mat.total_rad, file = paste(outputfolder["total"],"Radiance_",azimuth, day,"_",time,"_",reflec,"total.csv",sep = ""))
  write.csv(mat.up1, file = paste(outputfolder["up1"],"Radiance_",azimuth, day,"_",time,"_",reflec, "up1.csv",sep = ""))
  write.csv(mat.up2, file = paste(outputfolder["up2"],"Radiance_",azimuth, day,"_",time,"_",reflec,"up2.csv",sep = ""))
  write.csv(mat.grnd, file = paste(outputfolder["grnd"],"Radiance_",azimuth, day,"_",time,"_",reflec,"grnd.csv",sep = ""))
  write.csv(mat.surf, file = paste(outputfolder["surf"],"Radiance_",azimuth, day,"_",time,"_",reflec,"surf.csv",sep = ""))
  return(TRUE)
}
