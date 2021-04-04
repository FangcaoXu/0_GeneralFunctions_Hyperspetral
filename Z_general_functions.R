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

library(gtools)
#source("Z_global_variables.R")

create.results.dir <- function(base.dir=".", prefix="results_") {
  # Let's get the time to add to the output
  ctime <- format(Sys.time(),"%Y%m%d_%H%M%S")
  
  # Does dir already exist or not?
  if (!file.exists(base.dir) ) {
    dir.create(base.dir)
  } else if (!file.info(base.dir)$isdir) {
    file.rename(base.dir,paste(base.dir,ctime,sep=""))
  }
  
  # Create a directory for our stuff
  res.dir <- paste(base.dir,"/",prefix, ctime,"/",sep="")
  dir.create(res.dir)
  
  return(res.dir)
}

# filerename from LWIRRadiane Radiance
LWIRRadiance2Radiance <- function(files, pattern="LWIRRadiance", replacement="Radiance"){
  rename <- function(eachfile){
    file.rename(eachfile, sub(pattern, replacement, eachfile))
  }
  sapply(files, rename)
}

# files <- list.files('~/disk10TB/DARPA/MatrixCSV', '.csv', full.names = TRUE, recursive=TRUE)
# LWIRRadiance2Radiance(files)

# Normalize a value between 1..nmax
normalize1nmax <- function(x, na.rm=T, nmax=2) {
  norm <- ( (x - min(x, na.rm=na.rm)) / (max(x, na.rm=na.rm)-min(x, na.rm=na.rm)) ) 
  norm <- norm * (nmax-1) +1
  return (norm)
}

# ### Find outliers of boxplots
# wrongsimulatedData <- function(boxdata){
#   tmp <-boxdata[which(boxdata$path_trans %in% boxplot(boxdata$path_trans~boxdata$wavelength)$out),]
#   geometry.obs <- unique(tmp[,c(17:20)])
# }

# rename tables
rowcol.renames <-function(table, r.names, c.names){
  rownames(table) <- r.names
  colnames(table)<- c.names
  return(table)
}


# get the maximum value for a specific variable across all tables
minmax.crosstables <- function(tables, var=NULL){
  max.values <- vector('numeric')
  min.values <- vector('numeric')
  for(table in tables){
    if (is.null(var)){
      min.values <- c(min.values,min(table))
      max.values <- c(max.values,max(table))
    }else{
      id <- names(table) == var   # names() or colnames(); id: [FALSE, FALSE, TRUE...]
      if (any(id)) {
        min.values <- c(min.values,min(table[,id]))
        max.values <- c(max.values,max(table[,id]))
      }else{
        return(NULL)
      }
    }
  }
  min.value <- min(min.values)
  max.value <- max(max.values)
  return(c(min.value,max.value))
}

# calculate the residual of predicted matrix and observed matrix 
# Error is the difference between the observed value in a sample and the true value in the population. 
# Residual is the difference between the observed value and the predicted value.
diff.matrices <- function(predictmatrix, originalmatrix, square = F){
  if(square){
    return((predictmatrix-originalmatrix)^2 )
  }else{
    return(predictmatrix-originalmatrix)
  }
}

find.originalfiles <- function(predictfiles, originalfilefolder){
  # find original files
  # (.*?) ?: match at most 1 times which means extract first match
  components <- gsub(".*_([[:alnum:]]+)\\.csv","\\1",predictfiles)
  predictfiles.basename <-basename(predictfiles)
  originalfiles <- paste(originalfilefolder,folder.patterns[components], predictfiles.basename,sep = "")
}

find.originaltotalInput <- function(predictfiles, originalfilefolder){
  # find original files
  # (.*?) ?: match at most 1 times which means extract first match
  predictfiles.basename <-basename(predictfiles)
  components <- gsub(".*_([[:alnum:]]+)\\.csv","\\1",predictfiles.basename)
  predictfiles.basename <- unname(mapply(function(x,y) gsub(x,"total", y), components, predictfiles.basename))
  originalfiles <- paste(originalfilefolder,folder.patterns["total"], predictfiles.basename,sep = "")
}

# find the index of given data points nearest to the standard data point
findnearest <- function(given, standard){
  index <- which.min(sapply(given, function(y) sum((standard-y)^2)))
  return(index) # return index in given closest to standard
}


pretty.axislabels <- function(input, pretty){
  # set the cutting threshold
  input.pretty <- findnearest(input, pretty)
  at <- input[unique(input.pretty)]
  # get the true unique values (which has never been duplicated)  
  standard.pretty <- which(!(duplicated(input.pretty)|duplicated(input.pretty, fromLast=TRUE)))

  # check whether two are within threshold
  if(length(standard.pretty)!=length(input.pretty)){
    if(standard.pretty[1] >1){
      start.index <- standard.pretty[1]-1
    }else{
      start.index <- standard.pretty[1]
    }
    if(standard.pretty[length(standard.pretty)] < length(pretty)){
      stop.index <- standard.pretty[length(standard.pretty)] +1
    }
    else{
      stop.index <- standard.pretty[length(standard.pretty)]
    }
    labels <- pretty[start.index : stop.index]
  }else{
    labels <- pretty[standard.pretty]
  }
  return(list(at, labels))
}



# read predicted matrix file, which has 8 columns
read.predictedmatrix <- function(predictedfiles, wv=wv.default, an=an.default){
  # name the file.read list by each file name
  predictedfiles.name <- sub("\\.csv$", "",basename(predictedfiles))
  # read predicted files
  predictedfiles.read <- lapply(predictedfiles, function(x) read.csv(x, header = FALSE, stringsAsFactors = FALSE))
  predictedfiles.read <- lapply(predictedfiles.read, function(x) rowcol.renames(x, wv, an))
  names(predictedfiles.read) <- predictedfiles.name
  # read original files, skip first three rows and first column
  return(predictedfiles.read)
}
# read original matrix file, which has 13 columns
read.originalmatrix <- function(originalfiles, wv=wv.default, an=an.default, digits=4){
  wv <- round(wv, digits)
  originalfiles.name <- sub("\\.csv$", "",basename(originalfiles))
  originalfiles.read <- lapply(originalfiles, function(x) try(read.csv(x, row.names = 1, header = TRUE,                                                                      check.names=FALSE, stringsAsFactors = FALSE)))
  # only read rows and columns for given wavelength and angles
  wv.read <- as.numeric(rownames(originalfiles.read[[1]]))
  an.read <- as.numeric(colnames(originalfiles.read[[1]]))
  row.index <-  which(wv.read %in% wv)
  col.index <- which(an.read %in% an)
  originalfiles.read <- lapply(originalfiles.read, function(x) x[row.index, col.index])
  # name each file with their file basename
  names(originalfiles.read) <- originalfiles.name
  return(originalfiles.read)
}

# read predicted and original matrix files, which respectively have 8 columns and 13 columns
read.predorigmatrix <- function(predictfiles, originalfiles, wv=wv.default, an=an.default){
  # read predicted files
  predictfiles.read <- lapply(predictfiles, function(x) read.csv(x, header = FALSE, stringsAsFactors = FALSE))
  predictfiles.read <- lapply(predictfiles.read, function(x) rowcol.renames(x, wv, an))
  # name the file.read list by each file name
  predictfiles.name <- sub("\\.csv$", "",basename(predictfiles))
  names(predictfiles.read) <- predictfiles.name
  # read original files,
  originalfiles.read <- read.originalmatrix(originalfiles, wv, an)
  return(list(predictfiles.read, originalfiles.read))
}

# return the index of day, time and reflectivity based on the file based name
return.index <- function(filelists, range=1:4, names = c("Geometry", "Day","Time","Reflectivity")){
  f <- function(x){return(unlist(strsplit(x,"_"))[range])}
  ret <- lapply(basename(filelists), f)
  vars <- data.frame(matrix(unlist(ret),ncol=length(range),byrow=T), stringsAsFactors = FALSE)
  names(vars) <- names
  return(vars)
}


###### Check whether some date, time, reflectivity is missing in the simulated data
missing.data.simulation <- function(vars.byangles, an){ # Both angles here should be zeniths
  # "Input variable in missing.data.simulation function is a vector of zenith angles"
  angles.simulated <- lapply(vars.byangles, function(x) as.numeric(gsub(".*Radiance([0-9]+)", "\\1",x$Geometry)))
  missing.index <- which(sapply(angles.simulated, function(x) length(x))!=length(an))
  if(length(missing.index!=0)){
    missing.angles <- lapply(angles.simulated[missing.index], function(x) an[which(is.na(match(an, x)))]) 
    names(missing.angles) <- missing.index
  }else {
    missing.angles <- list()
  }
  return(missing.angles)
}

# filter out the matrix for the same variable and calculate their residuals
residual.prep <- function(predorigfiles.read, var="day", fix.var= c("time"=18,"reflectivity"=10,"angle"=60), wv=wv.default){
  predictfiles.read <- predorigfiles.read[[1]]
  originalfiles.read <- predorigfiles.read[[2]]
  filesbasename <- names(predictfiles.read)
  
  days <- as.numeric(gsub("Radiance_(.*?)_.*","\\1",filesbasename))
  times <- as.numeric(gsub(".*Radiance+_.*_(.*?)_.*","\\1",filesbasename))
  reflectivities <- as.numeric(gsub(".*_(.*?)_[[:alnum:]]+$","\\1",filesbasename))
  
  if(var=="day"){
    time.index <- which(times==fix.var["time"])
    reflec.index <- which(reflectivities==fix.var["reflectivity"])
    file.index <- intersect(time.index,reflec.index)
    predictfiles.read <- predictfiles.read[file.index]
    originalfiles.read<- originalfiles.read[file.index]
    
    col.index <- which(an==fix.var["angle"])
    # "[" whose first argument is the object being subsetted. Subsequent arguments are the index to that subset
    predictfiles.read <- lapply(predictfiles.read, function(x) x[,col.index,drop=FALSE]) 
    originalfiles.read<- lapply(originalfiles.read, function(x) x[,col.index,drop=FALSE]) 
  }
  
  else if(var=="time"){
    day.index <- which(days==fix.var["day"])
    reflec.index <- which(reflectivities==fix.var["reflectivity"])
    file.index <- intersect(day.index,reflec.index)
    predictfiles.read <- predictfiles.read[file.index]
    originalfiles.read<- originalfiles.read[file.index]
    
    col.index <- which(an==fix.var["angle"])
    predictfiles.read <- lapply(predictfiles.read, function(x) x[,col.index,drop=FALSE]) 
    originalfiles.read<- lapply(originalfiles.read, function(x) x[,col.index,drop=FALSE]) 
    
  }
  else if(var=="reflectivity"){
    day.index <- which(days==fix.var["day"])
    time.index <- which(times==fix.var["time"])
    file.index <- intersect(day.index,time.index)
    predictfiles.read <- predictfiles.read[file.index]
    originalfiles.read<- originalfiles.read[file.index]
    
    col.index <- which(an==fix.var["angle"])
    predictfiles.read <- lapply(predictfiles.read, function(x) x[,col.index,drop=FALSE])# don't convert one column to vector
    originalfiles.read<- lapply(originalfiles.read, function(x) x[,col.index,drop=FALSE]) 
    
  }
  else if(var=="angle"){
    day.index <- which(days==fix.var["day"])
    reflec.index <- which(reflectivities==fix.var["reflectivity"])
    time.index <- which(times==fix.var["time"])
    file.index <- Reduce(intersect, list(day.index,reflec.index,time.index))
    predictfiles.read <- predictfiles.read[file.index]
    originalfiles.read<- originalfiles.read[file.index]
    
  }else{
    stop('the dependent variable must be one from the day, time, reflectivity, and angle and the remaining ones should be fixed to a specfic value')
  }
  
  residuals <- data.frame(row.names = wv)
  for (i in 1:length(predictfiles.read)){ # 
    filebasename <- names(predictfiles.read)
    residuals <-cbind(residuals, diff.matrices(predictfiles.read[[i]],originalfiles.read[[i]])) 
  }
  return(residuals)
}

# prepare an array to store all matrix together for statistical calculation
statistic.prep <- function(input, check.residual=F){
  # marray: (256, 8, the number of files)
  if(check.residual){
    # input: predorigfiles
    matrixsize <- dim(input[[1]][[1]])
    filelengths <- length(input[[1]])
    marray <- array(NA, dim=c(matrixsize,filelengths))
    for (i in 1:(dim(marray)[3])){
      temp.pred <- as.matrix(input[[1]][[i]])
      temp.orig <- as.matrix(input[[2]][[i]])
      temp <- as.matrix(diff.matrices(temp.pred, temp.orig))
      marray[,,i] <- temp
    }
  }else{
    matrixsize <- dim(input[[1]])
    filelengths <- length(input)
    marray <- array(NA, dim=c(matrixsize,filelengths))
    for (i in 1:(dim(marray)[3])){
      temp <- as.matrix(input[[i]])
      marray[,,i] <- temp
    }
  }
  return(marray)
}

# calculate blackbody radiance
blackbody_radiance_calculation <- function(wv, temp=320, irradiance=F){
  print("Input wavelength unit: um \n temperature unit: K")
  h <- 6.6256*10^-34  # Planck constant, j*sec, W*sec^2  # watts, W: j/sec
  c <- 2.9979*10^8    # speed of light, m/sec #  # the unit of lambda and speed of light, 1 um ~ 10^-6 m 
  k <- 1.38*10^-23    # Bolzmann constant, j/K
  blackbody.emitted.radiance <- c()
  
  # https://www.spectralcalc.com/blackbody/CalculatingBlackbodyRadianceV2.pdf
  # https://www.spectralcalc.com/blackbody_calculator/blackbody.php
  if(irradiance){
    for(lambda in wv){
      blackbody.emitted.radiance <- c(blackbody.emitted.radiance, 2*pi*10^24*h*(c^2)*(lambda^-5)/(exp(10^6*h*c/(lambda*k*temp))-1))  # [W m^-2 um^-1]
    }
    print("Output Irradiance Unit: [W m^-2 um^-1] in blackbody_radiance_calculation function" )
  }
  else{
    for(lambda in wv){
      blackbody.emitted.radiance <- c(blackbody.emitted.radiance, 2*10^24*h*(c^2)*(lambda^-5)/(exp(10^6*h*c/(lambda*k*temp))-1))  # [W m^-2 sr^-1 um^-1]
    }
    print("Output Radiance Unit: [W m^-2 sr^-1 um^-1] in blackbody_radiance_calculation function")
  }
  return(blackbody.emitted.radiance)
}


# retrieve emissivity from radiative transfer equation
# modtran radiance unit: W cm-2 sr-1 um-1
emissivity.retrieval <- function(total, up1, up2, down, trans, temp=320, wv=wv.default, an=an.default){
  # length(which(trans<=10^-4)) # 14 wrong values
  surf.leaving <- (total - up1 - up2)/trans
  # calculate blackbody emittied radiance
  blackbody.emitted.radiance <- blackbody_radiance_calculation(wv, temp)
  blackbody.emitted.radiance <- blackbody.emitted.radiance*10^-4 # [W cm^-2 sr^-1 um^-1]
  # assume the blackbody radiance  won't change with different viewing angles
  blackbody.emitted.radiance <- matrix(rep(blackbody.emitted.radiance, 8), ncol=8)
  rownames(blackbody.emitted.radiance) <- wv
  colnames(blackbody.emitted.radiance) <- an
  # down*reflectivity + blackbody.emitted.radiance*(1-reflectivity) = surface leaving radiance
  reflectivity <- (surf.leaving - blackbody.emitted.radiance)/(down-blackbody.emitted.radiance)
}
