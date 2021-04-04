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

### folder paths
SimulationGeometry.dir <- "~/geolab_storage_V3/data/DARPA/SimulationGeometry"


### file and folder pattern
file.patterns= c("/TransmissionCSV/Radiance_{day}_{hh}_{reflect}_trans.csv",
                 "/GroundReflectedCSV/Radiance_{day}_{hh}_{reflect}_grnd.csv",
                 "/DownwellingCSV/Radiance_{day}_{hh}_{reflect}_down.csv",
                 "/PathThermalEmissionCSV/Radiance_{day}_{hh}_{reflect}_up1.csv",
                 "/PathThermalScatteringCSV/Radiance_{day}_{hh}_{reflect}_up2.csv",
                 "/SurfaceEmissionCSV/Radiance_{day}_{hh}_{reflect}_surf.csv",
                 "/SolarSingleScatteringCSV/Radiance_{day}_{hh}_{reflect}_solar1.csv",
                 "/SolarMultiScatteringCSV/Radiance_{day}_{hh}_{reflect}_solar2.csv",
                 "/TotalRadianceCSV/Radiance_{day}_{hh}_{reflect}_total.csv")
folder.patterns= c("/TransmissionCSV/","/GroundReflectedCSV/", "/DownwellingCSV/", "/PathThermalEmissionCSV/", "/PathThermalScatteringCSV/", 
                   "/SurfaceEmissionCSV/", "/SolarSingleScatteringCSV/","/SolarMultiScatteringCSV/", "/TotalRadianceCSV/")
# up1: path thermal emission; up2: path scattering
names(file.patterns) <- c("trans","grnd","down","up1","up2", "surf","solar1","solar2","total") 
names(folder.patterns) <- c("trans","grnd","down","up1","up2", "surf","solar1","solar2","total")

### load geometric csv files
# sensor elevation
fname <- file.path(SimulationGeometry.dir, "sensorElevation.csv")
eles <- read.csv(fname)[,-1]
eles$theta <- 90-eles$theta
# sensor location
fname <- file.path(SimulationGeometry.dir, "sensorLocation.csv")
locs <- read.csv(fname)[,-1]
locs <- locs[c(1,74:nrow(locs)),]
locs$theta <- 90-locs$theta
theta.list <- unique(locs$theta)
phi.list <- unique(locs$phi)

### boxplot color for different elevation angles
boxplot.color <-  c("#00A600FF","#1DB000FF","#3EBB00FF","#63C600FF","#8BD000","#C4E700","#FFFF00","#FFEB00","#FFD800","#FFC500","#FFB200","#FF9F00","#FF8C00")

### geometric variables
wv.default <- seq(7.5, 12, 0.0175)[-c(1,2)]
wv.BH.default<- seq(7.46, 13.5300, 0.0238)
# selected angles for the for radiance matrix 
an.full <- seq(30, 90, 5)
an.default <- c(30, 35, 40, 50, 60, 70, 80, 90)
an.BH.default<- c(30, 35, 40, 45, 50, 55, 60)

GMTtime <- c(2,6,10,14,18,22)
Localtime <- sapply(GMTtime, function(x) ifelse(x-4 <0, x+20, x-4))
reflectivity.default<-c(0,5,10,30,50,80,100)

