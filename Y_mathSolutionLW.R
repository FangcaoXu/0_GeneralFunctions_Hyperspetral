# persp is the base graphics function for creating wireframe surface plots
# The lattice graphics package has a function wireframe
library(colorRamps)
library(dplyr)
library(grid)
library(gridExtra)
library(reshape2)
library(stringr)
library(tidyverse)
library(viridis)
source("Z_general_functions.R")
### total assumption
### 1: lambertian surface; 2: clear sky; 3: gray-body object whose emissivity is constant for different lambdas; 
# 4: atmosphere temperature is constant over different layers

getwd()

# Section 1: Initiate parameters ----------------------------------------------------
h <- 6.6256 *10^-34  # Planck constant, j*sec
c <- 2.9979*10^8     # speed of light, m/sec
k <- 1.38*10^-23     # Bolzmann gas constant, j/K
s <- 5.67*10^-8      # Stefan–Boltzmann constant, W m-2 K-4
emiss.ls <- c(0.7, 0.8, 0.9, 1)       # emissivity
temp.ls <- c(0, 100,200,300,400)      # temperature, degree Kelvin
lambda.ls <- c(7,8,9,10,11,12,13)       # wavelength, um
theta.ls <-  seq(30,90,10)  # elevation angle, degree
range.ls <- seq(500, 10000, by=500) # range, meter
deg2rad <- pi/180


# Section 2: Profile atomosphere -----------------------------------------------------
# See Figure 3.17
# C.a = absorption cross section; m: number density of molecules; C.a cannot be measured directly
# layer.depth: defined atmosphere layer depth, [m]
# c.ppm: gas concentration parts per million [ppm] for each layer, 1 ppm = 1mg/L. For example, CO.2 (ppm) is around 400. 
# beta.cz: absorption per part per million per meter [ppm-1 m-1] for each layer
# Build the atmosphere profile from the surface to the outer space
c.ppm = rep(700, 100)
beta.cz = rep(0.1*10^-6, 100)
              
# Section 3: Radiance equation for LW -------------------------------------------------------------
# atmos.temp = set atmosphere temperature as constant, degree Kelvin
cal.transmission <- function(range, theta, layer.depth){
  nlayers<- range*sin(theta * deg2rad)/layer.depth
  delta.lambda <- 0
  if (floor(nlayers)>0){
    for(i in 1:floor(nlayers)){
      delta.lambda <- delta.lambda + c.ppm[i]*beta.cz[i] *layer.depth  # [m-1]
    }
  }
  delta.lambda <- delta.lambda + c.ppm[ceiling(nlayers)]*beta.cz[ceiling(nlayers)]*(nlayers-floor(nlayers))*layer.depth
  tau.lambda <- exp(-delta.lambda/sin(theta * deg2rad))
  return(tau.lambda)
}


# calculate atmosphere thermal emission radiance reaching at sensor from elevation theta
# [W m^-2 sr^-1 um^-1]
cal.atmosphere.emiss.oneangle <- function (lambda, theta, atmos.temp=200, layer.depth=500){
  BM.temp.lambda.atm <-  blackbody_radiance_calculation(lambda, atmos.temp)
  e.atm.lambda.layers <- NULL          # emissivity of the ith layer
  L.atm.lambda.layers<- NULL           # irradiance of the ith layer
  tau1.lambda.layers <- NULL
  tau1.lambda <- NULL
  for (i in 1: ceiling(2000/layer.depth)){
    delta.lambda.i <- c.ppm[i]*beta.cz[i] *layer.depth   # optical depth of the ith layer
    tau1.lambda.layers[i] <- exp(-delta.lambda.i/sin(theta * deg2rad)) # transmission rate of the ith layer at theta
    e.atm.lambda.layers[i] <-1-tau1.lambda.layers[i] 
    # transmission rate from the middle of ith layer to the surface
    # consider the middle point of each layer as its equivalent height 
    tau1.lambda[i] <- cal.transmission((i-0.5)*layer.depth/sin(theta * deg2rad), theta, layer.depth)
    # ith layer's irradiance reaching at surface vertically = 
    # black body irradiance * ith layer's emissivity * transmission rate to the ground 
    L.atm.lambda.layers[i] <- (BM.temp.lambda.atm/pi)*e.atm.lambda.layers[i]*tau1.lambda[i]   # Equation (3.78) (4.53)
  }
  L.atm.lambda <- sum(L.atm.lambda.layers)  # atmopshere radiance on the target
  return(list(tau1.lambda.layers, tau1.lambda, e.atm.lambda.layers, L.atm.lambda.layers, L.atm.lambda))
}

# Let's do some test -----------------------------------------------------
wv <- lambda.ls[2]
try <- lapply(theta.ls, function(x) cbind(data.frame(cal.atmosphere.emiss.oneangle(wv,x)[1:4]), 
                                          x, cal.atmosphere.emiss.oneangle(wv,x)[[5]]))
try <- lapply(try, setNames, c("tau1.layer","tau1.underlayer","emiss.layer", "radiance.layer", "theta", "total.radiance"))
try <- do.call(rbind, try)
# given the same range
plot(theta.ls, unique(try$total.radiance), xlab= expression("Elevation Angles"), 
     ylab = expression(paste("(W m-2 sr-1 ", mu, "m-1)")))

#
U.atm.lambda.oneangle <- function (wv, theta) {cal.atmosphere.emiss.oneangle(wv, theta)[[5]]}
U.atm.lambda.oneangle(wv, 60)
# integrate(E.atm.lambda.oneangle, lower = 30, upper =90)


# Section 4: Calculate Radiance and Transmission --------------------------
cal.radiance <- function(emiss, temp, lambda, theta, range, layer.depth=500){
  # calculate transmission tau
  # calculate the number of parallel atmosphere layers, each layer is 1000(m)
  # sum delta.lambda for each layer, Equation (3.108)
  tau2.lambda <- cal.transmission(range, theta, layer.depth)
  
  # claculate the self-emitted spectral radiance
  # http://www.spectralcalc.com/blackbody/units_of_wavelength.html
  # Black body thermal emission
  BM.temp.lambda <-blackbody_radiance_calculation(lambda, temp)
  BL.temp.lambda <- BM.temp.lambda/pi   # Equation (3.78), [W m-2 sr-1 um-1]
  # Gray body thermal emission
  L.T.lambda <- emiss * BL.temp.lambda * tau2.lambda
  
  # Downwelled self-emitted irradiance, Equation (4.53) & (4.55) 
  # The way to calculate the total atmosphere radiance (comments under Equation 4.49) is summing over 
  # (the emissiviry * irradiance of each atmosphere layer * transmission rate under this layer) for each theta and phi 
  # Then integrate the shape factor of the sky clear to the object (i.e., 2*pi).
  # Here we assume the atmosphere temperature and gas concentration are both constant for each layer
  # Approximately 63 % of LDR is due to emission within the lowest 100 m of the atmosphere, less than 5 % 
  
  # originates from layers above an altitude of 2 km altitude so we profile the atmophere atitude as 2km
  # The atmosphere transmission rate = 1- emissivity
  reflec <- 1-emiss
  
  # "vertically" reaching at the target*cos(theta*deg2rad); d\omega = sin(theta*deg2rad)d(theta)d(phi) from one angle
  # integrate from 0~90, 0~360 degree
  # Equations (3.71) (3.73) and (4.55)
  # E.downatm.lambda<- integrate(...)
  # L.E.lambda <- tau2.lambda*reflec*E.downatm.lambda/pi
  
  # Upwelled self-emitted irradiance
  L.U.lambda<- U.atm.lambda.oneangle(lambda, theta)
  
  # Total radiance
  L.total <- L.T.lambda + L.U.lambda
  return(list(L.T.lambda, L.U.lambda, L.total, tau2.lambda))
} 

mat <- expand.grid(emiss.ls, temp.ls, lambda.ls, theta.ls, range.ls) # data.frame
colnames(mat) <- c("emissivity","temperature","lambda","theta","range")
result.list <- list()
for ( i in 1:nrow(mat) ) {
  # print(paste(i,"of",nrow(mat)))
  tmp <- cal.radiance(mat[i,1],mat[i,2],mat[i,3],mat[i,4],mat[i,5])
  result.list[[i]]<- cbind(tmp[[1]], tmp[[2]],tmp[[3]], tmp[[4]])
}
result.df<- as.data.frame(do.call(rbind, result.list))
colnames(result.df) <- c("Target", "Atmosphere", "Total", "Transmission")

component <- c("T","U","total")  # Target self-emitted radiance, Upwelling Atmospheric self-emitted radiance or total radiance
radiance.check <- function(comp){
  index <- which(component == comp)
  radi.vec = result.df[,index]
  tau2.vec = result.df[,4]
  return(list(radi.vec, tau2.vec))
}


# Section 5: Check Correlations -------------------------------------------
# Which component of radicance do you wanna check?
radiance.comp <- radiance.check("U")
mat[,"radiance"] = radiance.comp[[1]]
mat[,"transmission"] = radiance.comp[[2]]
# correlation calculation
cormat <- round(cor(mat[,c(6,7)],mat[,1:5]),7)
# choose subset whose emissivity is 1
mat.emiss1 <-  mat[mat[,1]==1,][,-1]
mat.emiss1[,5] <- round(mat.emiss1[,5],7)   # round radiance to 7 digits
mat.emiss1[,6] <- round(mat.emiss1[,6],7)   # round transmission to 7 digits
mat.emiss1 <- unique(mat.emiss1) 

# Section 6: Define Colors for Plots --------------------------------------------------------
# pie(rep(1,720), col = rev(matlab.like2(n)), border = NA,label=NA)
drawggplot <- function (data, xlab, ylab, units,legend.posi, legend.labels, n, transp){
  ggplot(data, aes_string(x=xlab, y=ylab, color="id", group= "id"))+
    geom_point() + geom_line(alpha=transp) + scale_colour_gradientn(colours=rev(matlab.like2(n)),breaks = c(n/6,5*n/6),labels = legend.labels) + 
    labs(x=paste(str_to_sentence(xlab), units[1], sep = " "), y=str_to_sentence(paste(str_to_sentence(ylab), units[2], sep=" ")), color=ylab)+
    theme(axis.title.x = element_text(size=10, face="bold"), axis.title.y = element_text(size=10, face="bold"), 
          legend.position = legend.posi, legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 8),
          legend.background = element_rect(fill="transparent"))
}

draw <- function (xname, yname, units, position, labels, transp=1, tmp=mat.emiss1){ 
  if (xname == "temperature"){
    tmp <- data.frame(tmp %>% group_by(range,lambda,theta) %>% summarise(temperature = paste(temperature, collapse =","), 
                                                              radiance= paste(radiance,collapse =","), transmission = paste(transmission, collapse =",")))
    tmp$temperature <- lapply(tmp$temperature, function(x) as.numeric(unlist(str_split(x, ','))))
    tmp$radiance <- lapply(tmp$radiance, function(x) as.numeric(unlist(str_split(x, ','))))
    tmp$transmission <- lapply(tmp$transmission, function(x) as.numeric(unlist(str_split(x, ','))))
    n <- nrow(tmp)
    if (yname == "radiance"){
      tmp <- tmp[order(-tmp[,1], -tmp[,2],tmp[,3]),]      # order the variables based on their correlations with the radiance
      tmp$id <- 1:n               # must be numeric rather than character
      points.plot<- unnest(tmp)
      g<- drawggplot(points.plot, "temperature","radiance", units, position, labels, n, transp)
      }else{
        tmp <- tmp[order(-tmp[,1]),]    # order the variables based on their correlations with the transmission rate
        tmp$id <- 1:n               # must be numeric rather than character
        points.plot<- unnest(tmp)
        g<- drawggplot(points.plot, "temperature","transmission",units, position, labels, n, transp)
      }
  }
  if (xname == "lambda"){
    tmp <- data.frame(tmp %>% group_by(temperature,range,theta) %>% summarise(lambda = paste(lambda, collapse =","), 
                                                                                 radiance= paste(radiance,collapse =","), transmission = paste(transmission, collapse =",")))
    tmp$lambda <- lapply(tmp$lambda, function(x) as.numeric(unlist(str_split(x, ','))))
    tmp$radiance <- lapply(tmp$radiance, function(x) as.numeric(unlist(str_split(x, ','))))
    tmp$transmission <- lapply(tmp$transmission, function(x) as.numeric(unlist(str_split(x, ','))))
    n <- nrow(tmp)
    if (yname == "radiance"){
      tmp <- tmp[order(tmp[,1], -tmp[,2], tmp[,3]),]      # order the variables based on their correlations with the radiance
      tmp$id <- 1:n               # must be numeric rather than character
      points.plot<- unnest(tmp)
      g <- drawggplot(points.plot, "lambda","radiance", units, position,labels, n, transp)
    }else{
      tmp <- tmp[order(-tmp[,2]),]    # order the variables based on their correlations with the transmission rate
      tmp$id <- 1:n               # must be numeric rather than character
      points.plot<- unnest(tmp)
      g<- drawggplot(points.plot, "lambda","transmission",units, position,labels, n,transp)
    }
  }
  if (xname == "theta"){
    tmp <- as.data.frame(tmp %>% group_by(temperature,range,lambda) %>% summarise(theta = paste(theta, collapse =","), 
                                                                   radiance= paste(radiance,collapse =","), transmission = paste(transmission, collapse =",")))
    tmp$theta <- lapply(tmp$theta, function(x) as.numeric(unlist(str_split(x, ','))))
    tmp$radiance <- lapply(tmp$radiance, function(x) as.numeric(unlist(str_split(x, ','))))
    tmp$transmission <- lapply(tmp$transmission, function(x) as.numeric(unlist(str_split(x, ','))))
    n <- nrow(tmp)
    if (yname == "radiance"){
      tmp <- tmp[order(tmp[,1], -tmp[,2],-tmp[,3]),]      # order the variables based on their correlations with the radiance
      tmp$id <- 1:n               # must be numeric rather than character
      points.plot<- unnest(tmp)
      g<- drawggplot(points.plot, "theta","radiance",units, position,labels, n, transp)
    }else{
      tmp <- tmp[order(-tmp[,2]),]    # order the variables based on their correlations with the transmission rate
      tmp$id <- 1:n               # must be numeric rather than character
      points.plot<- unnest(tmp)
      g<- drawggplot(points.plot, "theta","transmission",units, position, labels, n, transp)
    }
  }
  if (xname == "range"){
    tmp <- as.data.frame(tmp %>% group_by(temperature,lambda,theta) %>% summarise(range = paste(range, collapse =","), 
                                                                    radiance= paste(radiance,collapse =","), transmission = paste(transmission, collapse =",")))
    tmp$range <- lapply(tmp$range, function(x) as.numeric(unlist(str_split(x, ','))))
    tmp$radiance <- lapply(tmp$radiance, function(x) as.numeric(unlist(str_split(x, ','))))
    tmp$transmission <- lapply(tmp$transmission, function(x) as.numeric(unlist(str_split(x, ','))))
    n <- nrow(tmp)
    if (yname == "radiance"){
      tmp <- tmp[order(tmp[,1], -tmp[,2], tmp[,3]),]      # order the variables based on their correlations with the radiance
      tmp$id <- 1:n               # must be numeric rather than character
      points.plot<- unnest(tmp)
      g<- drawggplot(points.plot, "range","radiance",units, position, labels,n,transp)
    }else{
      tmp$id <- 1:n               # must be numeric rather than character
      points.plot<- unnest(tmp)
      g<- drawggplot(points.plot, "range","transmission",units, position, labels, n, transp)
    }
  }
  return(g)
}


# Section 7: Draw Plots ---------------------------------------------------
### plot of radiance
# xlab, ylab, units for xlab and ylab, legend location, legend text 
g1 <- draw("temperature","radiance", c("(Kelvin)","(W m-2 sr-1 um-1)"), c(0.25, 0.8), 
           c("Larger Range, Lambda \nand Smaller Temperature", "Smaller Range, Lambda \nand Larger Temperature"), 0.1)
g2 <- draw("lambda","radiance", c("(um)","(W m-2 sr-1 um-1)"), "none", NULL)
g3 <- draw("theta","radiance",c("(degree)","(W m-2 sr-1 um-1)"), "none", NULL)
g4 <- draw("range","radiance",c("(m)","(W m-2 sr-1 um-1)"), "none", NULL)
jpeg("radiance.jpg", width = 6000, height = 6000, units = "px", quality = 100, res=600)
grid.arrange(g1, g2, g3, g4, ncol = 2, nrow=2,
             widths = c(5, 5), top = textGrob("Self-emitted Radiance of Graybody objects for Long-wave Length", gp=gpar(fontface="bold")))
dev.off()

jpeg("gg.jpg", width = 6000, height = 6000, units = "px", quality = 100, res=600)
  draw("lambda","radiance", c("(um)","(W m-2 sr-1 um-1)"),    c(.85, .85), 
            c("Larger Range \nand Smaller Temperature", "Smaller Range \nand Larger Temperature"), 0.1)
dev.off()


g5 <- draw("temperature","transmission",c("(Kelvin)",""),"none", NULL)
g6 <- draw("lambda","transmission",c("(um)",""),"none", NULL)
g7 <- draw("theta","transmission",c("(degree)",""),"none", NULL)
g8 <- draw("range","transmission",c("(m)",""),c(0.8, 0.8),c("Larger Range", "Smaller Range"))
jpeg("tau2.jpg", width = 6000, height = 6000, units = "px", quality = 100, res=600)
grid.arrange(g5, g6, g7, g8, ncol = 2, nrow=2,
             widths = c(5, 5), top = textGrob("Transmission Rate of Graybody objects for Long-wave Length",gp=gpar(fontface="bold")))
dev.off()


### plot of transmission rate
# define.cols <- function(xlab, ylab){
#   binlength <- length(ylab)/length(xlab)
#   rgb.palette <- rev(matlab.like2(binlength))
# }
# jpeg('tau2Lines.jpg',width = 10, height = 10, units = 'in', res = 300)
# par(mfrow=c(2,2))
# plot(mat.emiss1[order(mat.emiss1[,2],-mat.emiss1[,5]),c(2,7)], main="Order by decreasing range", xlab="Temperature (Kelvin)", ylab="Transmission Rate", col=define.cols(temp.ls, mat.emiss1[,7]))
# draw("temperature","transmission")
# plot(mat.emiss1[order(mat.emiss1[,3],-mat.emiss1[,5]),c(3,7)], main="Order by decreasing range", xlab="Lambda (um)", ylab="Transmission Rate", col=define.cols(lambda.ls, mat.emiss1[,7]))
# draw("lambda","transmission")
# plot(mat.emiss1[order(mat.emiss1[,4],-mat.emiss1[,5]),c(4,7)], main="Order by decreasing range", xlab="Theta (degree)", ylab="Transmission Rate", col=define.cols(theta.ls, mat.emiss1[,7]))
# draw("theta","transmission")
# plot(mat.emiss1[order(mat.emiss1[,5]),c(5,7)],main="Natural Ordering", xlab="Range (m)", ylab="Transmission Rate", col=define.cols(range.ls, mat.emiss1[,7]))
# draw("range","transmission")
# dev.off()



# Section 8: Plots of Hypersurface -----------------------------------------
filter<- mat[which(mat[,1]==1 & mat[,4]==0 & mat[,5]==3000),c(2,3,6,7)] # change theta
filter <- filter[order(filter[,1],filter[,2]),]
z <- acast(filter, temperature~lambda, value.var="radiance")
persp(temp.ls, lambda.ls, z, phi = 0, theta = 0,
        xlab = "Temperature (Kelvin)", ylab = "Lambda (um)",
        main = "Radiance Surface")

wireframe(Height ~ x*y, data = elevation.fit,
          xlab = "X Coordinate", ylab = "Y Coordinate",
          main = "Surface Radiance",
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -60, x = -60)
)





