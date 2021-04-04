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


# display.brewer.all()
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(grDevices)
library(viridis) 
library(dplyr)

source("Z_general_functions.R")
source("Z_global_variables.R")


### Functions for Generating Plots
# use facet wrap to generate the data
boxplotfunc <- function(plotdata,filename){
  theme.boxplot <- theme(axis.text.x = element_text(colour = "grey20", size = 12),
                         axis.text.y = element_text(colour = "grey20", size = 12),
                         axis.title.x = element_text(size=12, face="bold"),
                         axis.title.y = element_text(size=12, face="bold"),
                         strip.text = element_text(size = 12, face ="bold"))
  jpeg(filename, width= 3000, height=4500, units = "px", quality = 100, res=300)
  ggplot(plotdata, mapping = aes(x = wavelength, y = value, group= wavelength))+
    facet_wrap(variable~., scales = "free_y", labeller = as_labeller(label), ncol=1) + 
    geom_boxplot(fill = "white", outlier.color= "red", outlier.fill = "red", outlier.shape=8, outlier.size=1) + 
    stat_boxplot(geom ="errorbar") + geom_point(size = 1, color="black", stroke=0) + theme.boxplot + 
    labs(x=expression(paste("Wavelength (", mu, "m)")), y=expression(paste("(W cm-2 sr-1 ", mu, "m-1)")))
  dev.off()
}

# highly customized plots
# dates, integer are considered as continous data. Discrete values cannot be applied to gradient colors
# variable is the column name of the plotdata which decides the color of points and the legend
boxplotfunc.scan <- function(plotdata, variable, layout, imgwidth, imgheight, boxcolor="black",boxlwd=0.2, title=NULL, filename="", fixScatteringSolar=F, 
                                  legend.pos=c(0.005,0.99) ,legend.com = T, legend.continous=F){
  # whether the legend should be continous or not
  if(legend.continous){
    sc <- scale_colour_gradientn(name= variable, colours=boxplot.color)
    gp <- geom_point(aes_string(colour = variable), size = 0.5, stroke=0)
    theme.boxplot <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"), 
                           legend.position=c(0.005,0.99), legend.justification = c("left", "top"), legend.title=element_text(size=12, face = "bold"), legend.text = element_text(size=8),
                           legend.key.height= unit(0.4, "cm"), legend.key.width=unit(0.3, "cm"), legend.background = element_blank(),
                           legend.box.margin = margin(0, 0, 0, 0))
    
  }else{
    theme.boxplot <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"),
                           legend.position=legend.pos, legend.justification = c("left", "top"), legend.title=element_text(size=12, face = "bold"), legend.text = element_text(size=8),
                           legend.key.size = unit(0, "cm"), legend.background = element_blank(),
                           legend.box.margin = margin(0, 0, 0, 0))
    colorlength<- length(unique(plotdata[,variable]))
    sc <- scale_colour_manual(name= variable, values=colorRampPalette(c("green","yellow","orange"))(colorlength)) +  
      guides(colour = guide_legend(override.aes = list(size = 2),ncol = ifelse(colorlength>8, 3, 2), byrow = TRUE,reverse = TRUE))
    gp <- geom_point(aes_string(colour = sprintf("factor(%s)",variable)), size = 0.5, stroke=0)
  }
  # whether the legend should be plotted for every component
  if(legend.com){
    theme.boxplot1 <- theme.boxplot
  }else{
    theme.boxplot1 <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"),
                            legend.position="none")
  }
  # set the xlimits
  xlim <- c(min(plotdata$wavelength),max(plotdata$wavelength))
  # same ylimit for the surface_emiss and path_emiss
  ylimits1<- c(min(c(plotdata$path_emiss,plotdata$surface_emiss)),max(c(plotdata$path_emiss,plotdata$surface_emiss)))
  # same ylimit for the solar multiple scattering, solar single scattering and 
  # whether to set a universal same ylimits across different plots
  if(fixScatteringSolar){
    ylimits2<-  c(0,4.5e-09)
  }else{
    ylimits2<- c(min(c(plotdata$sing_scat,plotdata$path_multi_scat)),max(c(plotdata$sing_scat,plotdata$path_multi_scat)))
  }
  # plot the components
  g1<- ggplot(plotdata, mapping = aes(x = wavelength, y = total_rad, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks(4))+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd, outlier.shape = NA) + labs(title= "Total radiance at sensor")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd) + gp + theme_bw() +theme.boxplot+ sc+ expand_limits(x=xlim,y = 0) 
  g2<- ggplot(plotdata, mapping = aes(x = wavelength, y = path_emiss, group= wavelength))+ scale_y_continuous(limits=ylimits1, breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) +labs(title= "Path thermal emission radiance") + 
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc + expand_limits(x=xlim,y = 0) 
  g3<- ggplot(plotdata, mapping = aes(x = wavelength, y = surface_emiss, group= wavelength))+ scale_y_continuous(limits=ylimits1,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Surface thermal emission radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw()  + theme.boxplot1+ sc + expand_limits(x=xlim, y = 0) 
  g4<- ggplot(plotdata, mapping = aes(x = wavelength, y = grnd_rflt, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks(3))+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Ground reflected radiance")+ 
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc + expand_limits(x=xlim,y = 0) 
  g5<- ggplot(plotdata, mapping = aes(x = wavelength, y = path_thermal_scat, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Path thermal scattering radiance")+ 
    stat_boxplot(geom ="errorbar",color =boxcolor, lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc + expand_limits(x=xlim,y = 0)
  g6<- ggplot(plotdata, mapping = aes(x = wavelength, y = sing_scat, group= wavelength))+ scale_y_continuous(limits=ylimits2,breaks=pretty_breaks(2))+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Single scattering solar radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc + expand_limits(x=xlim,y = 0)
  g7<- ggplot(plotdata, mapping = aes(x = wavelength, y = path_multi_scat, group= wavelength))+ scale_y_continuous(limits=ylimits2,breaks=pretty_breaks(2))+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Multiple scattering solar radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc+ expand_limits(x=xlim,y = 0)
  
  plotgrobs <- list(g1,g2,g3,g4,g5,g6,g7)[sort(unique(as.vector(layout)))]
  if(filename==""){
    grid.arrange(grobs=plotgrobs, layout_matrix = layout, 
                 top =textGrob(title, gp = gpar(fontface = "bold", cex = 2)),
                 left =textGrob(expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1))) 
  }else{
    jpeg(filename, width= imgwidth, height=imgheight, units = "px", quality = 100, res=300)
    grid.arrange(grobs=plotgrobs, layout_matrix = layout, 
                 top =textGrob(title, gp = gpar(fontface = "bold", cex = 2)),
                 left =textGrob(expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1))
    ) 
    dev.off()
  }
}

### Fix phi, plot different theta, refer to comments in boxplotfunc.scan() above
boxplotfunctheta.scan <- function(plotdata, layout, imgwidth, imgheight, boxcolor="black",boxlwd=0.2, title=NULL, filename="", fixScatteringSolar=F, 
                                  legend.pos=c(0.005,0.99) ,legend.com = T){
  theme.boxplot <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"),
                         legend.position=legend.pos, legend.justification = c("left", "top"), legend.title=element_text(size=12, face = "bold"), legend.text = element_text(size=8),
                         legend.key.size = unit(0, "cm"), legend.background = element_blank(),
                         legend.box.margin = margin(0, 0, 0, 0))
  if(legend.com){
    theme.boxplot1 <- theme.boxplot
  }else{
    theme.boxplot1 <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"),
                            legend.position="none")
  }
  gp <- geom_point(aes(colour = factor(theta)), size = 0.5, stroke=0)
  sc <- scale_colour_manual(name= expression(theta), values=boxplot.color)
  #sc <- scale_colour_manual(name=NULL, values=boxplot.color)
  gl<- guides(colour = guide_legend(override.aes = list(size = 2), title.position="left", nrow = 2, byrow = TRUE,reverse = TRUE))
  xlim <- c(min(plotdata$wavelength),max(plotdata$wavelength))
  ylimits1<- c(min(c(plotdata$path_emiss,plotdata$surface_emiss)),max(c(plotdata$path_emiss,plotdata$surface_emiss)))
  if(fixScatteringSolar){
    ylimits2<-  c(0,4.7e-09)
  }else{
    ylimits2<- c(min(c(plotdata$sing_scat,plotdata$path_multi_scat)),max(c(plotdata$sing_scat,plotdata$path_multi_scat)))
  }
  g1<- ggplot(plotdata, mapping = aes(x = wavelength, y = total_rad, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks(4))+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd, outlier.shape = NA) + labs(title= "Total radiance at sensor")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd) + gp + theme_bw() +theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0) 
  g2<- ggplot(plotdata, mapping = aes(x = wavelength, y = path_emiss, group= wavelength))+ scale_y_continuous(limits=ylimits1, breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) +labs(title= "Path thermal emission radiance") + 
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot+ sc+ gl+ expand_limits(x=xlim,y = 0) 
  g3<- ggplot(plotdata, mapping = aes(x = wavelength, y = surface_emiss, group= wavelength))+ scale_y_continuous(limits=ylimits1,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Surface thermal emission radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw()  + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim, y = 0) 
  g4<- ggplot(plotdata, mapping = aes(x = wavelength, y = grnd_rflt, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks(3))+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Ground reflected radiance")+ 
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0) 
  g5<- ggplot(plotdata, mapping = aes(x = wavelength, y = path_thermal_scat, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Path thermal scattering radiance")+ 
    stat_boxplot(geom ="errorbar",color =boxcolor, lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0)
  g6<- ggplot(plotdata, mapping = aes(x = wavelength, y = sing_scat, group= wavelength))+ scale_y_continuous(limits=ylimits2,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Single solar scattering radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0)
  g7<- ggplot(plotdata, mapping = aes(x = wavelength, y = path_multi_scat, group= wavelength))+ scale_y_continuous(limits=ylimits2,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Multiple solar scattering radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0)
  
  plotgrobs <- list(g1,g2,g3,g4,g5,g6,g7)[sort(unique(as.vector(layout)))]
  if(filename==""){
    grid.arrange(grobs=plotgrobs, layout_matrix = layout, 
                 top =textGrob(title, gp = gpar(fontface = "bold", cex = 1.5)),
                 left =textGrob(expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1))) 
  }else{
    jpeg(filename, width= imgwidth, height=imgheight, units = "px", quality = 100, res=300)
    grid.arrange(grobs=plotgrobs, layout_matrix = layout, 
                 top =textGrob(title, gp = gpar(fontface = "bold", cex = 1.5)),
                 left =textGrob(expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1))
    ) 
    dev.off()
  }
}
### Fix phi, plot different theta, refer to comments in boxplotfunc.scan() above
boxplotfunctheta.digit <- function(plotdata, layout, imgwidth, imgheight, boxcolor="black",boxlwd=0.3, title=NULL, filename="", fixScatteringSolar=F, legend.com = T){
  theme.boxplot <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"),
                         legend.position=c(0.995,0.99), legend.justification = c("right", "top"), legend.title=element_text(size=12, face = "bold"), legend.text = element_text(size=8),
                         legend.key.size = unit(0, "cm"), legend.background = element_blank(),
                         legend.box.margin = margin(0, 0, 0, 0))
  if(legend.com){
    theme.boxplot1 <- theme.boxplot
  }else{
    theme.boxplot1 <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"),
                            legend.position="none")
  }
  gp <- geom_point(aes(colour = factor(theta)), size = 0.5, stroke=0)
  sc <- scale_colour_manual(name= expression(theta), values=boxplot.color)
  gl<- guides(colour = guide_legend(override.aes = list(size = 2),ncol = 3, byrow = TRUE,reverse = TRUE))
  xlim <- c(min(plotdata$wavenumber),max(plotdata$wavenumber))
  ylimits1<- c(min(c(plotdata$path_emiss,plotdata$surface_emiss)),max(c(plotdata$path_emiss,plotdata$surface_emiss)))
  if(fixScatteringSolar){
    ylimits2<-  c(0,4.4e-09)
  }else{
    ylimits2<- c(min(c(plotdata$sing_scat,plotdata$path_multi_scat)),max(c(plotdata$sing_scat,plotdata$path_multi_scat)))
  }
  g1<- ggplot(plotdata, mapping = aes(x = wavenumber, y = total_rad, group= wavenumber))+ scale_y_continuous(breaks=pretty_breaks(4))+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd, outlier.shape = NA) + labs(title= "Total radiance at sensor")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd) + gp + theme_bw() +theme.boxplot+ sc+ gl+ expand_limits(x=xlim,y = 0) 
  g2<- ggplot(plotdata, mapping = aes(x = wavenumber, y = path_emiss, group= wavenumber))+ scale_y_continuous(limits=ylimits1, breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) +labs(title= "Path thermal emission radiance") + 
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0) 
  g3<- ggplot(plotdata, mapping = aes(x = wavenumber, y = surface_emiss, group= wavenumber))+ scale_y_continuous(limits=ylimits1,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Surface thermal emission radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw()  + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0) 
  g4<- ggplot(plotdata, mapping = aes(x = wavenumber, y = grnd_rflt, group= wavenumber))+ scale_y_continuous(breaks=pretty_breaks(3))+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Ground reflected radiance")+ 
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0) 
  g5<- ggplot(plotdata, mapping = aes(x = wavenumber, y = path_thermal_scat, group= wavenumber))+ scale_y_continuous(breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Path thermal scattering radiance")+ 
    stat_boxplot(geom ="errorbar",color =boxcolor, lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0)
  g6<- ggplot(plotdata, mapping = aes(x = wavenumber, y = sing_scat, group= wavenumber))+ scale_y_continuous(limits=ylimits2,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Single solar scattering radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0)
  g7<- ggplot(plotdata, mapping = aes(x = wavenumber, y = path_multi_scat, group= wavenumber))+ scale_y_continuous(limits=ylimits2,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Multiple solar scattering radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+ sc+ gl+ expand_limits(x=xlim,y = 0)
  
  plotgrobs <- list(g1,g2,g3,g4,g5,g6,g7)[sort(unique(as.vector(layout)))]
  if(filename==""){
    grid.arrange(grobs=plotgrobs, layout_matrix = layout, 
                 top =textGrob(title, gp = gpar(fontface = "bold", cex = 2)),
                 left =textGrob(expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1))) 
  }else{
    jpeg(filename, width= imgwidth, height=imgheight, units = "px", quality = 100, res=300)
    grid.arrange(grobs=plotgrobs, layout_matrix = layout, 
                 top =textGrob(title, gp = gpar(fontface = "bold", cex = 2)),
                 left =textGrob(expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1))
    ) 
    dev.off()
  }
}

### Fix theta, plot different phi, refer to comments in boxplotfunc.scan() above
boxplotfuncphi.scan <- function(plotdata, layout, imgwidth, imgheight, boxcolor="black",boxlwd=0.2, title=NULL, filename="", legend.com = T){
  theme.boxplot <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"), 
                         legend.position=c(0.005,0.99), legend.justification = c("left", "top"), legend.title=element_text(size=12, face = "bold"), legend.text = element_text(size=8),
                         legend.key.height= unit(0.3, "cm"), legend.key.width=unit(0.2, "cm"), legend.background = element_blank(),
                         legend.box.margin = margin(0, 0, 0, 0))
  if(legend.com){
    theme.boxplot1 <- theme.boxplot
  }else{
    theme.boxplot1 <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"),
                            legend.position="none")
  }
  gp <- geom_point(aes(colour = phi), size = 0.5, stroke=0)
  sc <- scale_colour_gradientn(name= expression(phi),colours=boxplot.color,breaks=c(0,90,180,270,360))
  xlim <- c(min(plotdata$wavelength),max(plotdata$wavelength))
  ylimits1<- c(min(c(plotdata$path_emiss,plotdata$surface_emiss)),max(c(plotdata$path_emiss,plotdata$surface_emiss)))
  ylimits2<- c(min(c(plotdata$sing_scat,plotdata$path_multi_scat)),max(c(plotdata$sing_scat,plotdata$path_multi_scat)))
  g1<- ggplot(plotdata, mapping = aes(x = wavelength, y = total_rad, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks(4))+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd, outlier.shape = NA) + labs(title= "Total radiance at sensor")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd) + gp + theme_bw() +theme.boxplot+ sc+expand_limits(x=xlim,y = 0) 
  g2<- ggplot(plotdata, mapping = aes(x = wavelength, y = path_emiss, group= wavelength))+ scale_y_continuous(limits=ylimits1,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) +labs(title= "Path thermal emission radiance") + 
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot1+sc+expand_limits(x=xlim,y = 0) 
  g3<- ggplot(plotdata, mapping = aes(x = wavelength, y = surface_emiss, group= wavelength))+ scale_y_continuous(limits=ylimits1,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Surface thermal emission radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() +theme.boxplot1+sc+expand_limits(x=xlim,y = 0) 
  g4<- ggplot(plotdata, mapping = aes(x = wavelength, y = grnd_rflt, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks(2))+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Ground reflected radiance")+ 
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() +theme.boxplot1 +sc+expand_limits(x=xlim,y = 0)  
  g5<- ggplot(plotdata, mapping = aes(x = wavelength, y = path_thermal_scat, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Path thermal scattering radiance")+ 
    stat_boxplot(geom ="errorbar",color =boxcolor, lwd=boxlwd/2) + gp + theme_bw() +theme.boxplot1 +sc+expand_limits(x=xlim,y = 0)
  g6<- ggplot(plotdata, mapping = aes(x = wavelength, y = sing_scat, group= wavelength))+ scale_y_continuous(limits=ylimits2,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Single solar scattering radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() +theme.boxplot1 +sc+expand_limits(x=xlim,y = 0) 
  g7<- ggplot(plotdata, mapping = aes(x = wavelength, y = path_multi_scat, group= wavelength))+ scale_y_continuous(limits=ylimits2,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Multiple solar scattering radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() +theme.boxplot1 + sc+expand_limits(x=xlim,y = 0) 
  
  plotgrobs <- list(g1,g2,g3,g4,g5,g6,g7)[sort(unique(as.vector(layout)))]
  if(filename==""){
    grid.arrange(grobs=plotgrobs, layout_matrix = layout, 
                 top =textGrob(title, gp = gpar(fontface = "bold", cex = 1.5)),
                 left =textGrob(expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1))) 
  }else{
    jpeg(filename, width= imgwidth, height=imgheight, units = "px", quality = 100, res=300)
    grid.arrange(grobs=plotgrobs, layout_matrix = layout, 
                 top =textGrob(title, gp = gpar(fontface = "bold", cex = 1.5)),
                 left =textGrob(expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1))
    ) 
    dev.off()
  }
}
### Fix theta, plot different phi, refer to comments in boxplotfunc.scan() above
boxplotfuncphi.digit <- function(plotdata, layout, imgwidth, imgheight, boxcolor="black",boxlwd=0.3, title=NULL, filename=""){
  theme.boxplot <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"), 
                         legend.position=c(0.995,0.99), legend.justification = c("right", "top"), legend.title=element_text(size=12, face = "bold"), legend.text = element_text(size=8),
                         legend.key.height= unit(0.4, "cm"), legend.key.width=unit(0.3, "cm"), legend.background = element_blank(),
                         legend.box.margin = margin(0, 0, 0, 0))
  gp <- geom_point(aes(colour = phi), size = 0.5, stroke=0)
  sc <- scale_colour_gradientn(name= expression(phi),colours=boxplot.color,breaks=c(0,90,180,270,360))
  xlim <- c(min(plotdata$wavenumber),max(plotdata$wavenumber))
  ylimits1<- c(min(c(plotdata$path_emiss,plotdata$surface_emiss)),max(c(plotdata$path_emiss,plotdata$surface_emiss)))
  ylimits2<- c(min(c(plotdata$sing_scat,plotdata$path_multi_scat)),max(c(plotdata$sing_scat,plotdata$path_multi_scat)))
  g1<- ggplot(plotdata, mapping = aes(x = wavenumber, y = total_rad, group= wavenumber))+ scale_y_continuous(breaks=pretty_breaks(4))+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd, outlier.shape = NA) + labs(title= "Total radiance at sensor")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd) + gp + theme_bw() +theme.boxplot+ sc+expand_limits(x=xlim,y = 0) 
  g2<- ggplot(plotdata, mapping = aes(x = wavenumber, y = path_emiss, group= wavenumber))+ scale_y_continuous(limits=ylimits1,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) +labs(title= "Path thermal emission radiance") + 
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() + theme.boxplot+sc+expand_limits(x=xlim,y = 0) 
  g3<- ggplot(plotdata, mapping = aes(x = wavenumber, y = surface_emiss, group= wavenumber))+ scale_y_continuous(limits=ylimits1,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Surface thermal emission radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() +theme.boxplot+sc+expand_limits(x=xlim,y = 0) 
  g4<- ggplot(plotdata, mapping = aes(x = wavenumber, y = grnd_rflt, group= wavenumber))+ scale_y_continuous(breaks=pretty_breaks(2))+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Ground reflected radiance")+ 
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() +theme.boxplot +sc+expand_limits(x=xlim,y = 0)  
  g5<- ggplot(plotdata, mapping = aes(x = wavenumber, y = path_thermal_scat, group= wavenumber))+ scale_y_continuous(breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor, lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Path thermal scattering radiance")+ 
    stat_boxplot(geom ="errorbar",color =boxcolor, lwd=boxlwd/2) + gp + theme_bw() +theme.boxplot +sc+expand_limits(x=xlim,y = 0)
  g6<- ggplot(plotdata, mapping = aes(x = wavenumber, y = sing_scat, group= wavenumber))+ scale_y_continuous(limits=ylimits2,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Single solar scattering radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() +theme.boxplot +sc+expand_limits(x=xlim,y = 0) 
  g7<- ggplot(plotdata, mapping = aes(x = wavenumber, y = path_multi_scat, group= wavenumber))+ scale_y_continuous(limits=ylimits2,breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd/2, outlier.shape = NA) + labs(title= "Multiple solar scattering radiance")+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd/2) + gp + theme_bw() +theme.boxplot + sc+expand_limits(x=xlim,y = 0) 
  
  plotgrobs <- list(g1,g2,g3,g4,g5,g6,g7)[sort(unique(as.vector(layout)))]
  if(filename==""){
    grid.arrange(grobs=plotgrobs, layout_matrix = layout, 
                 top =textGrob(title, gp = gpar(fontface = "bold", cex = 2)),
                 left =textGrob(expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1))) 
  }else{
    jpeg(filename, width= imgwidth, height=imgheight, units = "px", quality = 100, res=300)
    grid.arrange(grobs=plotgrobs, layout_matrix = layout, 
                 top =textGrob(title, gp = gpar(fontface = "bold", cex = 2)),
                 left =textGrob(expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1))
    ) 
    dev.off()
  }
}

### Transmission rate, refer to comments in boxplotfunc.scan() above
boxplotfunctrans.comp.scan <- function(plotdata1, plotdata2, imgwidth, imgheight, boxcolor="black",boxlwd=0.2, title, filename=""){
  theme.boxplot <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"),
                         legend.position=c(0.005,0.99), legend.justification = c("left", "top"), legend.title=element_text(size=12, face = "bold"), legend.text = element_text(size=8),
                         legend.key.size = unit(0, "cm"), legend.background = element_blank(),
                         legend.box.margin = margin(0, 0, 0, 0))
  xlim1 <- c(min(plotdata1$wavelength),max(plotdata1$wavelength))
  xlim2 <- c(min(plotdata2$wavelength),max(plotdata2$wavelength))
  
  gp <- geom_point(aes(colour = factor(theta)), size = 0.5, stroke=0)
  sc <- scale_colour_manual(name= expression(theta), values=boxplot.color)
  gl<- guides(colour = guide_legend(override.aes = list(size = 2),ncol = 3, byrow = TRUE,reverse = TRUE))
  
  g1<- ggplot(plotdata1, mapping = aes(x = wavelength, y = path_trans, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd, outlier.shape = NA) + labs(title= expression("Tansmission Rate for"~theta~"= [30-90]"^degree~"with fixed range"))+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd) + gp + theme_bw() +theme.boxplot+ sc+ gl + expand_limits(x=xlim1,y = c(0,1))
  g2<- ggplot(plotdata2, mapping = aes(x = wavelength, y = path_trans, group= wavelength))+ scale_y_continuous(breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd, outlier.shape = NA) + labs(title= expression("Tansmission Rate for"~theta~"= [30-90]"^degree~"with fixed altitude"))+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd) + gp + theme_bw() +theme.boxplot+ sc+ gl +expand_limits(x=xlim2,y = c(0,1))
  if(filename==""){
    grid.arrange(g1,g2, ncol = 2, top =textGrob(title, gp = gpar(fontface = "bold", cex = 2)),
                 left =textGrob("Transmission Rate", rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1)))
  }else{
    jpeg(filename, width= imgwidth, height=imgheight, units = "px", quality = 100, res=300)
    grid.arrange(g1,g2, ncol = 2, top =textGrob(title, gp = gpar(fontface = "bold", cex = 2)),
                 left =textGrob("Transmission Rate", rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 1)))
    dev.off()
  }
}
### Transmission rate, refer to comments in boxplotfunc.scan() above
boxplotfunctrans.comp.digit <- function(plotdata1, plotdata2, imgwidth, imgheight, boxcolor="black",boxlwd=0.3, title, filename=""){
  theme.boxplot <- theme(axis.title=element_blank(), axis.text = element_text(size=12), plot.title = element_text(hjust=0.5, face = "bold",size=12, margin=margin(3,3,3,3)), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.3,0.3,0.3), "cm"),
                         legend.position=c(0.005,0.99), legend.justification = c("left", "top"), legend.title=element_text(size=12, face = "bold"), legend.text = element_text(size=8),
                         legend.key.size = unit(0, "cm"), legend.background = element_blank(),
                         legend.box.margin = margin(0, 0, 0, 0))
  
  xlim1 <- c(min(plotdata1$wavenumber),max(plotdata1$wavenumber))
  xlim2 <- c(min(plotdata2$wavenumber),max(plotdata2$wavenumber))
  gp <- geom_point(aes(colour = factor(theta)), size = 0.5, stroke=0)
  sc <- scale_colour_manual(name= expression(theta), values=boxplot.color)
  gl<- guides(colour = guide_legend(override.aes = list(size = 2),ncol = 3, byrow = TRUE,reverse = TRUE))
  
  g1<- ggplot(plotdata1, mapping = aes(x = wavenumber, y = path_trans, group= wavenumber))+ scale_y_continuous(breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd, outlier.shape = NA) + labs(title= expression("Tansmission Rate for"~theta~"= [30-90]"^degree~"with fixed range"))+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd) + gp + theme_bw() +theme.boxplot+ sc+ gl + expand_limits(x=xlim1,y = c(0,1))
  g2<- ggplot(plotdata2, mapping = aes(x = wavenumber, y = path_trans, group= wavenumber))+ scale_y_continuous(breaks=pretty_breaks())+
    geom_boxplot(fill = "white",color =boxcolor,lwd=boxlwd, outlier.shape = NA) + labs(title= expression("Tansmission Rate for"~theta~"= [30-90]"^degree~"with fixed altitude"))+
    stat_boxplot(geom ="errorbar",color =boxcolor,lwd=boxlwd) + gp + theme_bw() +theme.boxplot+ sc+ gl +expand_limits(x=xlim2,y = c(0,1))
  if(filename==""){
    grid.arrange(g1,g2, ncol = 2, top =textGrob(title, gp = gpar(fontface = "bold", cex = 2)),
                 left =textGrob("Transmission Rate", rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression("Wavenumber (cm-1)"), gp = gpar(fontface = "bold", cex = 1)))
  }else{
    jpeg(filename, width= imgwidth, height=imgheight, units = "px", quality = 100, res=300)
    grid.arrange(g1,g2, ncol = 2, top =textGrob(title, gp = gpar(fontface = "bold", cex = 2)),
                 left =textGrob("Transmission Rate", rot = 90, gp = gpar(fontface = "bold", cex = 1)),
                 bottom =textGrob(expression("Wavenumber (cm-1)"), gp = gpar(fontface = "bold", cex = 1)))
    dev.off()
  }
}


### Pie chart 
xu.pie <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, init.angle = if (clockwise) 90 else 0, 
                    density = NULL, angle = 45, col = NULL, border = NULL, lty = NULL, main = NULL, ...) {
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p), an=t2p)
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
      text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
           srt = ifelse(P$x < 0, P$an/pi*180+180, P$an/pi*180),
           adj = ifelse(P$x < 0, 1, 0), cex=0.5)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}

# offset sun location line over the pie chart
xu.pie.xy <- function( angle, r=1, offset= 2.5) {
  angle <- angle + offset
  return ( cbind( r*sin(angle*pi/180), r*cos(angle*pi/180) ) )
}

xu.filled.contour <-function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, length.out = ncol(z)), z, 
                              xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                              levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                              col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                              axes = TRUE, frame.plot = axes, plotKey=T, ...) {
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  
  # las: the style of axis labels.
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #record the expression given to be executed when the current function exits, 
  #useful for resetting graphical parameters or performing other cleanup actions
  on.exit(par(par.orig)) 
  
  if (plotKey==T) {
    # A double ("numeric") vector uses 8 bytes per element. An integer vector uses only 4 bytes per element. 1L means force to be integer
    w <- (3 + mar.orig[2L]) * par("csi") * 2.6   # the key plot width depends on the left margin they set
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w))) # matrix: [2,1]
    par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]+1 # right = left +1
    mar[2L] <- 0.5 # left becomes 0.5
    par(mar = mar) 
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
                yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col,border = "transparent")
    if (missing(key.axes)) {
      if (axes) 
        axis(4)
    }
    else key.axes
    box()
    if (!missing(key.title)) 
      key.title
  }
  mar <- mar.orig 
  mar[4L] <- 0.5    #right becomes 0.5
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1) # nrow values at the below
      Axis(y, side = 2) # ncol values at the left
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}


# generate matrix plot with a spectral scheme and given matrix values, wavelength versus angles
MODTRAN.matplot <- function(mat, col=NULL, pretty.wv = seq(7.5,12, 0.5),  pretty.an = seq(30,90,5),  filename="", 
                          plotKey=T, residual="F", zlim=NULL, zname=NULL, zunit=NULL){
  if (is.null(zlim)){
    zlim <- range(mat)
  }
  if(zlim[1]<0 & zlim[2] >0){
    zlim <- c(-max(abs(zlim[1]), abs(zlim[2])), max(abs(zlim[1]), abs(zlim[2])))
  }
  if (is.null(col)){
    col <- rev(brewer.pal(11,"Spectral"))
  } 
  
  if (is.null(zname)){
    if(residual=="F"){
      ztitle <- bquote(atop(NA,atop("Radiance","(W cm-2 sr-1"~mu*"m-1)")))
    }else if(residual=="S"){
      ztitle<- bquote(atop(NA,atop("Squared Residual","(W cm-2 sr-1"~mu*"m-1)"^2)))
    }else if(residual=="T"){
      ztitle <- bquote(atop(NA,atop("Residual","(W cm-2 sr-1"~mu*"m-1)")))
    }else{
      ztitle <- NULL
    }
    ztitle.line <-1.2
    ztitle.cex <- 0.8
  }else{
    if(residual=="F"){
      if(!is.null(zunit)){
        ztitle <- bquote(atop(.(zname), .(zunit)))
      }else{
        ztitle <- zname
      }
      
    }else if(residual=="S"){
      if(!is.null(zunit)){
        ztitle <- bquote(atop("Squared Residual", .(zunit)))
      }else{
        ztitle <- "Squared Residual"
      }
    }else if(residual=="T"){
      if(!is.null(zunit)){
        ztitle <- bquote(atop("Residual", .(zunit)))
      }else{
        ztitle <- "Residual"
      }
    }else{
      ztitle <- NULL
    }
    ztitle.line <-0.7
    ztitle.cex <- 0.6
  }
  wv <- as.numeric(rownames(mat))
  an <- as.numeric(colnames(mat))
  wv.todraw <- pretty.axislabels(wv, pretty.wv)
  an.todraw <- pretty.axislabels(an, pretty.an)
  # filled.contour: x is the rows of the matrix, y is the columns of the matrix so 
  # the plot is 90-degrees counterclockwise rotated to the matrix 
  if(filename==""){
    par(mar=c(2.5,2.5,1.7,0),xpd=NA,oma=c(0, 0, 0, 0))
    if(!is.null(ztitle)){
      xu.filled.contour(wv, an, mat, plot.title=title(xlab=expression(paste("Wavelength (", mu, "m)")),ylab="Angle (degrees)",line=1.3,cex.lab=0.7), color.palette=colorRampPalette(col),nlevels = 100, zlim=zlim, 
                        plot.axes={axis(1, at=wv.todraw[[1]], labels=wv.todraw[[2]], cex.axis=0.7, padj =-1.5); axis(2, at=an.todraw[[1]], label=an.todraw[[2]], cex.axis=0.7, hadj=0.4)}, 
                        key.title = title(main= ztitle, line=ztitle.line, cex.main = ztitle.cex, outer = FALSE, font.main = 1),
                        key.axes=axis(4, at=seq(zlim[1],zlim[2],length.out =5),labels=format(seq(zlim[1],zlim[2],length.out =5),scientific = TRUE,digits=2), cex.axis=0.7))
    }else{
      plot.new()
      text(x=0.45, y=0.5, labels="residual option can only be default F (False), T (True) or S (Squared)", cex=0.8)
    }
  }else{
    if(plotKey){
      jpeg(filename, width= 2000, height=1600, units = "px", quality = 100, res=300)
      par(mar=c(2.5,2.5,1.7,0),xpd=NA,oma=c(0, 0, 0, 0)) #if the mar is given here, the original xu.filled.contour will take this parameter
      # force the key axis labels to be scientific notation with two digits
      # padj: adjustment for tick label perpendicular to the reading direction
      if(!is.null(ztitle)){
        xu.filled.contour(wv, an, mat, plot.title=title(xlab=expression(paste("Wavelength (", mu, "m)")),ylab="Angle (degrees)",line=1.3,cex.lab=0.7), color.palette=colorRampPalette(col),nlevels = 100, zlim=zlim, 
                          plot.axes={axis(1, at=wv.todraw[[1]], labels=wv.todraw[[2]], cex.axis=0.7, padj =-1.5); axis(2, at=an.todraw[[1]], label=an.todraw[[2]], cex.axis=0.7, hadj=0.4)}, 
                          key.title = title(main= ztitle, line=ztitle.line, cex.main =ztitle.cex, outer = FALSE, font.main = 1),
                          key.axes=axis(4, at=seq(zlim[1],zlim[2],length.out =5),labels=format(seq(zlim[1],zlim[2],length.out =5),scientific = TRUE,digits=2), cex.axis=0.7))
      }else{
        plot.new()
        text(x=0.45, y=0.5, labels="residual option can only be default F (False), T (True) or S (Squared)", cex=0.8)
      }
      dev.off()
    }else{
      jpeg(filename, width= 600, height=600, units = "px", quality = 100, res=300)
      par(mar=c(0,0,0,0))# (bottom, left, top, right) 
      xu.filled.contour(wv, an, mat, color.palette=colorRampPalette(col),nlevels = 100, axes = FALSE, zlim=zlim, plotKey=F)
      dev.off()
    }
  }
}

# plot the original matrix and predicted matrix as well as their residuals
comparisonMatrixPlot <- function(predorigfiles.read, outputfolder, residual.opt="F", 
                                 residual.square.opt=F, zlim=NULL, zlim.residual=NULL, zname=NULL, zunit=NULL){
  predictfiles.read <- predorigfiles.read[[1]]
  originalfiles.read <- predorigfiles.read[[2]]
  filesbasename <- names(predictfiles.read)
  for (i in 1:length(predictfiles.read)){ # length(predictfiles.read)
    filebasename <- filesbasename[[i]]
    if(!is.null(zlim)){
      zlim <- zlim
    }else{
      zlim <- minmax.crosstables(list(predictfiles.read[[i]],originalfiles.read[[i]]))
    }
    MODTRAN.matplot(as.matrix(predictfiles.read[[i]]),zlim=zlim, filename=paste(outputfolder,filebasename,"_","predict", ".jpg",sep = ""),zname=zname, zunit=zunit)
    MODTRAN.matplot(as.matrix(originalfiles.read[[i]]),zlim=zlim, filename=paste(outputfolder,filebasename, "_","original", ".jpg",sep = ""),zname=zname, zunit=zunit)
    MODTRAN.matplot(as.matrix(diff.matrices(predictfiles.read[[i]],originalfiles.read[[i]], square=residual.square.opt)),
                    filename=paste(outputfolder,filebasename, "_","residual", ".jpg",sep = ""), residual=residual.opt, zlim =zlim.residual,zname=zname, zunit=zunit)
  }
}

# reform the matrix to check the residual over different variables by fixing other variables
# for example, check the residual for different days at a specific same time, viewing angle of the same reflectivity
# the mat input here is the output of residual.prep() function
matplot.residual.check <- function(mat, compDown=FALSE,  var, fix.var,  wv=wv.default, pretty.wv = seq(7.5,12, 0.5), col=NULL, 
                                   filename="",zlim=NULL, nlevels=20, plotKey=T) {
  #fix.var= c("time"=18,"reflectivity"=10,"angle"=60)
  if (is.null(zlim)){
    zlim <- range(mat)
  }
  if(zlim[1]<0 & zlim[2] >0){
    zlim <- c(-max(abs(zlim[1]), abs(zlim[2])), max(abs(zlim[1]), abs(zlim[2])))
  }
  if (is.null(col)){
    col <- rev(brewer.pal(11,"Spectral"))
  } 
  if(var=="day"){
    y <- 1:365
    ylab <- "Day"
    ylabels <- seq(1,365,30)
    title <- paste("Time: ",fix.var["time"], ":00", "; Reflectivity: ", fix.var["reflectivity"],  "%", "; Angle: ", fix.var["angle"]," degrees", sep="")
  }else if(var=="time"){
    y <- c(2,6,10,14,18,22)
    ylab <- paste("Time", "(:00)")
    ylabels <- y
    title <- paste("Day: ",fix.var["day"], "; Reflectivity: ", fix.var["reflectivity"],  "%", "; Angle: ", fix.var["angle"]," degrees", sep="")
  }else if(var == "reflectivity"){
    if(compDown) y <- c(5,10,30,50,80,100) else y <- c(0,5,10,30,50,80,100)
    ylab <- paste("Reflectivity", "(%)")
    ylabels <- y
    title <- paste("Day: ",fix.var["day"], "; Time: ", fix.var["time"], ":00", "; Angle: ", fix.var["angle"]," degrees",sep="")
  }else if(var == "angle"){
    y <- c(30, 35, 40, 50, 60, 70, 80, 90)
    ylab <- paste("Angle", "(degrees)")
    ylabels <-y
    title <- paste("Day: ",fix.var["day"], "; Time: ", fix.var["time"], ":00", "; Reflectivity: ", fix.var["reflectivity"], "%", sep="")
  }else {
    stop("'the dependent variable must be one from the day, time, reflectivity, and angle and the remaining ones should be fixed to a specfic value'")
  }
  
  wv.todraw <- pretty.axislabels(wv, pretty.wv)
  # filled.contour: x is the rows of the matrix, y is the columns of the matrix so 
  # the plot is 90-degrees counterclockwise rotated to the matrix 
  if(filename==""){
    par(mar=c(2.5,2.5,1.7,0),xpd=NA,oma=c(0, 0, 0, 0)) # mar (bottom, left, top, right)
    xu.filled.contour(wv, y, mat, plot.title=title(xlab=expression(paste("Wavelength (", mu, "m)")),ylab=ylab,line=1.5, cex.lab=0.7), color.palette=colorRampPalette(col),nlevels = nlevels, zlim=zlim, 
                        plot.axes={axis(1, at=wv.todraw[[1]], labels=wv.todraw[[2]], cex.axis=0.7, padj =-1.5); axis(2, ylabels, cex.axis=0.7, hadj=0.5)}, 
                        key.title = title(main= bquote(atop(NA,atop("Radiance Residual","(W cm-2 sr-1"~mu*"m-1)"))), line=1.2, cex.main =0.8, outer = FALSE),
                        key.axes=axis(4, at=seq(zlim[1],zlim[2],length.out =5),labels=format(seq(zlim[1],zlim[2],length.out =5),scientific = TRUE,digits=2), hadj=0.15, cex.axis=0.7))
    mtext(side=3, text=title, cex= 0.7, line=0.5, adj=0.4)
  }else{
    if(plotKey){
      jpeg(filename, width= 2000, height=1600, units = "px", quality = 100, res=300)
      par(mar=c(2.5,2.5,1.7,0),xpd=NA,oma=c(0, 0, 0, 0),mgp = c(0.7, 1, 0)) #if the mar is given here, the original xu.filled.contour will take this parameter
      # force the key axis labels to be scientific notation with two digits
      # padj: adjustment for tick label perpendicular to the reading direction
      xu.filled.contour(wv, y, mat, plot.title=title(xlab=expression(paste("Wavelength (", mu, "m)")), ylab=ylab, line=1.5, cex.lab=0.7),color.palette=colorRampPalette(col),nlevels = nlevels, zlim=zlim, 
                        plot.axes={axis(1, at=wv.todraw[[1]], labels=wv.todraw[[2]], cex.axis=0.7, padj =-1.5); axis(2, ylabels, cex.axis=0.7, hadj=0.5)}, 
                        key.title = title(main= bquote(atop(NA,atop("Radiance Residual","(W cm-2 sr-1"~mu*"m-1)"))), line=1.2, cex.main =0.8, outer = FALSE),
                        key.axes=axis(4, at=seq(zlim[1],zlim[2],length.out =5),labels=format(seq(zlim[1],zlim[2],length.out =5),scientific = TRUE,digits=2), hadj=0.15, cex.axis=0.7))
      mtext(side=3, text=title, cex= 0.7, line=0.5, adj=0.4)
      dev.off()
    }else{
      jpeg(filename, width= 600, height=600, units = "px", quality = 100, res=300)
      par(mar=c(0,0,0,0))# (bottom, left, top, right) 
      xu.filled.contour(wv, y, mat, color.palette=colorRampPalette(col),nlevels = nlevels, axes = FALSE, zlim=zlim, plotKey=F)
      dev.off()
    }
  }
}

# fix  matplot.residual
batch.matplot.residual.check <- function(input, comp ="down", check="day", plotKey=T, an=an.default){
  if(comp == 'down'){
    reflectivities <- c(5,10,30,50,80,100)
    compDown <- TRUE
  }else{
    reflectivities <- c(0,5,10,30,50,80,100)
    compDown <- FALSE
  }
  if (check=="day"){
    for (time in Localtime){
      for(reflectivity in reflectivities){
        for(angle in an){
          filename <- paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/daycheck_", time, "_", reflectivity, "_", angle, "_",comp,".jpg", sep="")
          day.check <- residual.prep(input, var="day", fix.var= c("time"=time,"reflectivity"=reflectivity,"angle"=angle))
          matplot.residual.check(as.matrix(day.check), var="day", fix.var= c("time"=time,"reflectivity"=reflectivity,"angle"=angle), filename= filename, nlevels=20, plotKey=T)
        }
      }
    }
  }
  else if(check=="time"){
    for (day in 1:365){
      for(reflectivity in reflectivities){
        for(angle in an){
          filename <- paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/timecheck_", day, "_", reflectivity, "_", angle, "_",comp, ".jpg", sep="")
          time.check <- residual.prep(input, var="time", fix.var= c("day"=day,"reflectivity"=reflectivity,"angle"=angle))
          matplot.residual.check(as.matrix(time.check), var="time", fix.var= c("day"=day,"reflectivity"=reflectivity,"angle"=angle), filename= filename, nlevels=20)
        }
      }
    }
  }
  else if(check=="reflectivity"){
    for (day in 1:365){
      for(time in Localtime){
        for(angle in an){
          filename <- paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/reflectivitycheck_", day, "_", time, "_", angle, "_", comp, ".jpg", sep="")
          reflectivity.check <- residual.prep(input, var="reflectivity", fix.var= c("day"=day,"time"=time,"angle"=angle))
          matplot.residual.check(as.matrix(reflectivity.check), compDown, var="reflectivity", fix.var= c("day"=day,"time"=time,"angle"=angle), filename= filename, nlevels=20)
        }
      }
    }
  }
  else if(check=="angle"){
    for (day in 1:365){
      for(time in Localtime){
        for(reflectivity in reflectivities){
          filename <- paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/anglecheck_", day, "_", time, "_", reflectivity, "_",comp, ".jpg", sep="")
          angle.check <- residual.prep(input, var="angle", fix.var= c("day"=day,"time"=time,"reflectivity"=reflectivity))
          matplot.residual.check(as.matrix(angle.check), var="angle",  fix.var= c("day"=day,"time"=time,"reflectivity"=reflectivity), filename= filename, nlevels=20)
        }
      }
    }
  }
  else{
    stop("Wrong variable is chosen for checking")
  }
}

# 3D array for average residual plots 
matplot.average <- function(marray, avg, col=NULL, filename="",zlim=NULL, nlevels=20, plotKey=T, plotContour=T, wv=wv.default, an=an.default, reflectivities=reflectivity.default ) {
  #fix.var= c("time"=18,"reflectivity"=10,"angle"=60)
  if (is.null(col)){
    col <- rev(brewer.pal(11,"Spectral"))
  }

  # 1 rows: wavelength; 2 columns: angles; 3 third dimension: reflectivities
  if(avg=="wavelength"){
    mat <- apply(marray,c(2,3),mean)
    mn <- apply(marray, c(2,3),min)
    mx <- apply(marray, c(2,3),max)
    x <- an
    xlab<- "Angle (degrees)"
    xlabels <- an
    y <- reflectivities
    ylab <- "Reflectivity (%)"
    ylabels <- reflectivities
  }else if(avg=="reflectivity"){
    mat <- apply(marray,c(1,2),mean)
    mn <- apply(marray, c(1,2),min)
    mx <- apply(marray, c(1,2),max)
    x <- wv
    xlab <- expression(paste("Wavelength (", mu, "m)"))
    xlabels <- seq(7.75,12,0.5)
    y <- an
    ylab <- "Angle (degrees)"
    ylabels <- an
  }else if(avg == "angle"){
    mat <- apply(marray,c(1,3),mean)
    mn <- apply(marray, c(1,3),min)
    mx <- apply(marray, c(1,3),max)
    x <- wv
    xlab <- expression(paste("Wavelength (", mu, "m)"))
    xlabels <- seq(7.75,12,0.5)
    y <- reflectivities
    ylab <- "Reflectivity (%)"
    ylabels <- reflectivities
  }else {
    stop("")
  }
  
  if (is.null(zlim)){
    zlim <- range(mat)
  }
  if(zlim[1]<0 & zlim[2] >0){
    zlim <- c(-max(abs(zlim[1]), abs(zlim[2])), max(abs(zlim[1]), abs(zlim[2])))
  }
  
  if(filename==""){
    if(plotKey){
      if(plotContour){
        par(mar=c(2.5,2.5,1.7,0),xpd=NA,oma=c(0, 0, 0, 0)) # mar (bottom, left, top, right)
        xu.filled.contour(x, y, mat, plot.title=title(xlab=xlab, ylab=ylab, line=1.5, cex.lab=0.7),
                          plot.axes={axis(1, xlabels, cex.axis=0.7, padj =-1.5); axis(2, ylabels, cex.axis=0.7, hadj=0.5);
                            contour(x,y, mn, add=T,col="blue"); contour(x,y, mx, add=T,col="green")},
                          key.title = title(main= bquote(atop(NA,atop("Radiance","(W cm-2 sr-1"~mu*"m-1)"))), line=1.2, cex.main =0.8, outer = FALSE),
                          key.axes=axis(4, at=seq(zlim[1],zlim[2],length.out =5),labels=format(seq(zlim[1],zlim[2],length.out =5),scientific = TRUE,digits=2), cex.axis=0.7),
                          color.palette=colorRampPalette(col),nlevels = nlevels, zlim=zlim, plotKey=T)
      }
      else{
        par(mar=c(2.5,2.5,1.7,0),xpd=NA,oma=c(0, 0, 0, 0)) # mar (bottom, left, top, right)
        xu.filled.contour(x, y, mat, plot.title=title(xlab=xlab, ylab=ylab, line=1.5, cex.lab=0.7),
                          plot.axes={axis(1, xlabels, cex.axis=0.7, padj =-1.5); axis(2, ylabels, cex.axis=0.7, hadj=0.5)},
                          key.title = title(main= bquote(atop(NA,atop("Radiance","(W cm-2 sr-1"~mu*"m-1)"))), line=1.2, cex.main =0.8, outer = FALSE),
                          key.axes=axis(4, at=seq(zlim[1],zlim[2],length.out =5),labels=format(seq(zlim[1],zlim[2],length.out =5),scientific = TRUE,digits=2), cex.axis=0.7),
                          color.palette=colorRampPalette(col),nlevels = nlevels, zlim=zlim, plotKey=T)
      }
    }
    else{
      if(plotContour){
        par(mar=c(2.5,2.5,1.7,0),xpd=NA,oma=c(0, 0, 0, 0))
        xu.filled.contour(x, y, mat, plot.title=title(xlab=xlab, ylab=ylab, line=1.5, cex.lab=0.7), 
                          plot.axes={axis(1, xlabels, cex.axis=0.7, padj =-1.5); axis(2, ylabels, cex.axis=0.7, hadj=0.5);
                            contour(x,y, mn, add=T,col="blue"); contour(x,y, mx, add=T,col="green")},
                          color.palette=colorRampPalette(col),nlevels = nlevels, zlim=zlim, plotKey=F)
      }
      else{
        par(mar=c(2.5,2.5,1.7,0),xpd=NA,oma=c(0, 0, 0, 0))
        xu.filled.contour(x, y, mat, plot.title=title(xlab=expression(paste("Wavelength (", mu, "m)")), ylab=ylab, line=1.5, cex.lab=0.7),
                          plot.axes={axis(1, xlabels, cex.axis=0.7, padj =-1.5); axis(2, ylabels, cex.axis=0.7, hadj=0.5)},
                          color.palette=colorRampPalette(col),nlevels = nlevels, zlim=zlim, plotKey=F)
      }
    }
  }
  else{
    if(plotKey){
      if(plotContour){
        jpeg(filename, width= 2000, height=1600, units = "px", quality = 100, res=300)
        par(mar=c(2.5,2.5,1.7,0),xpd=NA,oma=c(0, 0, 0, 0),mgp = c(0.7, 1, 0)) #if the mar is given here, the original xu.filled.contour will take this parameter
        # force the key axis labels to be scientific notation with two digits
        # padj: adjustment for tick label perpendicular to the reading direction
        xu.filled.contour(x, y, mat, plot.title=title(xlab=xlab, ylab=ylab, line=1.5, cex.lab=0.7),
                          plot.axes={axis(1, xlabels, cex.axis=0.7, padj =-1.5); axis(2, ylabels, cex.axis=0.7, hadj=0.5);
                            contour(x,y, mn, add=T,col="blue"); contour(x,y, mx, add=T,col="green")},
                          key.title = title(main= bquote(atop(NA,atop("Radiance","(W cm-2 sr-1"~mu*"m-1)"))), line=1.2, cex.main =0.8, outer = FALSE),
                          key.axes=axis(4, at=seq(zlim[1],zlim[2],length.out =5),labels=format(seq(zlim[1],zlim[2],length.out =5),scientific = TRUE,digits=2), cex.axis=0.7),
                          color.palette=colorRampPalette(col),nlevels = nlevels, zlim=zlim, plotKey=T)
        dev.off()
      }
      else{
        jpeg(filename, width= 2000, height=1600, units = "px", quality = 100, res=300)
        par(mar=c(2.5,2.5,1.7,0),xpd=NA,oma=c(0, 0, 0, 0),mgp = c(0.7, 1, 0)) #if the mar is given here, the original xu.filled.contour will take this parameter
        # force the key axis labels to be scientific notation with two digits
        # padj: adjustment for tick label perpendicular to the reading direction
        xu.filled.contour(x, y, mat, plot.title=title(xlab=xlab, ylab=ylab, line=1.5, cex.lab=0.7),
                          plot.axes={axis(1, xlabels, cex.axis=0.7, padj =-1.5); axis(2, ylabels, cex.axis=0.7, hadj=0.5)},
                          key.title = title(main= bquote(atop(NA,atop("Radiance","(W cm-2 sr-1"~mu*"m-1)"))), line=1.2, cex.main =0.8, outer = FALSE),
                          key.axes=axis(4, at=seq(zlim[1],zlim[2],length.out =5),labels=format(seq(zlim[1],zlim[2],length.out =5),scientific = TRUE,digits=2), cex.axis=0.7),
                          color.palette=colorRampPalette(col),nlevels = nlevels, zlim=zlim, plotKey=T)
        dev.off()
      }
    }
    else{
      jpeg(filename, width= 2000, height=1600, units = "px", quality = 100, res=300)
      par(mar=c(0,0,0,0))# (bottom, left, top, right) 
      xu.filled.contour(x, y, mat, color.palette=colorRampPalette(col),nlevels = nlevels, axes = FALSE, zlim=zlim, plotKey=F)
      dev.off()
    }
  }
}


# A neat function that given a vector of values, will return colors codes.
val2col <- function(x, min=NULL, max=NULL, na.rm=F, col) {
  if (is.null(min)) {min=min(x, na.rm=na.rm)}
  if (is.null(max)) {max=max(x, na.rm=na.rm)}
  
  x <- round( ( (x-min) / (max-min) ) * (length(col)-1)  + 1 )
  x[x<1] = 1
  x[x>length(col)] = length(col)
  return(col[x])
}


# plot retrieved reflectivity and each panel represents for one elevation angle
plot_retrievedreflectivity <- function(retrieved_reflectivity, true_reflectivity, filename){
  max.outlier <- max(abs(retrieved_reflectivity[which(retrieved_reflectivity < 0)])) 
  retrieved_reflectivity[which(retrieved_reflectivity < 0)] <- 0.25* retrieved_reflectivity[which(retrieved_reflectivity < 0)]/max.outlier # [-0.5, 1]
  if(length(true_reflectivity)==1){
    true_reflectivity <- cbind(wv, rep(true_reflectivity, 256))
  }
  # plot the retrieved reflectivity versus true reflectivity
  jpeg(filename, width= 2700, height=3200, units = "px", quality = 100, res=300)
  par(mfrow=c(4,2),mar=c(2,2,1,1), oma = c(2,2,0,0), mgp=c(0, 0.5, 0)) # mgp: margin line for the axis title, labels and line
  col<-boxplot.color[angles.index]
  for(i in 1:ncol(retrieved_reflectivity)){
    plot(wv, retrieved_reflectivity[,i], col=col[i], type="p",pch=19, cex=0.4, xlab = '', ylab = '',ylim=c(-0.25,1.1), yaxt="n")
    axis(2, at=c(-0.25, 0, 0.1, 0.3, 0.5, 0.8, 1), labels=c(formatC(max.outlier,format="e", digits = 1),0, 0.1, 0.3, 0.5, 0.8,1), lwd=0, lwd.tick=1) #cex.axis=0.9
    lines(true_reflectivity[,1], true_reflectivity[,2], lwd=0.4)
    legend(7.4, 1.15, bty='n', legend=paste(an[i], "°",sep=""), col=col[i], lty=1, seg.len=1, cex=1, x.intersp=0.4)
  }
  mtext(expression(paste("Wavelength (", mu, "m)")), side=1, cex = 1, outer = TRUE)
  mtext("Reflectivity [0,1]", side=2, cex = 1,  outer = TRUE)
  dev.off()
}


# plot retrieved reflectivity for different elevation angles in one panel as the boxplot
plot_retrievedreflectivity_boxplot <- function(retrieved_reflectivity, true_reflectivity=NULL, filename=NULL, wv=wv.default, an=an.default, pretty.wv = seq(7.5,12, 0.5),
                                              ylim1=c(0, 1.1), ylim2=NULL, returnRMSE=FALSE, ablineRMSE="mean"){
  # set all abnormal retrieved reflectivity less than 0 as 0
  retrieved_reflectivity <- as.matrix(retrieved_reflectivity)
  retrieved_reflectivity[which(retrieved_reflectivity < 0)] <- 0 # [-0.25, 1]
  wv.todraw <- pretty.axislabels(wv, pretty.wv) # return two variable: at, lables
  # if true refletivity is provided, plot retrieved reflectivity and its rmse 
  if(!is.null(true_reflectivity)){
    if(is.null(dim(true_reflectivity))){
      if(length(true_reflectivity)==1){
        true_reflectivity <- rep(true_reflectivity, length(wv))
      }
      else if(length(true_reflectivity)!=length(wv)){
        stop("true_reflectivity input is not at the same length as the wavelength input")
      }
    }else{
      stop("input true_reflectivity must be a vector across the spectrum")
    }
    # expand the true reflectivity to a matrix
    poly.mat <- matrix(rep(true_reflectivity, length(an)), ncol=length(an))
    rmse.median <- sqrt(apply((retrieved_reflectivity-poly.mat)^2, 1, median))
    if(returnRMSE) return(rmse.median)
    # decide the ylim for the RMSE
    if(is.null(ylim2)){
      ylim2 = c(0, max(rmse.median))
    }
    # decide which abline to draw
    abline.funcs <- c("min", "mean", "max")
    abline.legends <-  c("Minimum RMSE", "Average RMSE", "Maximum RMSE")
    if(is.na(match(ablineRMSE, abline.funcs)) ){
      if(!is.numeric(ablineRMSE)){
        stop("ablineRMSE must be numeric or one from min, mean, or max")
      }else{
        abline.at<- ablineRMSE
        abline.legend <- NULL
      }
    }else{
      abline.func <- abline.funcs[match(ablineRMSE, abline.funcs)]
      abline.at <- get(abline.func)(rmse.median)
      abline.legend <- abline.legends[match(ablineRMSE, abline.funcs)]
    }
    # plot the retrieved reflectivity versus true reflectivity
    if(is.null(filename)){      
      boxplot(t(retrieved_reflectivity), ylim=ylim1, xaxt="n", yaxt="n", outpch= 8, outcex=0.3, 
              main="Boxplot of Retrieved Reflectivity at Different Angles", 
              ylab="Reflectivity [0,1]")
      lines(1:length(wv), true_reflectivity, lwd=1, col="blue")
      axis(1, at = seq(1, length(wv), length.out=length(wv.todraw[[1]])), labels = wv.todraw[[2]], lwd=0, lwd.tick=1)
      axis(2, at = c(0, 0.1, 0.3, 0.5, 0.8, 1), labels=c(0, 0.1, 0.3, 0.5, 0.8, 1), lwd=0, lwd.tick=1)
      legend('top', bty='n', legend="True", col="blue", lty=1, seg.len=1, cex=1, x.intersp=0.4)
      # second plot for RMSE
      plot(wv, rmse.median, ylim=ylim2, col="red", type="l", xaxt="n", xlab = '', ylab = 'RMSE')
      abline(h=abline.at, col="gray")
      axis(1, at = wv.todraw[[1]], labels =wv.todraw[[2]], lwd=0, lwd.tick=1)
      axis(2, pos=0)
      legend("top", bty='n', legend=c("RMSE", abline.legend), col=c("red", "gray"), lty=1, seg.len=1, cex=1, x.intersp=0.4)
      mtext(expression(paste(bold("Wavelength ("), mu, bold("m)"))), side=1, cex = 1, outer = TRUE)
    }else{
      jpeg(filename, width= 2700, height=2400, units = "px", quality = 100, res=300)
      par(mfrow=c(2,1), mar=c(1.5,2.7,1.5,1), oma = c(1,0,0,0), mgp=c(1.4, 0.5, 0), font.lab = 2)
      boxplot(t(retrieved_reflectivity), ylim=ylim1, xaxt="n", yaxt="n", outpch= 8, outcex=0.3, 
              main="Boxplot of Retrieved Reflectivity at Different Angles", 
              ylab="Reflectivity [0,1]")
      lines(1:length(wv), true_reflectivity, lwd=1, col="blue")
      axis(1, at = seq(1, length(wv), length.out=length(wv.todraw[[1]])), labels = wv.todraw[[2]], lwd=0, lwd.tick=1)
      axis(2, at = c(0, 0.1, 0.3, 0.5, 0.8, 1), labels=c(0, 0.1, 0.3, 0.5, 0.8, 1), lwd=0, lwd.tick=1)
      legend('top', bty='n', legend="True", col="blue", lty=1, seg.len=1, cex=1, x.intersp=0.4)
      # second plot for RMSE
      plot(wv, rmse.median, ylim=ylim2, col="red", type="l", xaxt="n", xlab = '', ylab = 'RMSE')
      abline(h=abline.at, col="gray")
      axis(1, at = wv.todraw[[1]], labels =wv.todraw[[2]], lwd=0, lwd.tick=1)
      axis(2, pos=0)
      legend("top", bty='n', legend=c("RMSE", abline.legend), col=c("red", "gray"), lty=1, seg.len=1, cex=1, x.intersp=0.4)
      mtext(expression(paste(bold("Wavelength ("), mu, bold("m)"))), side=1, cex = 1, outer = TRUE)
      dev.off()
    }
  }
  else{
    if(is.null(filename)){
      boxplot(t(retrieved_reflectivity), ylim=ylim1, xaxt="n", yaxt="n", outpch= 8, outcex=0.2, 
              main="Boxplot of Retrieved Reflectivity at Different Angles", 
              ylab="Reflectivity [0,1]")
      axis(1, at = seq(1, length(wv), length.out=length(wv.todraw[[1]])), labels = wv.todraw[[2]], lwd=0, lwd.tick=1)
      axis(2, at = c(0, 0.1, 0.3, 0.5, 0.8, 1), labels=c(0, 0.1, 0.3, 0.5, 0.8, 1), lwd=0, lwd.tick=1)
    }else{
      jpeg(filename, width= 2700, height=1200, units = "px", quality = 100, res=300)
      par(mar=c(2.7,2.7,1.5,1), oma = c(0,0,0,0), mgp=c(1.6, 0.5, 0), font.lab = 1)
      boxplot(t(retrieved_reflectivity), ylim=ylim1, xaxt="n", yaxt="n", outpch= 8, outcex=0.2, 
              main="Boxplot of Retrieved Reflectivity at Different Angles", 
              ylab="Reflectivity [0,1]", xlab=expression(paste(bold("Wavelength ("), mu, bold("m)"))))
      axis(1, at = seq(1, length(wv), length.out=length(wv.todraw[[1]])), labels = wv.todraw[[2]], lwd=0, lwd.tick=1)
      axis(2, at = c(0, 0.1, 0.3, 0.5, 0.8, 1), labels=c(0, 0.1, 0.3, 0.5, 0.8, 1), lwd=0, lwd.tick=1)
      dev.off()
    }
  }
}



# true_reflectivity <- read.csv(reflectivities.exoscan[1])[,-1]
# plot retrieved reflectivity for different elevation angles in one panel as the boxplot
plot_randomerror_reflectivity_boxplot <- function(true_reflectivity, filename, wv=wv.default, an=an.default, pretty.wv = seq(7.5,12, 0.5)){
  if(is.null(dim(true_reflectivity))){
    if(length(true_reflectivity)==1){
      true_reflectivity <- rep(true_reflectivity, length(wv))
    }
    else if(length(true_reflectivity)!=length(wv)){
      stop("true_reflectivity input is not at the same length as the wavelength input")
    }
  }else{
    stop("input true_reflectivity must be a vector across the spectrum")
  }
  
  # expand the true reflectivity to a matrix
  poly.mat <- matrix(rep(true_reflectivity, length(an)), ncol=length(an))
  
  wv.error <- runif(length(wv), min = 0, max = 0.01)
  angle.error <- lapply(1:length(wv), function(x) rnorm(length(an), mean=true_reflectivity[x], sd = wv.error[x]))
  true_reflectivity.randomerror <- matrix(unlist(angle.error), nrow=length(wv), byrow=TRUE)
  true_reflectivity.randomerror[which(true_reflectivity.randomerror < 0)] <- 0
  true_reflectivity.randomerror[which(true_reflectivity.randomerror > 1)] <- 1
  
  if(max(true_reflectivity.randomerror) > 0.5){
    ylim=c(0, 1.1)
  }else if(max(true_reflectivity.randomerror > 0.3)){
    ylim=c(0, 0.5)
  }else{
    ylim=c(0, 0.3)
  }
  wv.todraw <- pretty.axislabels(wv, pretty.wv) # return two variable: at, lables
  # plot the retrieved reflectivity versus true reflectivity
  jpeg(filename, width= 4700, height=3200, units = "px", quality = 100, res=300)
  par(mar=c(2,3,2,1), oma = c(2,2,0,0), mgp=c(2, 0.5, 0))
  boxplot(t(true_reflectivity.randomerror), ylim=ylim,
          xaxt="n", yaxt="n", outpch= 8, outcex=0.5, 
          main="Boxplot for Reflectivity at Different Angles with random error added", ylab="Reflectivity [0,1]")
  lines(1:length(wv), true_reflectivity, lwd=1, col="blue")
  axis(1, at= seq(1, length(wv), length.out=length(wv.todraw[[1]])), labels = wv.todraw[[2]], lwd=0, lwd.tick=1)
  axis(2, at=c(0, 0.1, 0.3, 0.5, 0.8, 1), labels=c(0, 0.1, 0.3, 0.5, 0.8, 1), lwd=0, lwd.tick=1)
  legend('top', bty='n', legend="True", col="blue", lty=1, seg.len=1, cex=1, x.intersp=0.4)
  mtext(expression(paste("Wavelength (", mu, "m)")), side=1, cex = 1, outer = TRUE)
  dev.off()
}
