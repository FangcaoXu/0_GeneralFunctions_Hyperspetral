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


library(colorRamps)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(reshape2)
library(stringr)
library(scales)
library(tidyverse)
library(viridis) 
library(RColorBrewer)

# Section 1: Read the data -------------------------------------------------------------------------------------
result <- head(read.csv('/amethyst/s0/fbx5002/geolab_storage_V3/docs/DARPA/Results/Table/Radiance15cm-1.csv', skip=6, header = FALSE, stringsAsFactors = FALSE),-1)

# refer to the exact meaning of column names in NotationDoc in the DARPA project
colnames(result) <- c("freq", "path_trans", "path_emiss","path_thermal_scat","surface_emiss", "path_multi_scat",
                      "sing_scat","grnd_rflt","drct_rflt","total_rad","irrad_ref","irrad_obs","-natlog_path_trans","direct_emiss",
                      "ToA_irrad", "bbody_temp")
for(i in 1:ncol(result)){
  result[,i] <- as.numeric(result[,i])
}
result$wavelength <- 1e4/result$freq
result <- result[which(result$wavelength<=13 & result$wavelength>=7),-1]
result <- result[order(result$wavelength),]
result <- cbind(result[,16], result[,-16])
colnames(result)[1] <- "wavelength"
rownames(result) <- 1: nrow(result)

# For"W" wavenumber, the radiance unit is W cm-2 sr-1 /cm-1 and the irradiance unit is W cm-2/cm-1 ; 
# for "M", micron, the radiance unit is W cm-2 sr-1/um and the irradiance unit is W cm-2/um;
# for "N", nanometer, the radiance unit is uW cm-2 sr-1/nm and the irradiance unit is uW cm-2 sr-1/nm
# Set the grid layouts, palettes and themes for plots
layoutgrid<- matrix(c(1,1,2,3,1,1,4,5,6,6,7,7),nrow=3,ncol=4,byrow=TRUE)
layoutgrid1<- matrix(c(2,3,4,5,6,1,7,1),nrow=4,ncol=2,byrow=TRUE)
layoutgrid2<- matrix(c(1,1,2,3,4,5,6,7),nrow=4,ncol=2,byrow=TRUE)
palette1 <-brewer.pal(8, "Set1")[c(1:5,7,8)]
pie(rep(1,7), col =palette1, border = NA,label=NA)
palette2 <- c("#00FF3F","#00FF7F","#00FFBF","#00FFFF","#00BFFF","#007FFF","#003FFF","#0000FF")

# Section 2: Use plot() to plot the data -------------------------------------------------
# colfunc <- colorRampPalette(c("red","orange"))  colfunc(6)
# rgb(t(col2rgb("orange")), maxColorValue=255)


### Plot 1:
plot1<- function(plottable,layout){
  layout(layout)
  plot(plottable$wavelength, plottable$total_rad, xlab=expression(paste("Wavelength (", mu, "m)")), ylab=expression(paste("W cm-2 sr-1/", mu, "m")),
     main = "Total radiance at sensor", type='l',col="red3")
  plot(plottable$wavelength, plottable$path_emiss, xlab=expression(paste("Wavelength (", mu, "m)")), ylab=expression(paste("W cm-2 sr-1/", mu, "m")),
       main = "Path thermal emission radiance", type='l',col=palette1[1], ylim = c(0,max(max(plottable$path_emiss),max(plottable$surface_emiss))))
  plot(plottable$wavelength, plottable$surface_emiss, xlab=expression(paste("Wavelength (", mu, "m)")), ylab=expression(paste("W cm-2 sr-1/", mu, "m")),
     main="Surface thermal emission radiance", type='l',col=palette1[2], ylim = c(0,max(max(plottable$path_emiss),max(plottable$surface_emiss))))
  plot(plottable$wavelength, plottable$grnd_rflt, xlab=expression(paste("Wavelength (", mu, "m)")), ylab=expression(paste("W cm-2 sr-1/", mu, "m")),
     main="Ground reflected radiance", type='l',col=palette1[3])
  plot(plottable$wavelength, plottable$path_thermal_scat, xlab=expression(paste("Wavelength (", mu, "m)")), ylab=expression(paste("W cm-2 sr-1/", mu, "m")),
     main="Path thermal scattering radiance", type='l',col= palette1[4])
  plot(plottable$wavelength, plottable$sing_scat, xlab=expression(paste("Wavelength (", mu, "m)")), ylab=expression(paste("W cm-2 sr-1/", mu, "m")),
     main="Single scattering solar radiance", type='l',col=palette1[5], ylim = c(0,max(max(plottable$sing_scat),max(plottable$path_multi_scat))))
  plot(plottable$wavelength, plottable$path_multi_scat, xlab=expression(paste("Wavelength (", mu, "m)")), ylab=expression(paste("W cm-2 sr-1/", mu, "m")),
     main="Path multiple scattering solar radiance", type='l',col=palette1[6], ylim = c(0,max(max(plottable$sing_scat),max(plottable$path_multi_scat))))
}
plot1(result,layoutgrid)

# Section 3: Use ggplot to plot the data -------------------------------------------------
### Plot 2: Highly customized setting for each plot
theme1 <- theme(axis.title=element_blank(), axis.text = element_text(size=4), plot.title = element_text(hjust = 0.5,face = "bold",size=4, margin=margin(0,0,0,0)), 
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(1,5,0,0), "points")) # top,right,bottom,left

ggplot1 <- function(plottable, layout, theme, imgwidth, imgheight, filename, legend.pos, log=FALSE){
  if(log==TRUE){
    plottable[,c(2:8,10)] = log(plottable[c(2:8,10)] ,10)
  }
  # plot all data in one grid together with the total radiance
  plot.data <- melt(plottable, id.vars ="wavelength", 
                    measure.vars = c( "total_rad","path_emiss","surface_emiss", "grnd_rflt","path_thermal_scat","sing_scat","path_multi_scat"))
  ylimits1<- c(min(c(plottable$sing_scat,plottable$path_multi_scat)),max(c(plottable$sing_scat,plottable$path_multi_scat)))
  ylimits2<- c(min(c(plottable$path_emiss,plottable$surface_emiss)),max(c(plottable$path_emiss,plottable$surface_emiss)))
  
  g4 <- ggplot(plottable,aes(x=wavelength, y=grnd_rflt)) + geom_line(color=palette1[4], size=0.2) + theme_bw() + theme + theme(legend.position = "none") + scale_y_continuous(breaks=pretty_breaks(3))+
    annotate("text", x= 7, y=0.98*max(plottable$grnd_rflt), label = "c", size=1, fontface="bold") + labs(title= "Ground reflected radiance")
  g5 <- ggplot(plottable,aes(x=wavelength, y=path_thermal_scat)) + geom_line(color=palette1[5],size =0.2) + theme_bw() + theme + theme(legend.position = "none") + scale_y_continuous(breaks=pretty_breaks(3))+
    annotate("text", x= 7, y=0.98*max(plottable$path_thermal_scat), label = "d", size=1, fontface="bold") + labs(title= "Path thermal scattering radiance")
  if(log==TRUE){
    g1 <- ggplot(plot.data,aes(x=wavelength, y=value, color = variable)) + geom_line(size=0.2) + theme_bw() + theme + labs(title= "Total radiance at sensor and it's components")+
      theme(legend.position =legend.pos, legend.title = element_blank(), legend.text = element_text(size=3, face="bold"), 
            legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.7, 'lines'), legend.background = element_rect(fill="transparent")) + 
      guides(color=guide_legend(keyheight=0))+ scale_color_manual(values = palette1, labels = c("Sum of Figs a~f", "a","b","c","d","e","f"))
      
    g2 <- ggplot(plottable,aes(x=wavelength, y=path_emiss)) + geom_line(color=palette1[2], size=0.2) + theme_bw() + theme +  theme(legend.position = "none")+
      scale_y_continuous(breaks=pretty_breaks()) + annotate("text", x= 7, y=0.98*max(plottable$path_emiss), label = "a", size=1, fontface="bold")+
      labs(title= "Path thermal emission radiance")
    g3 <- ggplot(plottable,aes(x=wavelength, y=surface_emiss)) + geom_line(color=palette1[3], size=0.2) + theme_bw() + theme +  theme(legend.position = "none")+
      scale_y_continuous(breaks=pretty_breaks()) + annotate("text", x= 7, y=0.98*max(plottable$surface_emiss), label = "b", size=1, fontface="bold") + 
      labs(title= "Surface thermal emission radiance")
    g6 <- ggplot(plottable,aes(x=wavelength, y=sing_scat)) + geom_line(color=palette1[6], size=0.2) +theme_bw() + theme  + theme(legend.position = "none")+
      scale_y_continuous(limits=ylimits1,breaks=pretty_breaks(3)) +  annotate("text", x= 7, y=1.02*ylimits1[2], label = "e", size=1, fontface="bold")+ 
      labs(title= "Single scattering solar radiance")
    g7 <- ggplot(plottable,aes(x=wavelength, y=path_multi_scat)) + geom_line(color=palette1[7], size=0.2) + theme_bw() + theme + theme(legend.position = "none") +
      scale_y_continuous(limits=ylimits1, breaks=pretty_breaks(3)) + annotate("text", x= 7, y=1.02*ylimits1[2], label = "f", size=1, fontface="bold")+ 
      labs(title= "Multiple scattering solar radiance")
    }else{
      g1 <- ggplot(plot.data,aes(x=wavelength, y=value, color = variable)) + geom_line(size=0.2) + theme_bw() + theme + labs(title= "Total radiance at sensor and it's components")+
        theme(legend.position = legend.pos, legend.title = element_blank(), legend.text = element_text(size=3, face="bold"), 
              legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.7, 'lines'), legend.background = element_rect(fill="transparent")) + 
        guides(color=guide_legend(keyheight=0))+ scale_color_manual(values = palette1, labels = c("Sum of Figs a~f", "a","b","c","d","e","f"))
      g2 <- ggplot(plottable,aes(x=wavelength, y=path_emiss)) + geom_line(color=palette1[2], size=0.2) + theme_bw() + theme + theme(legend.position = "none") +
        scale_y_continuous(breaks=pretty_breaks(), limits=ylimits2) + annotate("text", x= 7, y=0.98*ylimits2[2], label = "a", size=1, fontface="bold")+
        labs(title= "Path thermal emission radiance")
      g3 <- ggplot(plottable,aes(x=wavelength, y=surface_emiss)) + geom_line(color=palette1[3], size=0.2) + theme_bw() + theme + theme(legend.position = "none") +
        scale_y_continuous(breaks=pretty_breaks(),limits=ylimits2) + annotate("text", x= 7, y=0.98*ylimits2[2], label = "b", size=1, fontface="bold") + 
        labs(title= "Surface thermal emission radiance")
      g6 <- ggplot(plottable,aes(x=wavelength, y=sing_scat)) + geom_line(color=palette1[6], size=0.2) +theme_bw() + theme  + theme(legend.position = "none") +
        scale_y_continuous(limits=ylimits1,breaks=pretty_breaks(3)) +  annotate("text", x= 7, y=0.98*ylimits1[2], label = "e", size=1, fontface="bold")+ 
        labs(title= "Single scattering solar radiance")
      g7 <- ggplot(plottable,aes(x=wavelength, y=path_multi_scat)) + geom_line(color=palette1[7], size=0.2) + theme_bw() + theme + theme(legend.position = "none") +
        scale_y_continuous(limits=ylimits1, breaks=pretty_breaks(3)) + annotate("text", x= 7, y=0.98*ylimits1[2], label = "f", size=1, fontface="bold")+ 
        labs(title= "Multiple scattering solar radiance")
    }
  jpeg(filename, width= imgwidth, height=imgheight, units = "px", quality = 100, res=300)
  grid.arrange(g1,g2,g3,g4,g5,g6,g7, ncol = 3,layout_matrix = layout, 
               top =textGrob("Total radiance at sensor and each component", gp = gpar(fontface = "bold", cex = 0.5)),
               left =textGrob(expression(paste("(W cm-2 sr-1/", mu, "m)")), rot = 90, gp = gpar(fontface = "bold", cex = 0.4)),
               bottom =textGrob(expression(paste("Wavelength (", mu, "m)")), gp = gpar(fontface = "bold", cex = 0.4))) 
  dev.off()
}
ggplot1(result, layoutgrid2, theme1, 900, 1200, "ggplot1.jpg",c(0.1,0.67))
ggplot1(result, layoutgrid1, theme1, 1200, 900, "ggplot1log.jpg", c(0.85,0.2), log=TRUE)


### Plot 3: Less customized but simple to draw all plots
# backticks `  is to delimit variable names in formulae
plot.data <- melt(result, measure.vars = c("path_trans", "path_emiss","surface_emiss", "grnd_rflt","path_thermal_scat","sing_scat",
                                           "path_multi_scat","drct_rflt","total_rad","irrad_ref","irrad_obs",
                                           "-natlog_path_trans","direct_emiss","ToA_irrad", "bbody_temp"))
theme2 <- theme(legend.position = "none", axis.title.y=element_blank(), axis.title.x=element_text(face = "bold",size=12), strip.text = element_text(size = 10, face ="bold"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background= element_blank())
label <- c(`total_rad`="Total radiance at sensor",`path_emiss` = "Path thermal emission radiance", `surface_emiss`="Surface thermal emission radiance",
           `grnd_rflt`="Ground reflected radiance", `path_thermal_scat`="Path thermal scattering radiance", `sing_scat`="Single scattering solar radiance", 
           `path_multi_scat`="Path multiple scattering solar radiance", `path_trans` = "Direct transmittance", `drct_rflt`="Direct reflected solar radiance", 
           `irrad_ref`="Reflected solar irradiance", `irrad_obs`="Solar irradiance at observer", `-natlog_path_trans`="-log(Path transmittance)", 
           `direct_emiss`="Ground directional emissivity", `ToA_irrad`="Top of Atmosphere solar irradiance", `bbody_temp`="Blackbody brightness temperature \nto emit the total radiance")
ggplot2 <- function(plottable, theme, label, imgwidth, imgheight, filename){
  grid1<- ggplot(plottable[which(plottable$variable == "total_rad"),],aes(x=wavelength, y=value)) + geom_line(color="red3", size=0.7) + labs(x=expression(paste("Wavelength (", mu, "m)")))+ 
    theme_bw() + theme + facet_wrap(variable~.,scales = "free_y", labeller = as_labeller(label), nrow = 2) +
    annotate("text", x= -Inf, y=Inf, hjust = -0.2, vjust = 4, label = "Sum of Figs a~f", size=3, fontface="bold")
  grid2 <- ggplot(plottable[which(plottable$variable %in% c("path_emiss","surface_emiss", "grnd_rflt", "path_thermal_scat","sing_scat", "path_multi_scat")),], aes(x=wavelength, y=value)) + 
    geom_line(aes(color=variable)) + labs(x=expression(paste("Wavelength (", mu, "m)")))+ theme_bw()+ theme + facet_wrap(variable~., scales = "free_y",labeller = as_labeller(label), nrow = 2) + 
    annotate("text", x= -Inf, y=Inf, hjust = -0.5, vjust = 2, label = letters[1:6], size=3, fontface="bold") + scale_color_manual(values=palette1[2:7])
  grid3 <- ggplot(plottable[-which(plottable$variable %in% c("total_rad","path_emiss","surface_emiss", "grnd_rflt", "path_thermal_scat","sing_scat", "path_multi_scat")),],aes(x=wavelength, y=value)) + 
    geom_line(aes(color=variable)) + labs(x=expression(paste("Wavelength (", mu, "m)")))+ theme_bw()+ theme + facet_wrap(variable~.,scales = "free_y", labeller = as_labeller(label), nrow = 2) + 
    scale_color_manual(values=rev(palette2))
  jpeg(filename, width =imgwidth, height =imgheight, units = "px", quality = 100, res=300)
  grid.arrange(grid1, grid2, grid3, ncol = 2, widths = c(3, 7),layout_matrix = cbind(c(1,3),c(2,3)))   # gtable
  dev.off()
}
ggplot2(plot.data,theme2, label, 6000,4000,"ggplot2")
