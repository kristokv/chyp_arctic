#R code to plot seasonal variation in Calanus hyperboreus C3 concentrations against growth season
#Created 2018 by Kristina Kvile (kristokv@gmail.com)

#Load R libraries:
library(maps)
library(mapproj)

#Load additional scripts:
source("arcticmap.r") #Script for plotting polar-centered map

#Load data:
load("../data/chyp_orig.rda") #Depth-specific abundance data (ind.m2 or ind.m3) for C. hyperboreus stages
#This dataset contains original data extracted from various sources. Some of the data are depth-integrated (surface-bottom),
#other data are for specific depth layers, with multiple layers per station.

#Columns:
#1: Sampling year
#2: Sampling month
#3: Sampling day (day of month)
#4: Sampling time (0-24, often not given)
#5: Julian day (1-365)
#6: Decimal latitude 
#7: Decimal longitude (negative values for western hemisphere)
#8: Distance from shelf (<500m depth) in km
#9: Bottom depth at station
#10: Upper sampling depth
#11: Lower sampling depth
#12: Sampling mesh
#13: Sampling gear
#14: Dataset (see data_sources_final.xlsx)
#15: Dataset ID (see data_sources_final.xlsx)
#16-22: C.hyperboreus abundance data (ind./m3) per copepodid stage C1-C5 and adult females and males (C6f/m)
#23-30: C.hyperboreus abundance data (ind./m2) per copepodid stage C1-C5 and adult females and males (C6f/m)

load("../data/phyto_daymean.rda") #Daily average phytoplankton concentrations/phytoplankton+microzooplankton in upper 60 m 
#for 3 different regions in 1987 and 2012 from BIOMAS: coupled pan-arctic Biology/Ice/Ocean Modeling and Assimilation System (Zhang et al. 2015)

lifestages_hyp<-c(paste0("Chypc",1:5),paste0("Chypc6",c("f","m"))) #Name developmental stages
j<-lifestages_hyp[3] #Select Chyp stage C3

#Region divisions:
latmin<-c(85,70,70)
latmax<-c(90,85,90)
shelfmin<-c(200,50,0)
shelfmax<-c(1000,1000,50)

#Plot:
legends<-c("(a)","(b)","(c)")
cbPalette <- c("magenta","sky blue","yellow")
par(mfrow=c(3,2),mar=c(0,1,0,0.5),oma=c(2.2,0,.5,2),ps=8,mgp=c(3,0.5,0))

#Loop through and plot for the 3 regions:
for(i in 1:3){ 
  #Extract C3 data (ind.m3) for region i:
  subdat<-chyp_orig[!is.na(chyp_orig$dist_shelf) & chyp_orig$dist_shelf<shelfmax[i] &
                    chyp_orig$dist_shelf>=shelfmin[i] & 
                    chyp_orig$lat>latmin[i] & 
                    chyp_orig$lat<=latmax[i] & 
                    !is.na(chyp_orig[,paste0(j,"_abunm3")]),] #& 
  #Divide into deep and shallow samples:
  subdat_shallow<-subdat[subdat$lowdep<=200,]
  subdat_deep<-subdat[ subdat$updep>=200,]
  #Map the distribution of data:
  arcticmap(addmap=FALSE,col="grey",map_labels = FALSE,map_grid = TRUE)
  with(mapproject(subdat$lon,subdat$lat,proj_str,parameters=proj_params,orientation=c(90, 0, rotation)),
       points(x, y,cex=0.25,pch=20,col="black"))
  legend("topleft",legends[i],bty="n")
  
  #Plot phytoplankton concentrations:
  phyto_daymean_1987<-phyto_daymean[,paste0("phyto_",1987,"_lat>",latmin[i],".shelf>",shelfmin[i])]
  phyto_daymean_2012<-phyto_daymean[,paste0("phyto_",2012,"_lat>",latmin[i],".shelf>",shelfmin[i])]
  plot(phyto_daymean_1987,type="l",col="white", axes=F, xlab=NA, ylab=NA,ylim=c(0,75000))
  polygon(x=c(1:365,365:1),y=c(phyto_daymean_2012,rep(0,365)),col=adjustcolor("#009E73",alpha.f=0.2),border=NA)
  polygon(x=c(1:365,365:1),y=c(phyto_daymean_1987,rep(0,365)),col=adjustcolor("#009E73",alpha.f=0.4),border=NA)
  axis(side = 4)
  
  #Plot C3 concentrations:
  par(new = T)
  plot(subdat_shallow$julday,log(subdat_shallow[,paste0(j,"_abunm3")]+1),pch=20,col="#D55E00",cex=0.25,
       xlim=c(1,365),ylim=c(0,7),axes=FALSE);box()
  Axis(side=2,at=seq(0,6,by=2),labels=seq(0,6,by=2))
  points(subdat_deep$julday,log(subdat_deep[,paste0(j,"_abunm3")]+1),pch=20,col="#0072B2",cex=0.25)
  
  #Add legends:
  if(i==1){legend("topleft",legend=c("<200m",">200m"),title="C3 concentration",bty="n",
                  pch=20,col=c("#D55E00","#0072B2"))}
  if(i==1){legend("topright",legend=c(1987,2012),title="Growth season",bty="n",
                  pch=20,col=c(adjustcolor("#009E73",alpha.f=0.4),adjustcolor("#009E73",alpha.f=0.2)))}
  if(i==3){Axis(side=1);mtext(side=1,"Day-of-year",line=1.1,cex=0.8)}
  if(i==2){mtext(side=2,expression(paste("C3 concentration (log"["e"],"ind.m"^"-3","+1)")),line=1,adj=1.25,cex=0.8)}
  }
mtext(side=4,expression(paste("Phytoplankton concentration (μgCm"^"-3",")")),outer=TRUE,line=0.75,cex=0.8)



