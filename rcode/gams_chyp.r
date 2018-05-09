#R code to fit GAMs of Calanus hyperboreus abundance data
#Created 2018 by Kristina Kvile (kristokv@gmail.com)

#Load R libraries:
library(maps)
library(mapproj)
library(sp)
library(mgcv)
library(plotrix)

#Load additional scripts:
source("image_scale.R") #Script for plotting image color scale
source("arcticmap.r") #Script for plotting polar-centered map
source(file="loyo.CV.r") #Script for calculating leave-one-year-out cross validation

#Load data:
load("../data/chyp_depthint.rda") #Depth-integrated abundance data (ind.m2) for C. hyperboreus stages
#This dataset contains standardized depth-integrated data (surface-bottom at least down to 100m) based on original data in chyp_orig.rda

#Columns:
#1: Sampling year
#2: Sampling month
#3: Sampling day (day of month)
#4: Sampling time (0-24, often not given)
#5: Julian day (1-365)
#6: Decimal latitude 
#7: Decimal longitude (negative values for western hemisphere)
#8: x-coordinate Azimuthal Equidistant (distance from north pole)
#9: y-coordinate Azimuthal Equidistant (distance from north pole)
#10: Distance from shelf (<500m depth) in km
#11: Bottom depth at station
#12: Upper sampling depth
#13: Lower sampling depth
#14: Total sampling depth
#15: Sampling mesh
#16: Sampling gear
#17: Dataset (see data_sources_final.xlsx)
#18: Dataset ID (see data_sources_final.xlsx)
#19: Sea surface salinity at the time and location of sampling according to PIOMAS (Pan-arctic Ice/Ocean Modeling and Assimilation System, Zhang and Rothrock, 2003)
#20: Sea ice concentration at the time and location of sampling according to PIOMAS (Pan-arctic Ice/Ocean Modeling and Assimilation System, Zhang and Rothrock, 2003)
#21-27: Depth-integrated C.hyperboreus abundance data (ind./m2) per copepodid stage C1-C5 and adult females and males (C6f/m)
#28-34: Depth-integrated abundance data (columns 21-27) corrected for varying body width/mesh size ratio 

load("../data/predframe.rda") #Dataframe for spatial predictions created based on PIOMAS model grid (Pan-arctic Ice/Ocean Modeling and Assimilation System, Zhang and Rothrock, 2003)
#Columns:
#1: Decimal latitude 
#2: Decimal longitude (negative values for western hemisphere)
#3: Bottom depth 
#4: Land mask
#5: Sea ice concentration averaged per model grid in June (1930-2016)
#6: Sea surface temperature averaged per model grid in June (1930-2016)
#7: Sea surface salinity averaged per model grid in June (1930-2016)
#8: x-coordinate Azimuthal Equidistant (distance from north pole)
#9: y-coordinate Azimuthal Equidistant (distance from north pole)
#10: Average total sampling depth in chyp_depthint.rda
#11: June
#12: Average sampling year in chyp_depthint.rda


chyp_depthint$totdep<-log(chyp_depthint$totdep) #Log-transform depth range to reduce effect of outliers
lifestages_hyp<-c(paste0("Chypc",1:5),paste0("Chypc6",c("f","m"))) #Name developmental stages


####Fit spatial GAMs####

#GAM predictor variables:
Yvar<-list(
    "s(x_meter,y_meter,k=30)",
    "s(month,bs='cc',k=5)",
    "s(totdep,k=5)",
    "s(depth,k=5)",
    "s(sss,k=5)",
    "s(ice,k=5)",
    "as.factor(year)"
    )

fams<-c("binomial","gaussian")

for(fam in fams){
  
  #Plotting outline for single additive effects:
  par(mfrow=c(length(lifestages_hyp),length(Yvar)-2),mar=c(0.5,1,0.5,1),oma=c(3,3.5,0,0))
  
  #Run models per life stage:
  for (j in lifestages_hyp){
      subdat<-chyp_depthint[!is.na(chyp_depthint[,paste0(j,"_abunm2_corr")]) & 
                           !is.na(chyp_depthint$ice) &
                           !is.na(chyp_depthint$sss),]
      
      #Convert data to 0/1 for binomial or log-transformed positive values for gaussian:
      if(fam=="binomial"){
        subdat[,paste0(j,"_abunm2_corr")][subdat[,paste0(j,"_abunm2_corr")]>0]<-1 
      }
      if(fam=="gaussian"){
        subdat<-subdat[subdat[,paste0(j,"_abunm2_corr")]>0,]
        subdat[,paste0(j,"_abunm2_corr")]<-log(subdat[,paste0(j,"_abunm2_corr")])
      }
      
      #Fit model using abundance corrected for different mesh size:
      Xvar<-paste0(j,"_abunm2_corr~")
      mod_full<-gam(formula(paste(Xvar, paste(Yvar, collapse="+"))),data=subdat,family=fam)
      
      #Calculate deviance explained and loyoCV for full model (factor effect of year not included for loyoCV):
      #Note: loyoCV takes a long time to run
      devframe<-array(dim=c(1,1+length(Yvar),2),dimnames=list("Full",c("Tot",Yvar),c("Dev.","CV")))
      devframe[1,"Tot","Dev."]<-100*round(as.numeric(summary(mod_full)["dev.expl"]),3)
      devframe[1,"Tot","CV"]<-loyo.CV(form=as.formula(paste(Xvar,paste(Yvar[-length(Yvar)],collapse="+"))),
            type="gam",data=subdat,grouping.variable="year",family=fam,output="mean")
      #Find partial deviance explained and loyo cv by removing single terms at a time:
      for (i in (c(1:length(Yvar)))){
        devframe[1,unlist(Yvar[i]),"Dev."]<-round(devframe[1,"Tot","Dev."]-100*as.numeric(summary(gam(formula(paste(Xvar, paste(Yvar[-i], collapse="+"))),
                                         data=subdat,family=fam))["dev.expl"]),3)
        if(i<length(Yvar)){
        devframe[1,unlist(Yvar[i]),"CV"]<-loyo.CV(form=as.formula(paste(Xvar,paste(Yvar[-c(i,length(Yvar))],collapse="+"))),
              type="gam",data=subdat,grouping.variable="year",family=fam,output="mean")
        }
      }
    #Get predicted abundances in space with other Yvars at mean value
    pred=predict.gam(mod_full,type="response",newdata=predframe,se=F)
    predframe[,paste("pred_mod",fam,j,sep="_")]<-pred
  
    #Plot single additive model terms
    varnames<-c("","month","totdep","depth","sss","ice")
    labnames<-c("","Month","log(Depth interval)","Bottom depth","Salinity","Ice concentration")
    for(i in 2:length(Yvar[-length(Yvar)])){
      plot.gam(mod_full,select=i,xlab="",ylab="",rug=TRUE,shade=TRUE,scale=0,main="",
               axes=FALSE,frame=TRUE,xlim=range(chyp_depthint[,varnames[i]],na.rm=TRUE))
      Axis(side=1,labels = FALSE); Axis(side=2)
      if(j==lifestages_hyp[7]){Axis(side=1)}
      abline(h=0,lty=2) #0-effect isoline
      #Mark if partial deviance>1% with thicker panel width:
      devpart<-devframe[1,paste0(Yvar[i]),"Dev."]
      if(devpart>=1){ box(lwd=3)}
      #Last labeling:
      if(j==lifestages_hyp[7]){mtext(side=1,labnames[i],line=2,cex=0.75)}
      if(i==2){mtext(side=2,gsub("hypc", "",j, fixed = TRUE),line=2,cex=0.75,outer=FALSE)}
    }
    } 
   mtext(side=2,"Additive effect on response variable",line=2.5,cex=0.75,outer=TRUE)
}
 
####Plot predicted abundances in space####

fams<-c("binomial","gaussian","combined") #"combined" used to plot product of gaussian and binomial

for(fam in fams){
  #Merge predictions for all stages into 1 column to define common color scale:
  pred_all<-data.frame(matrix(ncol=4,nrow=0));colnames(pred_all)<-c("lon","lat","stage","pred")
  
  for (j in c(lifestages_hyp)){
      if(fam!="combined"){
        pred_sub<-data.frame(lon=predframe$lon,lat=predframe$lat,stage=j,
                             pred=predframe[,paste("pred_mod",fam,j,sep="_")]);
        pred_all<-rbind(pred_all,pred_sub);
        }else{
        #Calculate the product of binomial and gaussian model predictions:
        pred_gaus<-predframe[,paste("pred_mod_gaussian",j,sep="_")]
        pred_binom<-predframe[,paste("pred_mod_binomial",j,sep="_")]
        pred_sub<-data.frame(lon=predframe$lon,lat=predframe$lat,stage=j,
                             pred=pred_gaus*pred_binom)
        pred_all<-rbind(pred_all,pred_sub)
      }
  }
  
  #Define color scale:
  rbPal <- colorRampPalette(c("#253494","#2c7fb8","#41b6c4","#a1dab4","#ffffcc"))(10)
  pred_all$col<- rbPal[as.numeric(cut(pred_all$pred,breaks = 10))]
  
  #Plot predictions per stage:
  layout(matrix(c(1:8,5:7,9),nrow=3,byrow=TRUE),heights=c(0.5,0.25,0.25))
  par(mar=rep(0,4),oma=c(0,1,1,1),ps=10,mgp=c(3,0.5,0))
  for (j in c(lifestages_hyp)){
    pred_sub<-pred_all[pred_all$stage==j,]
    #Plot color contour:
    with(mapproject(pred_sub$lon,pred_sub$lat, proj_str, parameters=proj_params, orientation=c(90, 0, rotation)),
         plot(x,y,col=pred_sub$col,cex=0.25,axes=FALSE,xlab="",ylab=""))
    #Add map:
    #arcticmap(add=T) 
    #Add datapoints to plot:
    subdat<-chyp_depthint[!is.na(chyp_depthint[,paste0(j,"_abunm2_corr")]) & 
                          !is.na(chyp_depthint$ice) &
                          !is.na(chyp_depthint$sss),]
    #Plot zeroes as crosses (not for Gaussian):
    if(fam!="gaussian"){
      with(mapproject(subdat$lon[subdat[,paste0(j,"_abunm2")]==0],subdat$lat[subdat[,paste0(j,"_abunm2")]==0],
                      proj_str,parameters=proj_params,orientation=c(90, 0, rotation)),
           points(x, y,pch=4))}
    #Plot positive values as circles (binomial):
    if(fam=="binomial"){
      subdat[,paste0(j,"_abunm2_corr")][subdat[,paste0(j,"_abunm2_corr")]>0]<-1;
      with(mapproject(subdat$lon[subdat[,paste0(j,"_abunm2")]>0],subdat$lat[subdat[,paste0(j,"_abunm2")]>0],
                      proj_str,parameters=proj_params,orientation=c(90, 0, rotation)),
           points(x, y,pch=21,bg="grey",cex=1))
    }
    #Plot positive values with circle size depending on observation (not binomial):
    if(fam!="binomial"){
      subdat<-subdat[subdat[,paste0(j,"_abunm2_corr")]>0,];
      subdat[,paste0(j,"_abunm2_corr")]<-log(subdat[,paste0(j,"_abunm2_corr")]);
      with(mapproject(subdat$lon,subdat$lat,
          proj_str,parameters=proj_params,orientation=c(90, 0, rotation)),
          points(x, y,pch=21,bg="grey",cex=0.25*subdat[,paste0(j,"_abunm2_corr")]))
    }
    #Developmental stage label:
    mtext(side=3,gsub("hypc","",j),outer=FALSE,cex=1.2,adj=0.6,line=-.75)
  }
  
  #Bubble plot legend:
  par(mar=c(2,7,1,2))
  plot(1, type="n",axes=FALSE,xlim=c(0,4),ylim=c(0,10),xlab="",ylab="")
  if(fam!="binomial"){
  vals<-as.vector(as.matrix(chyp_depthint[,c(28:34)]))
  vals<-log(vals[!is.na(vals) & vals>0])
  vals<-vals[vals>=0] #Remove a few negative values
  points(x=rep(2.5,5),y=seq(1,by=2,length=5),pch=21,bg="dark grey",
         cex=0.25*round(seq(min(vals),max(vals),length=5),1))
    if(fam=="combined"){
      points(2.5,1,pch=4);
      text(rep(0.5,5),seq(1,by=2,length=5),format(round(seq(min(vals),max(vals),length=5),1),trim=T))}
    #Do not include zeroes for gaussian:
    if(fam=="gaussian"){
      text(rep(0.5,5),seq(2.5,by=2,length=4),format(round(seq(min(vals),max(vals),length=5),1)[-1],trim=T))}
  mtext(side=2,"Observations",outer=FALSE,line=3.5)
  mtext(side=2,expression(paste("log"["e"],"ind.",m^-2)),outer=FALSE,line=1.75)
  }
  if(fam=="binomial"){
    points(x=2.5,y=2.5,pch=21,bg="dark grey",cex=1);
    points(2,1,pch=4);
    text(0.5,1,0)
    text(0.5,2.5,1)
    }
    
  #Image scale:
  Col_sorted<-pred_all$col[order(pred_all$pred)]
  par(mar=c(2,7,1,3))
  image.scale(sort(unique(pred_all$pred)),col=unique(Col_sorted),horiz=FALSE) 
  if(fam=="binomial"){
    mtext(side=2,"Predcited presence",outer=FALSE,line=3.5)
    mtext(side=2,"(probability)",outer=FALSE,line=1.75)
  }else{
    mtext(side=2,"Predictions",outer=FALSE,line=3.5)
    mtext(side=2,expression(paste("log"["e"],"ind.",m^-2)),outer=FALSE,line=1.75)
  }
  
}


