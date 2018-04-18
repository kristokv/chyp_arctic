#R code to estimate Calanus hyperboreus energy requirements to develop from N3-C3
#Created 2018 by Kristina Kvile (kristokv@gmail.com)

#Load data:
pred_c<-get(load(file="../data/carb_weight_stages.rda")) #Carbon weights per dev. stage 
#Age 2:6 = naupliar stage N2:N6, 7:9 = copepodite stages C1:C3
#Data for stages N3 and N4 from Jung-Madsen et al. 2013: Early development of Calanus hyperboreus nauplii: Response to a changing ocean. Limnol Oceanogr 58(6):2109–2121.
#Data for stages N6-C3 from Ashjian et al. 2003: Annual cycle in abundance, distribution, and size in relation to hydrography of important copepod species in the western Arctic Ocean. Deep Sea Res Part I Oceanogr Res Pap 50(10–11):1235–1261.
#We assumed weight of N2=N3 unfed and N5=average of N4 and N6

load("../data/phyto_daymean.rda") #Daily average phytoplankton concentrations/phytoplankton+microzooplankton in upper 60 m 
#for 3 different regions in 1987 and 2012 from BIOMAS: coupled pan-arctic Biology/Ice/Ocean Modeling and Assimilation System (Zhang et al. 2015)
load("../data/filtrations.rda") #Observed filtration rates per stage from Campbell et al. 2009: Mesozooplankton prey preference and grazing impact in the western Arctic Ocean. Deep Res Part II Top Stud Oceanogr 56(17):1274–1289.

####Setup analysis####

#1. Calculate metabolism based on carbon weight
#Function from Maps et al 2014: A metabolic approach to dormancy in pelagic copepods helps explaining inter- and intra-specific variability in life-history strategies. J Plankton Res 36(1):18–30.
#Standardized to 0dg for active phase individuals based on meta-analysis

#Per individual per second:
meta_µg_C_s<-7*(10^-7)*pred_c$C_µg^(0.76)
#Per individual per day:
pred_c$meta_µg_C_day<-meta_µg_C_s*60*60*24

#2. Calculate relative development time per stage based on Behlradek functions
#Based on Ji R, et al. 2012: Life history and biogeography of Calanus copepods in the Arctic Ocean: An individual-based modeling study. Prog Oceanogr 96(1):40–56.
#D=a(T+α)^β
coefs<-c(1461,3485,1907,1799,2113,2427,2856,3588) #Stage-specific coefficients for a
alfa<-13.66 #α
b<-(-2.05) #β
#Development at 0 degrees:
pred_c$dev_0_dg<-NA
pred_c$dev_0_dg[1:7]<-coefs[1:7]*(0+alfa)^b
pred_c$dev_0_dg[8]<-(coefs[8]*(0+alfa)^b)/2 #We consider that C3 only develops halfway before diapause
#Total and relative development time N3-C3 (relative assumed constant):
tot_dev_0_dg<-sum(pred_c$dev_0_dg[2:8])
pred_c$rel_dev_stages<-NA
pred_c$rel_dev_stages[2:8]<-pred_c$dev_0_dg[2:8]/tot_dev_0_dg


#3. Calculate development time in days based on duration of growth season in northern Arctic for 1987 and 2012:

#Region divisions:
latmin<-c(85,70,70)
latmax<-c(90,85,90)
shelfmin<-c(200,50,0)
shelfmax<-c(1000,1000,50)
i<-1 #Northernmost region

for(y in c(1987,2012)){
  #Get daily phytoplankton concentrations:
  phyto_daymean_y<-phyto_daymean[,paste0("phyto_",y,"_lat>",latmin[i],".shelf>",shelfmin[i])]
  #Define growth seasons as day with phytoplankton concentrations >500 μg phytoplankton C/m3, 
  startday<-which(phyto_daymean_y>500)[1]
  endday<-max(which(phyto_daymean_y>500))
  growthseason_length<-endday-startday
  #Duration per stage if developing from N3 to halfway through C3 in growth season:
  pred_c[,paste0("dev_growthseason_",y)]<-NA
  pred_c[,paste0("dev_growthseason_",y)][2:8]<-growthseason_length*pred_c$rel_dev_stages[2:8]
}

#4. Find carbon needed for growth:
#Increase in carbon weight from one stage to the next:
diff_c<-diff(pred_c$C_µg)
pred_c<-pred_c[pred_c$age>=3,] #Remove stage N2
pred_c<-pred_c[,-c(4)] #Remove unnecessary columns

#Daily carbon needed for growth for stage i=(weighti-weight(i-1))/dev.time i:
for(y in c(1987,2012)){
  pred_c[,paste0("growth_µg_C_day_",y)]<-diff_c/pred_c[,paste0("dev_growthseason_",y)]
}


#5. Take into consideration filtration rates:
#Range of observed filtration rates:
filtration<-seq(min(filtrations$F_m3_ind_day),max(filtrations$F_m3_ind_day),by=0.000001)


####How many days are the requirements exceeded?####
layout(matrix(1:2,ncol=1),heights=c(0.62,0.38))
par(mar=c(0.2,0,0,0),oma=c(1.5,1.5,0.1,0.1),ps=8,mgp=c(3,0.2,0))

plot(x=NULL,y=NULL,ylim=c(360,5),xlim=c(min(filtration),max(filtration-0.00002)),axes=FALSE,frame=TRUE)
Axis(side=2,at=c(1,100,200,300,365))

for(y in c(1987,2012)){
  #Energy available per day:
  phyto_daymean_y<-phyto_daymean[,paste0("phyto_",y,"_lat>",latmin[i],".shelf>",shelfmin[i])]
  #Cumulative energy available per day:
  phyto_cumsum<-cumsum(phyto_daymean_y)
  startday<-which(phyto_daymean_y>500)[1]
  endday<-max(which(phyto_daymean_y>500))
  growthseason_length<-endday-startday
  #Plot outline with growth season:
  if(y==1987){
    lty=1;f=0.4
    }else{
    lty<-2 ;f=0.2   }
  rect(-1,startday,1,endday,col=adjustcolor("#009E73",alpha.f=f),border=NA);box()
  #Estimate total energy required summing all stages:
  #Growth+metabolism for full development time of stage
  tot_C_req<-sum((pred_c$meta_µg_C_day+pred_c[,paste0("growth_µg_C_day_",y)])*
                   pred_c[,paste0("dev_growthseason_",y)])
  #Total energy requirement adjusted per filtration rate:
  tot_C_req_filtration<-tot_C_req/filtration
  #The day the requirements are met per filtration rate:
  day_reach<-sapply(tot_C_req_filtration,function(x) min(which(phyto_cumsum>=x)))
  day_reach[is.infinite(day_reach)]<-NA
  points(filtration,day_reach,type="l",lty=lty)
  #Mark below which filtration rate requirement cannot be met:
  points(max(filtration[day_reach==max(day_reach,na.rm=TRUE)],na.rm=TRUE),
         max(day_reach,na.rm=TRUE),pch=4,cex=0.8)
  #4: With microzooplankton included as food source:
  zoop_daymean_y<-phyto_daymean[,paste0("zoop_",y,"_lat>",latmin[i],".shelf>",shelfmin[i])]
  zoop_cumsum<-cumsum(zoop_daymean_y)
  day_reach<-sapply(tot_C_req_filtration,function(x) min(which(zoop_cumsum>=x)));day_reach[is.infinite(day_reach)]<-NA
  points(filtration,day_reach,type="l",col="dark grey",lty=lty)
  points(max(filtration[day_reach==max(day_reach,na.rm=TRUE)],na.rm=TRUE),
         max(day_reach,na.rm=TRUE),pch=4,col="dark grey",cex=0.8)
  }
#Add legends:
legend("topleft",legend=c("1987","2012","With microzooplankton food"),y.intersp=0.7,
                   lty=c(1,2,1),col=c("black","black","dark grey"),bty="n")
legend("topright",legend=c(1987,2012),title="Growth season",bty="n",y.intersp=0.7,
                pch=20,col=c(adjustcolor("#009E73",alpha.f=0.4),adjustcolor("#009E73",alpha.f=0.2)))
mtext(side=2,outer=FALSE,"Day-of-year energy intake > requirement",line=1)

#Plot filtration rates per stage:
plot(x=NULL,y=NULL,ylim=c(1,5.2),xlim=c(min(filtration),max(filtration-0.00002)),axes=FALSE,frame=TRUE)
stripchart(filtrations$F_m3_ind_day~filtrations$Stage,add=T,pch=20,cex=.5)
axis(2, at=1:5, labels=unique(filtrations$Stage))
axis(1,at=c(0,2e-4,4e-4,6e-4),labels=10^4*c(0,2e-4,4e-4,6e-4) )
abline(h=6e-6,lty=2);abline(h=0.00034,lty=2)
mtext(side=1,outer=TRUE,expression(paste("Filtration rate (10"^"-4 ","m"^"3 ","ind."^"-1","day"^"-1",")")),line=1)
mtext(side=2,outer=FALSE,"Stage",line=1)


