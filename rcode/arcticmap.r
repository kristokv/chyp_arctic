#R code to create polar centered map
#Created 2018 by Kristina Kvile (kristokv@gmail.com)

## Set the projection string and parameters. 
proj_str = 'orthographic'
proj_params = c() #no parameters for this projection
rotation = 0 # deg CW longitude to rotate the projection. 0 means E is down.

#Generate the basemap.
arcticmap<-function(addmap=TRUE,Xlim=c(-180, 180),Ylim=c(60, 90),col="black",map_labels=TRUE,map_grid=TRUE){
  if(addmap==TRUE){
  map('world',add=TRUE,projection=proj_str, parameters=proj_params,fill=T,border=col, col=col, orient=c(90, 0, rotation))
  }else{
  map('world',add=FALSE,projection=proj_str, parameters=proj_params,fill=T,border=col, col=col,xlim=Xlim, ylim=Ylim,orient=c(90, 0, rotation))
  }
  
  if(map_grid==TRUE){
#Add a grid. Plot first to get the full 360 deg lon
#Latitudes:
map.grid(c(0, 360, 50, 90), nx=0, ny=5, lty=1, col="grey", labels=FALSE,lwd=0.5,pretty=FALSE);
#Longitudes
map.grid(c(0, 360, 40, 90), nx=7, ny=0, lty=1, col="grey", labels=FALSE,lwd=0.5,pretty=FALSE)
}
if(map_labels==TRUE){
#Add grid labels:
#y-label:
xy = data.frame(x=seq(0,0,length.out=6),y=seq(40, 90,length.out=6))
xy$lab<-paste0(xy$y,"°N")
for (i in 3:length(xy$y)){
  lab<-xy$lab[i]
  with(mapproject(xy$x, xy$y, proj_str, parameters=proj_params, orientation=c(90, 0, rotation)),
       text(x[i],y[i], cex=0.75, label=lab))
}
#x-label (note: drawing on 40°N makes some labels disappear):
xy = data.frame(x=seq(0,360,length.out=7),y=seq(45,45,length.out=7))
xy$lab<-c(paste0(xy$x[1:4],"°E"),paste0(-(xy$x[5:7]-360),"°W"))
for (i in 1:(length(xy$y)-1)){
  lab<-xy$lab[i]
  with(mapproject(xy$x, xy$y, proj_str, parameters=proj_params, orientation=c(90, 0, rotation)),
       text(x[i],y[i], cex=0.75, label=lab))
}
}
}
