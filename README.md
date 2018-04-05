#  chyp_arctic 

<b>R code (rcode folder) and input files (data folder) for analyses of <i>Calanus hyperboreus</i> distribution and energy requirements.</b>

The rcode folder contains the following R scripts:<br>
gams_chyp.r: Fit GAMs of <i>Calanus hyperboreus</i> abundance data<br> 
c3_seasonal.r: Plot seasonal variation in <i>Calanus hyperboreus</i> C3 concentrations against growth season<br>
energy_budget.r: Estimate <i>Calanus hyperboreus</i> energy requirements to develop from N3-C3<br>
... plus additional files to plot polar-centered maps, draw color image scales and calculate leave-one-year-out cross validation<br>


The data folder contains the following R data files:<br>
chyp_orig.rda: Original <i>Calanus hyperboreus</i> abundance data in ind.m2 or ind.m3 <br>   
chyp_depthint.rda: Depth integrated <i>Calanus hyperboreus</i> abundance data in ind.m2<br>     
carb_weight_stages.rda: Stage-specific body weights (see energy_budget.r)<br>           
filtrations.rda: Stage-specific filtration rates (see energy_budget.r)<br>     
phyto_daymean.rda: Averaged daily phytoplankton/phyto+microzooplankton concentrations (see energy_budget.r/c3_seasonal.r)<br>          
predframe.rda: Dataframe for GAM predictions (see gams_chyp.r)<br> 