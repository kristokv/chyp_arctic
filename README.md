#  chyp_arctic 

<b>R code (rcode folder) and input files (data folder) for analyses of <i>Calanus hyperboreus</i> distribution and energy requirements.</b>

The <b>rcode</b> folder contains the following R scripts:<br>
gams_chyp.r: Fit GAMs of <i>Calanus hyperboreus</i> abundance data<br> 
c3_seasonal.r: Plot seasonal variation in <i>Calanus hyperboreus</i> C3 concentrations against growth season<br>
energy_budget.r: Estimate <i>Calanus hyperboreus</i> energy requirements to develop from N3-C3<br>
... plus additional files to plot polar-centered maps, draw color image scales and calculate leave-one-year-out cross validation<br>


The <b>data</b> folder contains the following R data files:<br>
chyp_orig.rda: Original <i>Calanus hyperboreus</i> abundance data in ind.m2 or ind.m3<br>
chyp_depthint.rda: Depth integrated <i>Calanus hyperboreus</i> abundance data in ind.m2<br>
carb_weight_stages.rda: Stage-specific body weights (see <i>energy_budget.r</i>)<br>    
filtrations.rda: Stage-specific filtration rates (see <i>energy_budget.r</i>)<br>
phyto_daymean.rda: Averaged daily phytoplankton/phyto+microzooplankton concentrations (see <i>energy_budget.r</i>/<i>c3_seasonal.r</i>)<br>
predframe.rda: Dataframe for GAM predictions (see <i>gams_chyp.r</i>)<br> 