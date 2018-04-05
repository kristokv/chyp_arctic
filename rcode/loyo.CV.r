#R code to calculate Leave-One-Year-Out (LOYO) cross validation
#Created by LC Stige, UiO


## Loyo (leave-one-year-out) cv:
loyo.CV<-function(form,				# form = model formula
			type, 			            # type = model type ("gam", "threshold.gam", etc.)
			data, 			            # data = name of data.frame
			grouping.variable="Year", # Grouping variable, e.g. year: predictions are made for one year at a time, based on a model fitted to data without the given year  
			family="gaussian",      # family argument
      random=NULL,cor=NULL, 	# Arguments used in gamm models
			a=.2,b=.8,d=.4,nthd=30,threshold.name=NULL,steps="even",	# Arguments used in threshold models
			thd.v1.name=NULL, thd.v2.name=NULL,d.00 = 0.2, d.10 = 0.2, d.01= 0.2, d.11 = 0.2,                  
			output="mean",		      # If "mean", sqrt of mean squared prediction error is printed. If "x", predicted and observed values are printed.
			outlier="none"		      # If "null", the categorical "Outlier" variable in the prediction data set is set to 0.
	) {
	X <- data.frame(NULL)
	response.name <- all.vars(form)[1]
	data$obs <- eval(parse(text = paste("data", response.name,sep = "$")))
	data$group <- eval(parse(text = paste("data", grouping.variable,sep = "$")))
	yrs=levels(as.factor(as.character(data$group)))
	for(yr in as.character(yrs)){
	  DataTrain <- data[data$group!=yr,] #Leave the year out
	  DataNew <- data[data$group==yr,] #Dataset for predition = data for year yr
	if(dim(DataNew)[1]>0){
        if(outlier=="null"){	DataTrain$Outlier <- as.factor(as.character(DataTrain$Outlier))
					DataNew$Outlier <- 0
					levels(DataNew$Outlier)<-levels(DataTrain$Outlier)}
	if(type=="lm"){                             # lm-model
	  mod.new <- lm(form, data=DataTrain)
	  pred<- predict(mod.new, newdata=DataNew)     }
	if(type=="gamm"){                           # GAMM with no thresholds
	  mod.new <- gamm(form,family=family,random=random,cor=cor, data=DataTrain)
	  pred<- predict.gam(mod.new$gam, newdata=DataNew, type="response")     }
	if(type=="gam"){                             # no estimation of thresholds
	  mod.new <- gam(form,family=family, data=DataTrain)
	  pred<- predict.gam(mod.new, newdata=DataNew, type="response")     }
	if(type=="threshold.gam"){                             # estimation of 1 threshold
	  mod.new <- threshold.gam(form,family=family, data=DataTrain,a=a,b=b,nthd=nthd,threshold.name=threshold.name,steps=steps)
	  DataNew$r <- mod.new$mr
	  pred<- predict.gam(mod.new$res, newdata=DataNew, type="response")     }
	if(type=="bivariate.thd.gam"){                             # estimation of bivariate threshold
	  mod.new <- bivariate.thd.gam(form,family=family, data=DataTrain, a=a, thd.v1.name=thd.v1.name, thd.v2.name=thd.v2.name)
	  DataNew$thd.v1 <- eval(parse(text = paste("DataNew", thd.v1.name,sep = "$")))
	  DataNew$thd.v2 <- eval(parse(text = paste("DataNew", thd.v2.name,sep = "$")))
	  DataNew$score <- (DataNew$thd.v1 - mod.new$point1[1]) * (mod.new$point2[2] - mod.new$point1[2]) -
                (DataNew$thd.v2 - mod.new$point1[2]) * (mod.new$point2[1] - mod.new$point1[1])
	  pred<- predict.gam(mod.new$res, newdata=DataNew, type="response")     }
	if(type=="threshold.gamm"){                  # estimation of gamm with 1 threshold
	  mod.new <- threshold.gamm(form,family=family, random=random, cor=cor, data=DataTrain,a=a,b=b,nthd=nthd,threshold.name=threshold.name,steps=steps)
	  DataNew$r <- mod.new$mr
	  pred<- predict.gam(mod.new$res$gam, newdata=DataNew, type="response")     }
	if(type=="two.thr.gam"){                     # estimation of 2 thresholds
	  mod.new <- two.thr.gam(form,family=family, data=DataTrain,a=a,b=b,d=d,nthd=nthd,threshold.name=threshold.name)
	  DataNew$r1 <- mod.new$mr1
	  DataNew$r2 <- mod.new$mr2
	  pred<- predict.gam(mod.new$res, newdata=DataNew, type="response")     }
	if(type=="two.tvar.gam"){                     # estimation of 2 thresholds in 2 diff. variables
	  mod.new <- two.tvar.gam(form,family=family, data=DataTrain,a=a,b=b,d.00 = d.00, d.10 = d.10, d.01= d.01, d.11 = d.11,nthd=nthd,thd.v1.name=thd.v1.name, thd.v2.name=thd.v2.name,steps=steps)
	  DataNew$r1 <- mod.new$mr1
	  DataNew$r2 <- mod.new$mr2
	  pred<- predict.gam(mod.new$res, newdata=DataNew, type="response")     }
	if(type=="bivariate.thd.gamm"){                             # estimation of bivariate threshold gamm
	  mod.new <- bivariate.thd.gamm(form,family=family, random=random, cor=cor, data=DataTrain, a=a, thd.v1.name=thd.v1.name, thd.v2.name=thd.v2.name)
	  DataNew$thd.v1 <- eval(parse(text = paste("DataNew", thd.v1,sep = "$")))
	  DataNew$thd.v2 <- eval(parse(text = paste("DataNew", thd.v2,sep = "$")))
	  DataNew$score <- (DataNew$thd.v1 - mod.new$res$point1[1]) * (mod.new$res$point2[2] - mod.new$res$point1[2]) -
                (DataNew$thd.v2 - mod.new$res$point1[2]) * (mod.new$res$point2[1] - mod.new$res$point1[1])
	  pred<- predict.gam(mod.new$res$gam, newdata=DataNew, type="response")     }
	if(type=="threshold.lm"){                    # estimation of 1 threshold in linear model
	  mod.new <- threshold.lm(form, data=DataTrain,a=a,b=b,nthd=nthd,threshold.name=threshold.name,steps=steps)
	  DataNew$r <- mod.new$mr
	  pred<- predict(mod.new$res, newdata=DataNew)     }
	if(type=="two.thr.lm"){                      # estimation of 2 thresholds in linear model
	  mod.new <- two.thr.lm(form, data=DataTrain,a=a,b=b,d=d,nthd=nthd,threshold.name=threshold.name)
	  DataNew$r1 <- mod.new$mr1
	  DataNew$r2 <- mod.new$mr2
	  pred<- predict(mod.new$res, newdata=DataNew)     }
	if(type=="three.thr.lm"){                      # estimation of 3 thresholds in linear model
	  mod.new <- three.thr.lm(form, data=DataTrain,a=a,b=b,d12=d12,d23=d23,nthd=nthd,threshold.name=threshold.name)
	  DataNew$r1 <- mod.new$mr1
	  DataNew$r2 <- mod.new$mr2
	  DataNew$r3 <- mod.new$mr3
	  pred<- predict(mod.new$res, newdata=DataNew)     }
	X <- rbind.data.frame(X, data.frame(Group=DataNew$group, Obs=DataNew$obs, Pred=pred))
}  }
if(output=="mean"){out <- signif(mean(apply(cbind(X$Obs,X$Pred),1,function(x){(x[2]-x[1])^2})),4)} # Mean squared prediction error
if(output=="mean"&is.element(family,c("poisson","quasipoisson"))){out <- signif(mean(apply(cbind(log(X$Obs+1),log(X$Pred+1)),1,function(x){(x[2]-x[1])^2}))^.5,4)}
if(output=="misclass.rate" & is.element(family,c("binomial","quasibinomial"))){
  out <- signif(mean(apply(cbind(X$Obs,X$Pred),1,function(x){abs((1*(x[2]>0.5))-(1*(x[1]>0)))})),4)}  # Binomial misclassification rate
if(output=="x"){ out <- X }
out
	}
