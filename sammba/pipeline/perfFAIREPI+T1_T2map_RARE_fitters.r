library(oro.nifti)
library(minpack.lm)
library(parallel)


#lambda is the assumed blood tissue partition coefficient of water in ml/g
#if multiplier=1, CBF units are ml g-1 ms-1
#6000000 converts to ml per 100g per min

#fitter notes:
#mean(S0s) seems to work, 1 does not, bizarrely 0 does, max(S0s) does, maybe this would be better?
#fit equation does not work if bias is within abs. seems to work with nls using port or even just the default Gauss-Newton 
#(which does not accept bounds), maybe they would be better and/or faster?

#T1guess=prod(TIs)^(1/length(TIs)) does not work at high T1
#UPDATE: not sure if tested that with the newer robust LM equation below
#this code is not adaptive to T1, tend to just use max(TIs), need to make it more flexible


NIfTIreader=function(NIfTIfile){
  #ought to include error function in case TRs or TEs is the wrong length but whatever
  readNIfTI(file.path(normalizePath(NIfTIfile)), reorient=F)
}


#not extensively tested. deals with the weird order of the DICOM-converted data in T1_T2map_RARE
#in the T1_T2map_RARENIfTI, the third dimension is for TEs not Z. The fourth dimension loops Z within TRs
#the first array function converts the 4D NIfTI to 5D, aperm then swaps TEs and Z, before a second array function returns it to 4D 
T1_T2map_RAREreorder=function(T1_T2map_RARENIfTI,TRs,TEs){
  array(aperm(array(T1_T2map_RARENIfTI,dim=c(dim(T1_T2map_RARENIfTI)[1:3],length(TEs),length(TRs))),perm=c(1,2,4,3,5)),dim=dim(T1_T2map_RARENIfTI))
}


#generate list of co-ordinates for lapply/mclapply/mclapply2 to loop through
coordlistgen=function(NIfTI) {
  dims3=dim(NIfTI)[1:3]
  allcoords=array(0,c(dim=prod(dims3),4))
  allcoords[,1]=1:dim(allcoords)[1]
  allcoords[,2]=rep(1:dims3[1],dims3[2]*dims3[3])
  allcoords[,3]=rep(rep(1:dims3[2],each=dims3[1]),dims3[3])
  allcoords[,4]=rep(1:dims3[3],each=dims3[1]*dims3[2])
  allcoords=data.frame(allcoords)
  names(allcoords)=c("index","xd","yd","zd")
  split(allcoords,allcoords$index)
}


perfT1fitter=function(S0s,TIs,T1guess) {
	tryCatch(summary(nlsLM(S0s ~ bias + abs ( M0 * (1-2*exp(-TIs/T1)) ), start=list(bias=0, M0=mean(S0s), T1=T1guess), lower=c(0,0,0))), error=function(cond) return(NA))
}


perfFAIREPIfitter=function(T1blood,lambda,TIs,multiplier,T1guess,selS0s,nonselS0s) {

	sel=perfT1fitter(selS0s,TIs,T1guess)
	nonsel=perfT1fitter(nonselS0s,TIs,T1guess)	

	if (!is.list(sel) | !is.list(nonsel)) {
		NA
  }
	else {
		T1sel=sel[["coefficients"]]["T1","Estimate"]
		T1nonsel=nonsel[["coefficients"]]["T1","Estimate"]
		rCBF=100*(T1nonsel-T1sel)/T1nonsel
		CBF=multiplier * lambda * (T1nonsel/T1blood) * ( (1/T1sel) - (1/T1nonsel) )
		list(sel=sel,nonsel=nonsel,rCBF=rCBF,CBF=CBF)
	}
	
}


T1T2fitter=function(fitT1orT2,S0s,TRsorTEs,T1orT2guess) {
  if (fitT1orT2=="T1") T1orT2fit=tryCatch(summary(nlsLM(S0s ~ bias + M0 * (1-exp(-TRsorTEs/T1)), start=list(bias=0, M0=mean(S0s), T1=T1orT2guess), lower=c(0,0,0))), error=function(cond) return(NA))
  if (fitT1orT2=="T2") T1orT2fit=tryCatch(summary(nlsLM(S0s ~ bias + M0 * exp(-TRsorTEs/T2), start=list(bias=0, M0=mean(S0s), T2=T1orT2guess), lower=c(0,0,0))), error=function(cond) return(NA))
  if (!is.list(T1orT2fit)) NA else T1orT2fit
}


perfFAIREPIvoxel=function(coord,selselector,nonselselector,perfFAIREPINIfTI,T1blood,lambda,TIs,multiplier,T1guess) {
  perfFAIREPIfitter(T1blood,lambda,TIs,multiplier,T1guess,perfFAIREPINIfTI[coord[["xd"]],coord[["yd"]],coord[["zd"]],selselector],perfFAIREPINIfTI[coord[["xd"]],coord[["yd"]],coord[["zd"]],nonselselector])
}


T1_T2map_RAREvoxel=function(coord,selector,T1_T2map_RARENIfTI,fitT1orT2,TRsorTEs,T1orT2guess) {
  T1T2fitter(fitT1orT2,T1_T2map_RARENIfTI[coord[["xd"]],coord[["yd"]],coord[["zd"]],selector],TRsorTEs,T1orT2guess)
}


T1_T2map_RAREcalcparamextract=function(T1_T2map_RAREcalc,param,paramtype) {
  tryCatch(T1_T2map_RAREcalc[["coefficients"]][param,paramtype], error=function(cond) return(NA) )
}


T1_T2map_RAREtoNIfTI=function(NIfTI,TRs,TEs,T1guess,T2guess,mc.cores) {
  
  NIfTI=T1_T2map_RAREreorder(NIfTIreader(NIfTI),TRs,TEs)
  coordlist=coordlistgen(NIfTI)
  
  T1calc=mclapply(coordlist,T1_T2map_RAREvoxel,selector=seq(1,length(TRs)*length(TEs),length(TEs)),T1_T2map_RARENIfTI=NIfTI,fitT1orT2="T1",TRsorTEs=TRs,T1orT2guess=T1guess,mc.cores=mc.cores)
  T2calc=mclapply(coordlist,T1_T2map_RAREvoxel,selector=1:length(TEs),T1_T2map_RARENIfTI=NIfTI,fitT1orT2="T2",TRsorTEs=TEs,T1orT2guess=T2guess,mc.cores=mc.cores)

  #this is disgusting, find a better way
  T1bias      =sapply(T1calc, T1_T2map_RAREcalcparamextract, "bias", "Estimate")
  T1biasSE    =sapply(T1calc, T1_T2map_RAREcalcparamextract, "bias", "Std. Error")
  T1M0        =sapply(T1calc, T1_T2map_RAREcalcparamextract, "M0",   "Estimate")
  T1M0SE      =sapply(T1calc, T1_T2map_RAREcalcparamextract, "M0",   "Std. Error")
  T1T1        =sapply(T1calc, T1_T2map_RAREcalcparamextract, "T1",   "Estimate")
  T1T1SE      =sapply(T1calc, T1_T2map_RAREcalcparamextract, "T1",   "Std. Error")
  T2bias      =sapply(T2calc, T1_T2map_RAREcalcparamextract, "bias", "Estimate")
  T2biasSE    =sapply(T2calc, T1_T2map_RAREcalcparamextract, "bias", "Std. Error")
  T2M0        =sapply(T2calc, T1_T2map_RAREcalcparamextract, "M0",   "Estimate")
  T2M0SE      =sapply(T2calc, T1_T2map_RAREcalcparamextract, "M0",   "Std. Error")
  T2T2        =sapply(T2calc, T1_T2map_RAREcalcparamextract, "T2",   "Estimate")
  T2T2SE      =sapply(T2calc, T1_T2map_RAREcalcparamextract, "T2",   "Std. Error")
  
  as.nifti(array(c(T1bias,T1biasSE,T1M0,T1M0SE,T1T1,T1T1SE,T2bias,T2biasSE,T2M0,T2M0SE,T2T2,T2T2SE),dim=c(dim(NIfTI)[1:3],12)))
  
}


perfFAIREPIcalcparamextract=function(perfFAIREPIcalc,selornonsel,param,paramtype) {
  tryCatch(perfFAIREPIcalc[[selornonsel]][["coefficients"]][param,paramtype], error=function(cond) return(NA) )
}


perfFAIREPItoNIfTI=function(NIfTI,T1blood,lambda,TIs,multiplier,T1guess,mc.cores) {
  
  NIfTI=NIfTIreader(NIfTI)
  coordlist=coordlistgen(NIfTI)
  
  selselector=seq(1,dim(NIfTI)[4],2)
  nonselselector=seq(2,dim(NIfTI)[4],2)
    
  perfcalc=mclapply(coordlist,perfFAIREPIvoxel,selselector=selselector,nonselselector=nonselselector,perfFAIREPINIfTI=NIfTI,T1blood=T1blood,lambda=lambda,TIs=TIs,multiplier=multiplier,T1guess=T1guess,mc.cores=mc.cores)

  selbias      =sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "bias", "Estimate")
  selbiasSE    =sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "bias", "Std. Error")
  selM0        =sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "M0",   "Estimate")
  selM0SE      =sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "M0",   "Std. Error")
  selT1        =sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "T1",   "Estimate")
  selT1SE      =sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "T1",   "Std. Error")
  nonselbias   =sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "bias", "Estimate")
  nonselbiasSE =sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "bias", "Std. Error")
  nonselM0     =sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "M0",   "Estimate")
  nonselM0SE   =sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "M0",   "Std. Error")
  nonselT1     =sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "T1",   "Estimate")
  nonselT1SE   =sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "T1",   "Std. Error")
  rCBF         =sapply(perfcalc, function(x) tryCatch(x[["rCBF"]], error=function(cond) return(NA) ))
  CBF          =sapply(perfcalc, function(x) tryCatch(x[["CBF"]], error=function(cond) return(NA) ))
  
  as.nifti(array(c(selbias,selbiasSE,selM0,selM0SE,selT1,selT1SE,nonselbias,nonselbiasSE,nonselM0,nonselM0SE,nonselT1,nonselT1SE,rCBF,CBF),dim=c(dim(NIfTI)[1:3],14)))
  
}


perfFAIREPIROI=function(perfFAIREPIdata,selS0sselector,nonselS0sselector,T1blood,lambda,TIs,multiplier,T1guess) {
  selS0s=perfFAIREPIdata[selS0sselector,]$NZMean
  nonselS0s=perfFAIREPIdata[nonselS0sselector,]$NZMean
  perfFAIREPIfitter(T1blood,lambda,TIs,multiplier,T1guess,selS0s,nonselS0s)
}
