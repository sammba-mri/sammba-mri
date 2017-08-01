library(oro.nifti)
library(minpack.lm)
library(parallel)
library(plyr)

#lambda is the assumed blood tissue partition coefficient of water in ml per g
#multiplier converts to desired units (usually ml per 100g per min)

#fitter notes:
#mean(S0s) seems to work, 1 does not, bizarrely 0 does, max(S0s) does, maybe
#this would be better?
#fit equation does not work if bias is within abs. seems to work with nls using
#port or even just the default Gauss-Newton (which does not accept bounds),
#maybe they would be better and/or faster?

#T1guess = prod(TIs) ^ (1 / length(TIs)) does not work at high T1
#UPDATE: not sure if tested that with the newer robust LM equation below
#this code is not adaptive to T1, tend to just use max(TIs), need to make it
#more flexible

RefTablegenerator = function(MouseIDs, StudyCages, Treatments, MRISessions,
                             DOD, outputfile
                             ) {
  
  MouseIDs = read.table(MouseIDs, header = T, sep = "\t")
  StudyCages = read.table(StudyCages, header = T, sep = "\t")
  Treatments = read.table(Treatments, header = T, sep = "\t")
  MRISessions = read.table(MRISessions, header = T, sep = "\t")
  DOD = read.table(DOD, header = T, sep = "\t")
  
  RefTable = merge(MouseIDs, StudyCages)
  RefTable = merge(RefTable, Treatments)
  RefTable = merge(RefTable, MRISessions)
  RefTable = merge(RefTable, DOD)
  
  RefTable <- RefTable[!duplicated(RefTable), ]
  
  RefTable$ScanDate = substring(basename(as.character(RefTable$DICOMdir)), 1, 8)
  RefTable$ScanTime = substring(basename(as.character(RefTable$DICOMdir)), 10, 15)
  
  RefTable = RefTable[, c("MouseID", "altID", "DOB", "DOD", "Gender",
                          "Genotype", "StudyCage", "Study", "Treatment",
                          "Session", "DICOMdir", "Notes", "ScanDate",
                          "ScanTime", "Weight"
                          )
                      ]
  
  write.table(RefTable, file = paste(normalizePath(dirname(outputfile)),
                                     basename(outputfile), sep="/"
                                     ),
              quote = FALSE, sep = "\t", row.names = FALSE
              )
  
}


NIfTIreader = function(NIfTIfile) {
  #no error function in case TRs or TEs is the wrong length, whatevs
  readNIfTI(file.path(normalizePath(NIfTIfile)), reorient=F)
}


#not extensively tested. deals with the weird order of the DICOM-converted data
#in T1_T2map_RARE. in the T1_T2map_RARENIfTI, the third dimension is for TEs
#not Z. The fourth dimension loops Z within TRs. the first array function
#converts the 4D NIfTI to 5D, aperm then swaps TEs and Z, before a second array
#function returns it to 4D 
T1_T2map_RAREreorder = function(T1_T2map_RARENIfTI, TRs, TEs) {
  array(aperm(array(T1_T2map_RARENIfTI,
                    dim = c(dim(T1_T2map_RARENIfTI)[1:3],
                            length(TEs),
                            length(TRs)
                            )
                    ),
                    perm = c(1,2,4,3,5)
              ),
          dim = dim(T1_T2map_RARENIfTI)
        )
}


#generate list of co-ordinates for lapply/mclapply/mclapply2 to loop through
coordlistgen=function(NIfTI) {
  dims3 = dim(NIfTI)[1:3]
  allcoords = array(0, c(dim = prod(dims3), 4))
  allcoords[, 1] = 1:dim(allcoords)[1]
  allcoords[, 2] = rep(1:dims3[1], dims3[2] * dims3[3])
  allcoords[, 3] = rep(rep(1:dims3[2], each = dims3[1]), dims3[3])
  allcoords[, 4] = rep(1:dims3[3], each = dims3[1] * dims3[2])
  allcoords=data.frame(allcoords)
  names(allcoords) = c("index", "xd", "yd", "zd")
  split(allcoords, allcoords$index)
}


perfT1fitter = function(S0s, TIs, T1guess) {
  tryCatch(summary(nlsLM(S0s ~ bias + abs(M0 * (1 - 2 * exp(-TIs / T1))),
                         start = list(bias = 0, M0 = mean(S0s), T1 = T1guess),
                         lower = c(0, 0, 0)
                         )
                   ),
           error = function(cond) return(NA)
           )
}


perfFAIREPIfitter = function(T1blood,
                             lambda,
                             TIs,
                             multiplier,
                             T1guess,
                             selS0s,
                             nonselS0s
							 ) {

  sel=perfT1fitter(selS0s,TIs,T1guess)
  nonsel=perfT1fitter(nonselS0s,TIs,T1guess)    

  if (!is.list(sel) | !is.list(nonsel)) {
    NA
  }
  else {
    T1sel = sel[["coefficients"]]["T1","Estimate"]
    T1nonsel = nonsel[["coefficients"]]["T1","Estimate"]
    rCBF = 100 * (T1nonsel - T1sel) / T1nonsel
    CBF = multiplier *
          lambda *
          (T1nonsel / T1blood) *
          ((1 / T1sel) - (1 / T1nonsel))
    list(sel = sel, nonsel = nonsel, rCBF = rCBF, CBF = CBF)
    }
    
}


T1T2fitter = function(fitT1orT2, S0s, TRsorTEs, T1orT2guess) {
  if (fitT1orT2 == "T1") {
    T1orT2fit = tryCatch(summary(nlsLM(S0s ~ bias +
                                             M0 * (1 - exp(-TRsorTEs / T1)),
                                       start = list(bias = 0,
                                                    M0 = mean(S0s),
                                                    T1 = T1orT2guess
                                                    ),
                                       lower = c(0, 0, 0)
                                       )
                                 ),
                         error = function(cond) return(NA)
                         )
  }
  if (fitT1orT2 == "T2") {
    T1orT2fit = tryCatch(summary(nlsLM(S0s ~ bias + M0 * exp(-TRsorTEs / T2),
                                       start = list(bias = 0,
                                                    M0 = mean(S0s),
                                                    T2 = T1orT2guess
                                                    ),
                                       lower = c(0, 0, 0)
                                       )
                                 ),
                         error=function(cond) return(NA)
                         )
  if (!is.list(T1orT2fit)) NA else T1orT2fit
  }
}


perfFAIREPIvoxel = function(coord,
                            selselector,
                            nonselselector,
                            perfFAIREPINIfTI,
                            T1blood,
                            lambda,
                            TIs,
                            multiplier,
                            T1guess
                            ) {
  perfFAIREPIfitter(T1blood,
                    lambda,
                    TIs,
                    multiplier,
                    T1guess,
                    perfFAIREPINIfTI[coord[["xd"]],
                                     coord[["yd"]],
                                     coord[["zd"]],
                                     selselector
                                     ],
                    perfFAIREPINIfTI[coord[["xd"]],
                                     coord[["yd"]],
                                     coord[["zd"]],
                                     nonselselector
                                     ]
                    )
}

T1_T2map_RAREvoxel = function(coord,
                              selector,
                              T1_T2map_RARENIfTI,
                              fitT1orT2,
                              TRsorTEs,
                              T1orT2guess
                              ) {
  T1T2fitter(fitT1orT2,
             T1_T2map_RARENIfTI[coord[["xd"]],
                                coord[["yd"]],
                                coord[["zd"]],
                                selector
                                ],
             TRsorTEs,
             T1orT2guess
             )
}


T1_T2map_RAREcalcparamextract = function(T1_T2map_RAREcalc, param, paramtype) {
  tryCatch(T1_T2map_RAREcalc[["coefficients"]][param,paramtype],
           error=function(cond) return(NA)
           )
}


T1_T2map_RAREtoNIfTI = function(NIfTI, TRs, TEs, T1guess, T2guess, mc.cores) {
  
  NIfTI = T1_T2map_RAREreorder(NIfTIreader(NIfTI), TRs, TEs)
  coordlist = coordlistgen(NIfTI)
  
  T1calc = mclapply(coordlist,
                    T1_T2map_RAREvoxel,
                    selector = seq(1, length(TRs) * length(TEs), length(TEs)),
                    T1_T2map_RARENIfTI = NIfTI,
                    fitT1orT2 = "T1",
                    TRsorTEs = TRs,
                    T1orT2guess = T1guess,
                    mc.cores = mc.cores
                    )
                    
  T2calc = mclapply(coordlist,
                    T1_T2map_RAREvoxel,
                    selector = 1:length(TEs),
                    T1_T2map_RARENIfTI = NIfTI,
                    fitT1orT2 = "T2",
                    TRsorTEs = TEs,
                    T1orT2guess = T2guess,
                    mc.cores = mc.cores
                    )

  #this is disgusting, find a better way
  T1bias   = sapply(T1calc, T1_T2map_RAREcalcparamextract, "bias", "Estimate")
  T1biasSE = sapply(T1calc, T1_T2map_RAREcalcparamextract, "bias", "Std. Error")
  T1M0     = sapply(T1calc, T1_T2map_RAREcalcparamextract, "M0",   "Estimate")
  T1M0SE   = sapply(T1calc, T1_T2map_RAREcalcparamextract, "M0",   "Std. Error")
  T1T1     = sapply(T1calc, T1_T2map_RAREcalcparamextract, "T1",   "Estimate")
  T1T1SE   = sapply(T1calc, T1_T2map_RAREcalcparamextract, "T1",   "Std. Error")
  T2bias   = sapply(T2calc, T1_T2map_RAREcalcparamextract, "bias", "Estimate")
  T2biasSE = sapply(T2calc, T1_T2map_RAREcalcparamextract, "bias", "Std. Error")
  T2M0     = sapply(T2calc, T1_T2map_RAREcalcparamextract, "M0",   "Estimate")
  T2M0SE   = sapply(T2calc, T1_T2map_RAREcalcparamextract, "M0",   "Std. Error")
  T2T2     = sapply(T2calc, T1_T2map_RAREcalcparamextract, "T2",   "Estimate")
  T2T2SE   = sapply(T2calc, T1_T2map_RAREcalcparamextract, "T2",   "Std. Error")
  
  as.nifti(array(c(T1bias,
                   T1biasSE,
                   T1M0,
                   T1M0SE,
                   T1T1,
                   T1T1SE,
                   T2bias,
                   T2biasSE,
                   T2M0,
                   T2M0SE,
                   T2T2,
                   T2T2SE
                   ),
                 dim = c(dim(NIfTI)[1:3], 12)
                 )
           )
  
}


perfFAIREPIcalcparamextract = function(perfFAIREPIcalc,
                                       selornonsel,
                                       param,
                                       paramtype
									   ) {
  tryCatch(perfFAIREPIcalc[[selornonsel]][["coefficients"]][param,paramtype],
           error=function(cond) return(NA)
           )
}


perfFAIREPItoNIfTI = function(NIfTI,
                              T1blood,
                              lambda,
                              TIs,
                              multiplier,
                              T1guess,
                              mc.cores
							  ) {
  
  NIfTI = NIfTIreader(NIfTI)
  coordlist = coordlistgen(NIfTI)
  
  selselector = seq(1, dim(NIfTI)[4], 2)
  nonselselector = seq(2, dim(NIfTI)[4], 2)
    
  perfcalc = mclapply(coordlist,
                      perfFAIREPIvoxel,
                      selselector = selselector,
                      nonselselector = nonselselector,
                      perfFAIREPINIfTI = NIfTI,
                      T1blood = T1blood,
                      lambda = lambda,
                      TIs = TIs,
                      multiplier = multiplier,
                      T1guess = T1guess,
                      mc.cores = mc.cores
					  )

  selbias      = sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "bias", "Estimate")
  selbiasSE    = sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "bias", "Std. Error")
  selM0        = sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "M0",   "Estimate")
  selM0SE      = sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "M0",   "Std. Error")
  selT1        = sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "T1",   "Estimate")
  selT1SE      = sapply(perfcalc, perfFAIREPIcalcparamextract, "sel",    "T1",   "Std. Error")
  nonselbias   = sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "bias", "Estimate")
  nonselbiasSE = sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "bias", "Std. Error")
  nonselM0     = sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "M0",   "Estimate")
  nonselM0SE   = sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "M0",   "Std. Error")
  nonselT1     = sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "T1",   "Estimate")
  nonselT1SE   = sapply(perfcalc, perfFAIREPIcalcparamextract, "nonsel", "T1",   "Std. Error")
  rCBF         = sapply(perfcalc, function(x) tryCatch(x[["rCBF"]], error = function(cond) return(NA) ))
  CBF          = sapply(perfcalc, function(x) tryCatch(x[["CBF"]], error = function(cond) return(NA) ))
  
  as.nifti(array(c(selbias,
                   selbiasSE,
                   selM0,
                   selM0SE,
                   selT1,
                   selT1SE,
                   nonselbias,
                   nonselbiasSE,
                   nonselM0,
                   nonselM0SE,
                   nonselT1,
                   nonselT1SE,
                   rCBF,CBF
                   ),
                 dim=c(dim(NIfTI)[1:3], 14)
                 )
           )
}


perfFAIREPIROI = function(perfFAIREPIdata,
                          selS0sselector,
                          nonselS0sselector,
                          T1blood,
                          lambda,
                          TIs,
                          multiplier,
                          T1guess
                          ) {
  selS0s = perfFAIREPIdata[selS0sselector, ]$NZMean
  nonselS0s = perfFAIREPIdata[nonselS0sselector, ]$NZMean
  perfFAIREPIfitter(T1blood, lambda, TIs, multiplier, T1guess, selS0s, nonselS0s)
}


#####FORMERLY IN A SEPARATE DOCUMENT#####


readRBMdata = function(RBMfname) {
  iRBMdata = read.table(text = system(paste("/usr/lib/afni/bin/3dhistog -int", RBMfname),
                                      intern = TRUE, ignore.stderr = TRUE), 
                        col.names = c("label", "count", "CumFreq"))
  iRBMdata$CumFreq = NULL
  iRBMdata$fname = RBMfname
  iRBMdata
}


#using nz does not cause any bias, just excludes noisy data
readperfdata = function(perffname, mask) {
  if (mask == "no") mask = paste(dirname(perffname), 
                                 gsub("perfFAIREPI", 
                                      gsub("_NaMe.nii.gz", basename(perffname), 
                                           replacement = "_M0_N3.nii.gz"
                                           ), 
                                      replacement = "atlas_Na1_Op_perfFAIREPI"
                                      ), 
                                 sep = "/"
                                 )
  read.delim(text = system(paste("/usr/lib/afni/bin/3dROIstats -mask", mask,
                                 "-nzmean -nzvoxels", perffname),
                           intern = TRUE, ignore.stderr = TRUE
                           )
             )
}


readrsdata = function(rsfname) {
  mask = paste(dirname(rsfname), gsub("rs",
                                      gsub("_TsAv_NaMe.nii.gz",
                                           basename(rsfname),
                                           replacement = "_TsAvAvN3.nii.gz"
                                           ),
                                      replacement = "atlas_Na1_Op_rs"
                                      ),
               sep = "/"
               )
  read.delim(text=system(paste("/usr/lib/afni/bin/3dROIstats -mask",
                               mask, rsfname
                               ),
                         intern = TRUE, ignore.stderr = TRUE)
             )
}


cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}


RBMdataproc = function(MRIsessionspath, filepattern, regions, RefTable) {
  
  RBMfnames = list.files(path = MRIsessionspath, pattern = filepattern,
                         recursive = T, full.names = T
                         )
  RBMdata = llply(RBMfnames, readRBMdata, .progress = "text")
  
  RBMdata = do.call(rbind, RBMdata)
  
  RBMdata$region = RBMdata$label
  RBMdata$side = RBMdata$label
  
  for (n in 1:nrow(regions)) {
    label = regions$label[n]
    RBMdata$region[which(RBMdata$label == label)] = as.character(regions$structure[n])
    RBMdata$side[which(RBMdata$label == label)] = as.character(regions$side[n])
  }
  
  RBMdata = RBMdata[RBMdata$side == "left" |
                      RBMdata$side == "right" |
                      RBMdata$side == "both",
                    ] #gets rid of non-existent VOInumbers
  
  RefTable$dname = basename(as.character(RefTable$DICOMdir))
  RBMdata$dname = gsub("/", "",
                       gsub("/analyses", "",
                            gsub(filepattern, "",
                                 gsub(MRIsessionspath, "", RBMdata$fname)
                                 )
                            )
                       )
  
  RBMdata = merge(RefTable, RBMdata)
  
  RBMdata$fname = NULL
  
  RBMdataleftrightboth = ddply(.data = RBMdata[RBMdata$side != "both",],
                               c(names(RefTable), "region"),
                               summarise, count = sum(count)
                               )
  RBMdataleftrightboth$label = NA
  RBMdataleftrightboth$side = "both"
  RBMdata = rbind(RBMdata, RBMdataleftrightboth)
  
  RBMdatatotalboth = ddply(.data = RBMdata[RBMdata$side == "both",],
                           c(names(RefTable)), summarise, count=sum(count)
                           )
  RBMdatatotalboth$label = NA
  RBMdatatotalboth$region = "Total"
  RBMdatatotalboth$side = "both"
  RBMdatatotalleft = ddply(.data = RBMdata[RBMdata$side == "left",],
                           c(names(RefTable), "side"), summarise,
                           count=sum(count)
                           )
  RBMdatatotalleft$label = NA
  RBMdatatotalleft$region = "Total"
  RBMdatatotalright=ddply(.data = RBMdata[RBMdata$side == "right",],
                          c(names(RefTable), "side"), summarise,
                          count=sum(count)
                          )
  RBMdatatotalright$label = NA
  RBMdatatotalright$region = "Total"
  
  RBMdata = rbind(RBMdata, RBMdatatotalboth, RBMdatatotalleft, RBMdatatotalright)
  
  RBMdata
  
}


perfdataproc = function(MRIsessionspath, filepattern, mask, averageacross,
                        regions, RefTable, T1blood, lambda, TIs, multiplier,
                        T1guess, fTotalexcludedregions
                        ) {
  
  perffnames = list.files(path = MRIsessionspath, pattern = filepattern,
                          recursive = T, full.names = T
                          )
  perfdata = llply(perffnames, readperfdata, mask, .progress="text")
  
  perfdata = do.call(rbind.fill, perfdata)
  
  perfdata = reshape(perfdata, varying = 3:ncol(perfdata), timevar = "label",
                     direction = "long", sep="_"
                     )
  perfdata$id = NULL
  perfdata = split(perfdata, list(perfdata$File, perfdata$label))
  
  selS0sselector = seq(1, 2*length(TIs), 2)
  nonselS0sselector = seq(2, 2*length(TIs), 2)
  
  perfcalc=llply(perfdata, perfFAIREPIROI, selS0sselector, nonselS0sselector,
                 T1blood, lambda, TIs, multiplier, T1guess, .progress = "text") #parallel does not work yet
  
  #this is disgusting, need a better way to extract this data
  CBF = as.numeric(sapply(perfcalc, function(x) 
    tryCatch(x[["CBF", exact = FALSE]], 
             error = function(cond) return(NA)
             )
                          )
                   )
  File = as.character(sapply(perfdata, function(x) 
    tryCatch(x[["File", exact = FALSE]][1],
             error = function(cond) return(NA) 
             )
                             )
                      )
  label = as.character(sapply(perfdata, function(x)
    tryCatch(x[["label", exact = FALSE]][1],
             error = function(cond) return(NA)
             )
                              )
                       )
  vcount = as.numeric(sapply(perfdata, function(x)
    tryCatch(x[["NZcount", exact = FALSE]][1],
             error = function(cond) return(NA)
             )
                             )
                      )
  
  pd = data.frame(File, label, vcount, CBF)
  
  pd = pd[!is.na(pd$vcount) & !is.na(pd$CBF),]
  pd$File = as.character(pd$File)
  pd = subset(pd, ave(seq(label), label, FUN = length) >= 
                length(levels(as.factor(pd$File)))
              ) #remove labels not present in every acquisition
  
  pd$region = as.character(pd$label)
  pd$side = as.character(pd$label)
  
  for (n in 1:nrow(regions)) {
    label = regions$label[n] #label reused here, may change if above disgusting code is improved
    pd$region[which(pd$label == label)] = as.character(regions$structure[n])
    pd$side[which(pd$label == label)] = as.character(regions$side[n])
  }
  
  RefTable$dname = basename(as.character(RefTable$DICOMdir))
  pd$dname = gsub("/", "",
                  gsub("/analyses", "",
                       gsub(filepattern, "",
                            gsub(MRIsessionspath, "", pd$File)
                            )
                       )
                  )
  
  pd = merge(RefTable, pd)
  
  #voxel-weighted CBF; to combine left and right need to weight by VOI size as
  #left and right are not always symmetrical
  #also needed for any other VOI combinations that may be inserted after here
  pd$VCwCBF = pd$vcount * pd$CBF
  
  pdBLR = ddply(.data = pd[pd$side != "both",],
                c(names(RefTable), "dname", "File", "region"), summarise, 
                vcount = sum(vcount), CBF = mean(CBF), VCwCBF = sum(VCwCBF)
                )
  pdBLR$label = NA
  pdBLR$side = "both"
  pd = rbind(pd, pdBLR)
  
  pdBT = ddply(.data = pd[pd$side == "both",],
               c(names(RefTable), "dname", "File"), summarise,
               vcount = sum(vcount), CBF = mean(CBF), VCwCBF = sum(VCwCBF)
               )
  pdBT$label = NA
  pdBT$region = "Total"
  pdBT$side = "both"
  pdLT = ddply(.data = pd[pd$side == "left",],
               c(names(RefTable), "dname", "File", "side"), summarise,
               vcount = sum(vcount), CBF = mean(CBF), VCwCBF = sum(VCwCBF)
               )
  pdLT$label = NA
  pdLT$region = "Total"
  pdRT = ddply(.data = pd[pd$side == "right",],
               c(names(RefTable), "dname", "File", "side"), summarise,
               vcount = sum(vcount), CBF = mean(CBF), VCwCBF = sum(VCwCBF)
               )
  pdRT$label = NA
  pdRT$region = "Total"
  
  pdBfT = ddply(.data = pd[pd$side == "both" & !(pd$region %in% fTotalexcludedregions),],
                c(names(RefTable), "dname", "File"), summarise,
                vcount = sum(vcount), CBF = mean(CBF), VCwCBF = sum(VCwCBF)
                )
  pdBfT$label = NA
  pdBfT$region = "fTotal"
  pdBfT$side = "both"
  pdLfT=ddply(.data = pd[pd$side == "left" & !(pd$region %in% fTotalexcludedregions),],
              c(names(RefTable), "dname", "File", "side"), summarise,
              vcount = sum(vcount), CBF = mean(CBF), VCwCBF = sum(VCwCBF)
              )
  pdLfT$label = NA
  pdLfT$region = "fTotal"
  pdRfT = ddply(.data = pd[pd$side == "right" & !(pd$region %in% fTotalexcludedregions),],
                c(names(RefTable), "dname", "File", "side"), summarise,
                vcount = sum(vcount), CBF = mean(CBF), VCwCBF = sum(VCwCBF)
                )
  pdRfT$label = NA
  pdRfT$region = "fTotal"
  
  pd = rbind(pd, pdBT, pdLT, pdRT, pdBfT, pdLfT, pdRfT)
  pd$Session = as.factor(pd$Session)
  pd$wCBF = pd$VCwCBF / pd$vcount
  
  #average across all perfusion measures per session
  if (averageacross=="yes") {
    pd = ddply(pd, c(names(RefTable), "label", "region", "side"), summarise,
               vcount = mean(vcount),
               CBF = mean(CBF),
               VCwCBF = mean(VCwCBF),
               wCBF = mean(wCBF)
               )
  }
  #NEED: to add all the lost columns into above (or somehow make them automatically produce NAs),
  #then rbind to perfcalc. Call Session=0
  
  pd
  
}


rsdataproc = function(MRIsessionspath, filepattern, regions, RefTable) {
  rsnames = list.files(path = MRIsessionspath, pattern = filepattern,
                       recursive = T, full.names = T)
  rsdata = llply(rsfnames, readrsdata, .progress="text")
  rsdata = do.call(rbind.fill, rsdata)
  rsdata$dname = gsub("/", "",
                      gsub("/analyses", "",
                           gsub(filepattern, "",
                                gsub(MRIsessionspath, "", rsdata$File)
                                )
                           )
                      )
  RefTable$dname = basename(as.character(RefTable$DICOMdir))
  merge(RefTable, rsdata)
}
