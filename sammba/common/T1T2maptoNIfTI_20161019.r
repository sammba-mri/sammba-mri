args<-commandArgs(TRUE)

source("/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/Code/perfFAIREPI+T1_T2map_RARE_fitters_20161019.r")

#NIfTI=args[1]
#TRs=args[2]
#TEs=args[3]
#T1guess=args[4]
#T2guess=args[5]
#mc.cores=args[6]
#outputfile=args[7]

writeNIfTI(T1_T2map_RAREtoNIfTI(args[1],as.numeric(unlist(strsplit(args[2],","))),as.numeric(unlist(strsplit(args[3],","))),as.numeric(args[4]),as.numeric(args[5]),as.numeric(args[6])),args[7])
