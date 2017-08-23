#!/usr/bin/Rscript

#inputfile=args[1]
#T1blood=args[2]
#lambda=args[3]
#TIs=args[4]
#multiplier=args[5]
#T1guess=args[6]
#mc.cores=args[7]
#outputfile=args[8]

args<-commandArgs(TRUE)

source('/home/nadkarni/git/sammba-mri/sammba/common/variousRfunctions.R')

writeNIfTI(perfFAIREPItoNIfTI(args[1],
                              as.numeric(args[2]),
                              as.numeric(args[3]),
                              as.numeric(unlist(strsplit(args[4], ","))),
                              as.numeric(args[5]),
                              as.numeric(args[6]),
                              as.numeric(args[7])
                              ),
           args[8]
           )
