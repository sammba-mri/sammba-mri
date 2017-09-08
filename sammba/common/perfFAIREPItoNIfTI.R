#!/usr/bin/Rscript

#args[1] = perffile
#args[2] = T1blood
#args[3] = lambda
#args[4] = multiplier
#args[5] = T1guess
#args[6] = mc.cores
#args[7] = outputfile

args<-commandArgs(TRUE)

source('/home/nadkarni/git/sammba-mri/sammba/common/variousRfunctions.R')

writeNIfTI(perfFAIREPItoNIfTI(args[1], as.numeric(args[2]), as.numeric(args[3]),
                              as.numeric(args[4]), as.numeric(args[5]),
                              as.numeric(args[6])),
           args[7])
