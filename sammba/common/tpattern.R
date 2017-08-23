#!/usr/bin/Rscript

##$ACQ_obj_order=( 5 )
#0 2 4 1 3
##$ACQ_obj_order=( 30 )
#0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29
#assume equal distribution of slice timings in the tpatterns

args<-commandArgs(TRUE)

TR=as.numeric(args[1])
slices=as.numeric(args[2])

interslicetime=TR/slices
evenslices=seq(0,slices-1,2)
oddslices=seq(1,slices-1,2)

write(c(evenslices*interslicetime,oddslices*interslicetime),file=args[3],ncol=slices)