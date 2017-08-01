library(ggplot2)
library(plyr)
library(corrplot)
library(nlme)
library(plotrix)


source("/home/nadkarni/git/sammba-mri/sammba/common/variousRfunctions.R")
projectdir = "/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/AmylNet"


###MAKE REFTABLE###
setwd("/home/Pmamobipet/Tvx-Manips-MD_/MD_1602-AmylNet-Garin/MRIanalyses/texteAmylNet")
RefTablegenerator("MouseIDs.txt",
                  "StudyCages.txt",
                  "Treatments.txt",
                  "MRISessions6.txt",
                  "DOD.txt",
                  paste(projectdir, "RefTable.txt", sep = "/")
                  )


###LOAD REGION NAMES###
regions = read.table("/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/c57_brain_atlas_labels.txt",
                     header=T, sep="\t"
)


setwd(projectdir)


###READ REFTABLE###
RefTable = read.table("RefTable.txt", header=T, sep = "\t")
RefTable$Treatment = sub("APPPS1D9", "APP", RefTable$Treatment, fixed = TRUE)
RefTable$Treatment = factor(RefTable$Treatment, levels=c("WT", "APP", "0"))
RefTable$Session = as.factor(RefTable$Session)
RefTable$Genotype = sub("APPPS1D9", "APP", RefTable$Genotype, fixed = TRUE)
RefTable$Genotype = factor(RefTable$Genotype, levels=c("wt", "APP"))


setwd("analysis20170725")


## RBM. ASSUMES ALL ANATOMICAL SCANS WERE OF SAME FOV AND MATRIX!!!

RBMdata = RBMdataproc("MRIsessions", "atlas_Na1.nii.gz", regions, RefTable)
RBMdata$volume = RBMdata$count * 0.15625 * 0.104167 * 0.15625

pdf(file = "RBM.pdf", width = 15* sqrt(2), height = 15)
  qplot(Treatment, volume, data = subset(RBMdata, side == "both"),
        geom = "boxplot", outlier.shape = NA, ymin = 0, colour = Session
        ) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5,
                                               dodge.width = 0.75
                                               )
               ) +
    facet_wrap(~ region, scales = "free") +
    ylab(expression(volume ~ (mm ^ {3})))
  qplot(Treatment, volume,
        data = subset(RBMdata,
                      side == "both" & 
                      {Session == "0" | Session == "1" | Session == "2"}
                      ),
        geom = "boxplot", outlier.shape = NA, ymin = 0, colour = Session
        ) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5,
                                               dodge.width = 0.75
                                               )
               ) +
    facet_wrap( ~ region, scales = "free") +
    ylab(expression(volume ~ (mm ^ {3})))
dev.off()

pdf(file = "RBMMLRS.pdf", width = 15 * sqrt(2), height = 15)
  qplot(Treatment, volume, data = subset(RBMdata,
                                         ScanDate > 20170600 &
                                         ScanDate > 20170700 &
                                         side == "both"
                                         ),
        geom = "boxplot", outlier.shape = NA, ymin = 0, colour = Genotype
        ) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5,
                                               dodge.width = 0.75
                                               )
               ) +
    facet_wrap( ~ region, scales = "free") +
    ylab(expression(volume ~ (mm ^ {3})))
dev.off()

write.table(RBMdatavoltab <- ddply(.data = subset(RBMdata,
                                                  region == "Total" &
                                                  side == "both"
                                                  ),
                                   c("Study", "Treatment", "Session"),
                                   summarise,
                                   N = length(unique(DICOMdir)),
                                   mean = mean(volume),
                                   sd = sd(volume),
                                   se = std.error(volume),
                                   median = median(volume),
                                   halfIQR = 0.5 * IQR(volume)
                                   ),
            file = "RBMTotaltab.txt", quote = F, sep = "\t", row.names = F
            )
write.table(RBMdata, file = "RBMdata.txt", quote = F, sep = "\t", row.names = F)


##perfusion

TIs=c(35, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1600)
T1guess = 1600
fTotalexcludedregions = c("facial nerve (cranial nerve 7)",
                          "fourth ventricle",
                          "lateral olfactory tract",
                          "lateral ventricle",
                          "third ventricle"
                          )

perfdata = perfdataproc("MRIsessions", "perfFAIREPI_n[0-9]_NaMe.nii.gz",
                        "no", "yes", regions, RefTable, 2800, 0.9, TIs,
                        6000000, 1600, fTotalexcludedregions
                        )

pdf(file = "perfusion.pdf", width = 15 * sqrt(2), height = 15)
  qplot(Treatment, CBF, data = subset(perfdata, side == "both"),
        geom = "boxplot", outlier.shape = NA, ymin = 0, colour = Session) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5,
                                               dodge.width = 0.75)
               ) +
    facet_wrap( ~ region, scales = "free") +
    ylab(expression(CBF ~ (ml ~ 100 * g ^ {-1} ~ min ^ {-1}) ))
  qplot(Treatment, CBF, data = subset(perfdata, side == "both" & 
                                        {Session == "0" |
                                         Session == "1" |
                                         Session == "2"
                                         }
                                      ),
        geom = "boxplot", outlier.shape = NA, ymin = 0, colour = Session) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5,
                                               dodge.width = 0.75
                                               )
               ) +
    facet_wrap(~ region, scales = "free") +
    ylab(expression(CBF ~ (ml ~ 100 * g ^ {-1} ~ min ^ {-1}) ))
dev.off()

perfdatan = perfdataproc("MRIsessions", "perfFAIREPI_n[0-9]_NaMe.nii.gz",
                         "no", "no", regions, RefTable, 2800, 0.9, TIs,
                         6000000, 1600, fTotalexcludedregions
                         )
perfdatan$perforder = gsub("_NaMe.nii.gz", "",
                           gsub("perfFAIREPI_n", "", basename(perfdatan$File))
                           )

write.table(perfdatatab <- ddply(.data = subset(perfdata, 
                                                region == "fTotal" &
                                                  side == "both"
                                                ), 
                                 c("Study", "Treatment", "Session"),
                                 summarise,
                                 N = length(unique(DICOMdir)),
                                 mean = mean(wCBF),
                                 sd = sd(wCBF),
                                 se = std.error(wCBF),
                                 median = median(wCBF),
                                 halfIQR = 0.5 * IQR(wCBF)
                                 ),
            file = "perfTotaltab.txt", quote = F, sep = "\t", row.names = F
            )

write.table(perfdata, file = "perfdata.txt",
            quote = F, sep = "\t", row.names = F
            )
write.table(perfdatan, file = "perfdatan.txt",
            quote = F, sep = "\t", row.names = F
            )

png(file = "diagnostic.png", width = 1500 * sqrt(2), height = 1500, res = 300)
  qplot(Genotype, wCBF, data = subset(perfdata,
                                      ScanDate > 20170600 &
                                        region == "fTotal" &
                                        side == "both"
                                      ),
        geom = "boxplot", outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.25)) +
  geom_text(aes(label = MouseID), size = 1, position = position_nudge(x = 0.5))
dev.off()


#rsCorr to modify
#rsdata=rsdataproc("analysis20170725/MRIsessions","rs_n[0-9]_Av_NaMeBlErvent.nii.gz",regions,RefTable)

#rsdata$rsorder=gsub("rs_n","",gsub("_Av_NaMeBlErvent.nii.gz","",basename(as.character(rsdata$File))))
#rsdata=cbind(cbind(rsdata[,1:(ncol(RefTable)+3)],rsdata["rsorder"],rsdata[,(ncol(RefTable)+4):(ncol(rsdata)-1)]))

#rsdatanames=substr(names(rsdata)[(ncol(RefTable)+5):ncol(rsdata)],6,8)
#for (n in 1:nrow(regions)) {
#  rsdatanames[match(regions$label[n],rsdatanames)]=paste(regions$side[n],regions$structure[n],sep=" ")
#}
#rsdatanames=gsub("both ","",rsdatanames)
#names(rsdata)[(ncol(RefTable)+5):ncol(rsdata)]=rsdatanames
