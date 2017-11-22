library(ggplot2)
library(readxl)

anatROIs = function(anatfile){
  afnitbl = read.table(text = system2("/usr/lib/afni/bin/3dhistog",
                                      c("-int", anatfile),
                                      stdout = TRUE, stderr = FALSE), 
                       col.names = c("label", "count", "CumFreq"))
  afnitbl$CumFreq = NULL
  afnitbl$fname = anatfile
  afnitbl$volume = afnitbl$count * 0.2 ^ 3
  afnitbl
}

anatomydir = "/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/MLOct2017/anatomy"

brains = lapply(list.files(anatomydir, pattern = '*brainmask.nii.gz', full.names = T), anatROIs)
atlases = lapply(list.files(anatomydir, pattern = '*atlas.nii.gz', full.names = T), anatROIs)
brains = do.call(rbind, brains)
atlases = do.call(rbind, atlases)
brains$type = "whole"
atlases$type = "regional"
h = rbind(brains, atlases)

labelfile = "/home/Pmamobipet/Tvx-Manips-MD_/MD_1704_RsLemurs/atlas/Region_pour_python.xlsx"
labeltable = read_xlsx(labelfile)
names(labeltable) = c("label", "blah2", "blah3", "blah4", "side", "blah6", "blah7", "region")
x = merge(h, labeltable)

projectdir = "/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/Tri_data"

IDs = read_xlsx(paste(projectdir, "Trilemur_ok.xlsx", sep = "/"), sheet = 1)
MRIsessions = read_xlsx(paste(projectdir, "Trilemur_ok.xlsx", sep = "/"), sheet = 2)
all = merge(IDs, MRIsessions)
all$age = all$date - all$DOB
all$group = "filler"
all$group[all$age < 2250] = "young"
all$group[all$age >= 2250] = "old"
all$group = factor(all$group, levels = c("young", "old"))
all$MRIsession = basename(all$PVdir)
x$MRIsession = sapply(strsplit(basename(as.character(x$fname)), "__"), "[[", 1)
newall = merge(all, x)

whole = subset(newall, type == "whole" & label == 1)

qplot(age, volume, data = whole)
qplot(group, volume, data = whole, geom = "boxplot", outlier.shape = NA) + geom_jitter(width = 0.2)

regional = subset(newall, type == "regional")

qplot(group, volume, data = subset(regional, label == 1),
      geom = "boxplot", outlier.shape = NA) +
  geom_jitter(width = 0.2)


a = qplot(group, volume, data = regional, geom = "boxplot", outlier.shape=NA) +
  geom_jitter(width = 0.2) +
  facet_wrap(~region, scales = "free_y")

pdf(file = "anat.pdf", width = 15* sqrt(2), height = 15)
a
dev.off()

for (regionname in labeltable$region) {
  capture.output(regionname, file = "grouptests.txt", append = TRUE)
  capture.output(t.test(volume ~ group, data = subset(regional, region == regionname)), file = "grouptests.txt", append = TRUE)
  capture.output(wilcox.test(volume ~ group, data = subset(regional, region == regionname)), file = "grouptests.txt", append = TRUE)  
}


png("whole.png", width = 10, height = 10, units = "cm", res = 300)
  qplot(group, volume, data = whole, geom = "boxplot", outlier.shape = NA) +
    geom_jitter(aes(colour = Sex), width = 0.2) +
    ylab(expression(volume ~ (mm ^ {3}))) +
    ggtitle("whole brain volume")
dev.off()

png("basalforebrain.png", width = 10, height = 5, units = "cm", res = 300)
qplot(group, volume, data = subset(regional, label %in% c(24, 25)),
      geom = "boxplot", outlier.shape = NA) +
    geom_jitter(aes(colour = Sex), width = 0.2) +
    facet_wrap(~name) +
    ylab(expression(volume ~ (mm ^ {3}))) +
    ggtitle("basal forebrain volume")
dev.off()  

png("putamen.png", width = 10, height = 5, units = "cm", res = 300)
qplot(group, volume, data = subset(regional, label %in% c(34, 35)),
      geom = "boxplot", outlier.shape = NA) +
  geom_jitter(aes(colour = Sex), width = 0.2) +
  facet_wrap(~name) +
  ylab(expression(volume ~ (mm ^ {3}))) +
  ggtitle("putamen volume")
dev.off()  


