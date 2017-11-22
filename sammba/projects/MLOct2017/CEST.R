library(ggplot2)
library(readxl)

CESTdir = "/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/MLOct2017/CEST"

youngdir = paste(CESTdir, "young", sep = "/")
olddir = paste(CESTdir, "old", sep = "/")

CESTROIs = function(CESTfile, mask){
  afnitbl = read.delim(text = system2("/usr/lib/afni/bin/3dROIstats",
                                      c("-mask", mask, "-nzmean", "-nzvoxels", CESTfile),
                                      stdout = TRUE, stderr = FALSE))
  afnitbl
}

#mask = "../../Lemur-Atlas-Apr2Feb-cortexRL-label2_200.nii.gz"
mask = paste(CESTdir, "dilated.nii.gz", sep = "/")

old = lapply(list.files(olddir, full.names = T), CESTROIs, mask)
young = lapply(list.files(youngdir, full.names = T), CESTROIs, mask)
old = do.call(rbind, old)
young = do.call(rbind, young)
old$group = "old"
young$group = "young"
x = rbind(old, young)
x$group = factor(x$group, levels = c("young", "old"))
z = reshape(x, varying = names(x)[grep('_', names(x))], timevar = "intensity",
            direction = "long", sep= "_")

projectdir = "/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/Tri_data"

IDs = read_xlsx(paste(projectdir, "Trilemur_ok.xlsx", sep = "/"), sheet = 1)
MRIsessions = read_xlsx(paste(projectdir, "Trilemur_ok.xlsx", sep = "/"), sheet = 2)
all = merge(IDs, MRIsessions)
all$age = all$date - all$DOB
all$MRIsession = basename(all$PVdir)
z$MRIsession = sapply(strsplit(basename(as.character(z$File)), "__"), "[[", 1)
newall = merge(all, z)
mininewall = subset(newall, intensity == 1)

#sanity check
qplot(group, age, data = mininewall)

qplot(age, Mean, data = mininewall)

png(paste(CESTdir, "gluCEST.png", sep = "/"), width = 10, height = 10, units = "cm", res = 300)
qplot(group, Mean * 100, data = mininewall, geom = "boxplot",
      outlier.shape=NA) + geom_jitter(aes(colour = Sex), width = 0.2) +
  ylab("gluCEST (%)")
dev.off()


t.test(Mean ~ group, data=mininewall)
wilcox.test(Mean ~ group, data = mininewall)
summary(lm(Mean ~ age, data = mininewall))

a = qplot(group, Mean, data = newall, geom = "boxplot", outlier.shape=NA) +
      geom_jitter(width = 0.2) +
      facet_wrap(~intensity, scales = "free_y")

pdf(file = "CEST.pdf", width = 15* sqrt(2), height = 15)
a
dev.off()

func1 = qplot(age, data = subset(all, date < "2017-06-01"),
              geom = "histogram", fill = Sex, xmin = 0, binwidth = 365.23/12) +
  scale_x_continuous(breaks = seq(0, 2000, 365.25),
                     minor_breaks = seq(0, 2000, 365.25/12),
                     labels = 0:5) +
  scale_y_continuous(breaks = 0:5, minor_breaks = NULL) +
  xlab("age (years)")

png("func1.png", width = 11, height = 7, units = "cm", res = 300)
func1
dev.off()


func2 = qplot(age, data = subset(all, date > "2017-06-01"),
              geom = "histogram", fill = Sex, xmin = 0, binwidth = 365.23/12) +
  scale_x_continuous(breaks = seq(0, 4100, 365.25),
                     minor_breaks = seq(0, 4100, 365.25/12),
                     labels = 0:11) +
  scale_y_continuous(breaks = 0:11, minor_breaks = NULL) +
  xlab("age (years)")

png("func2.png", width = 11, height = 7, units = "cm", res = 300)
func2
dev.off()


MRIsessions = read_xlsx("Liste-Microc?bes-Atlas.xlsx", sheet = 1)
MRIsessions$age = MRIsessions$`Date IRM` - MRIsessions$`Date de naissance`
MRIsessions$Gender = MRIsessions$"Sexe\r\n(? v?rifier)"

anathist = qplot(age, data = MRIsessions, geom = "histogram", fill = Gender,
                 xmin = 0, binwidth = 365.23/12) +
             scale_x_continuous(breaks = seq(0, 2000, 365.25),
                                minor_breaks = seq(0, 2000, 365.25/12),
                                labels = 0:5) +
             scale_y_continuous(breaks = 0:5, minor_breaks = NULL) +
             xlab("age (years)")

png("atlashist.png", width = 11, height = 7, units = "cm", res = 300)
anathist
dev.off()






