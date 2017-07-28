library(ggplot2)
library(plotrix)

source("/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/Code/RBMperfrs_20170712.r")
setwd("/home/Pplateforme/Plate-forme_RMN/NachiketNadkarni/DCLK3/processed_20170727")

atlasnames = list.files(pattern = "atlas_Na.nii.gz", recursive=T)
RBMdata = llply(atlasnames, readRBMdata, .progress = "text")
RBMdata = do.call(rbind, RBMdata)
regions = read.table("/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/c57_brain_atlas_labels.txt",
                     header = T, sep = "\t"
                     )
RBMdata$region = RBMdata$label
RBMdata$side = RBMdata$label
for (n in 1:nrow(regions)) {
  label = regions$label[n]
  RBMdata$region[which(RBMdata$label == label)] = as.character(regions$structure[n])
  RBMdata$side[which(RBMdata$label == label)] = as.character(regions$side[n])
}
RBMdata = subset(RBMdata, side == "left" | side == "right" | side == "both") #gets rid of non-existent VOInumbers

RBMdataleftrightboth = ddply(.data = subset(RBMdata, side != "both"),
                             c("fname", "region"),
                             summarise,
                             count=sum(count)
                             )
RBMdataleftrightboth$label = NA
RBMdataleftrightboth$side = "both"
RBMdata = rbind(RBMdata, RBMdataleftrightboth)

RBMdatatotalboth = ddply(.data=subset(RBMdata, side == "both"),
                         "fname",
                         summarise,
                         count=sum(count)
                         )
RBMdatatotalboth$label = NA
RBMdatatotalboth$region = "Total"
RBMdatatotalboth$side = "both"
RBMdatatotalleft = ddply(.data=subset(RBMdata, side == "left"),
                         c("fname", "side"),
                         summarise,
                         count=sum(count)
                         )
RBMdatatotalleft$label = NA
RBMdatatotalleft$region = "Total"
RBMdatatotalright = ddply(.data = subset(RBMdata, side == "right"),
                          c("fname", "side"),
                          summarise,
                          count=sum(count)
                          )
RBMdatatotalright$label = NA
RBMdatatotalright$region = "Total"

RBMdata = rbind(RBMdata, RBMdatatotalboth, RBMdatatotalleft, RBMdatatotalright)
RBMdata$volume = RBMdata$count * 0.1 *0.1 *0.1
Genotypes=data.frame(Genotype = c("wt", "DCLK3", "wt", "DCLK3", "DCLK3", "wt", "DCLK3", "wt", "wt", "DCLK3", "DCLK3", "DCLK3", "wt", "wt", "wt", "wt", "DCLK3", "DCLK3"),
                     Scan_Mon = c("Feb", "Feb", "Feb", "Feb", "Feb", "Feb", "Jul",  "Jul", "Jul", "Jul", "Jul", "Jul", "Jul", "Jul", "Jul", "Jul", "Jul", "Jul"),
                     fname=atlasnames
                     )
Genotypes$Genotype=factor(Genotypes$Genotype,levels=c("wt","DCLK3"))
RBMdata=merge(RBMdata,Genotypes)

graph1 = qplot(Genotype, volume, data = subset(RBMdata, side == "both"),
               geom = "boxplot", outlier.shape = NA, ymin = 0) + 
         geom_jitter(aes(colour = fname, shape = Scan_Mon), width = 0.25) +
         facet_wrap(~region, scales = "free") + ylab(expression(volume ~ (mm ^ {3})))

graph2 = qplot(Genotype, volume, data = subset(RBMdata, side == "both"),
               geom = "boxplot", outlier.shape = NA) + 
         geom_jitter(aes(colour = fname, shape = Scan_Mon), width = 0.25) +
         facet_wrap(~region, scales = "free") + ylab(expression(volume ~ (mm ^ {3})))

pdf(file="RBM.pdf",width=22.5*sqrt(2),height=22.5)
  graph1
  graph2
dev.off()

png(filename="RBM.png",width=6000*sqrt(2),height=6000, res=300)
  graph2
dev.off()

wilcox.test(volume ~ Genotype, data = subset(RBMdata, side == "both" & region == "hippocampus"))
wilcox.test(volume ~ Genotype, data = subset(RBMdata, side == "both" & region == "cerebral cortex: entorhinal cortex"))
wilcox.test(volume ~ Genotype, data = subset(RBMdata, side == "both" & region == "midbrain"))
wilcox.test(volume ~ Genotype, data = subset(RBMdata, side == "both" & region == "cerebellar cortex"))

write.table(RBMdatavoltab <- ddply(.data = subset(RBMdata,
                                                  region == "Total" & side=="both"
                                                  ),
                                   c("Genotype"), 
                                   summarise,
                                   N = length(unique(fname)),
                                   mean = mean(volume),
                                   sd = sd(volume),
                                   se = std.error(volume),
                                   median = median(volume),
                                   halfIQR = 0.5 * IQR(volume)
                                   ),
            file = "RBMTotaltab.txt",
            quote = F,
            sep = "\t",
            row.names = F
            )
write.table(RBMdata,
            file = "RBMdata.txt",
            quote = F,
            sep = "\t",
            row.names = F
            )

