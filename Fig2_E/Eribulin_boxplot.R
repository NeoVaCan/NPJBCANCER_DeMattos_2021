library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(ggthemes)
library(lemon)


IO=read.table(file="log10_transformed_normalized_counts_trans.txt", header = TRUE, sep = "\t")

table(IO$Visit)
####ORR all significant
IO

IOboxplot <- IO[,c("Visit", "ORR", "ALDH1A1", "EVI2A", "PTGER4","LHFP", "ADRA2A", "VEGFA", "FANCA", "ORC6", "KIFC1", "ANGPTL4")]

x.m<-melt(IOboxplot, id.var=c("Visit", "ORR"))
x.m

pdf("orr_eribulina_boxplot.pdf", useDingbats = FALSE, height = 5, width = 12)
ggplot(x.m, aes(x = Visit, y = value, fill=ORR)) + ylab("Log10 transformed gene expression")+
  geom_boxplot(lwd = 0.2,  outlier.size=0.5)+
  scale_x_discrete(limits=c(
    "V1",
    "V2",
    "V3"),expand = c(0, 0.1))+
  scale_fill_manual(values = c(
    'Responders' = 'blue',
    'Non-responders' = 'red')) +
  stat_compare_means(size=4,label.y = 5, method = "wilcox.test", label = "p.signif", hide.ns = TRUE)+
  theme_classic()+facet_rep_wrap(~factor(x.m$variable), repeat.tick.labels = FALSE,  ncol= 5)+
  theme(legend.position="bottom", legend.text= element_text(size=11), legend.title=element_blank(), axis.title.x=element_blank(), text = element_text(size=13), strip.text = element_text(size=11), axis.text.y=element_text( size = 10))
dev.off()

###Recurrence all significant
IO

IOboxplot2 <- IO[,c("Visit", "Recurrence", "FABP5", "YBX1", "TUBB6","PLOD1", "CXCL8", "STC2", "TSPAN13", "SCUBE2", "MAGEA1", "ABCC8")]

x.m<-melt(IOboxplot2, id.var=c("Visit", "Recurrence"))
x.m


pdf("recurrence_eribulina_boxplot.pdf", useDingbats = FALSE, height = 5, width = 12)
ggplot(x.m, aes(x = Visit, y = value, fill=Recurrence)) + ylab("Log10 transformed gene expression")+
  geom_boxplot(lwd = 0.2,  outlier.size=0.5)+
  scale_x_discrete(limits=c(
    "V1",
    "V2",
    "V3"),expand = c(0, 0.1))+
  scale_fill_manual(values = c(
    'Non-recurrence' = 'grey',
    'Recurrence' = 'purple')) +
  stat_compare_means(size=4,label.y = 5, method = "wilcox.test",label = "p.signif", hide.ns = TRUE)+
  theme_classic()+facet_rep_wrap(~factor(x.m$variable), repeat.tick.labels = FALSE,  ncol= 5)+
  theme(legend.position="bottom", legend.text= element_text(size=11), legend.title=element_blank(), axis.title.x=element_blank(), text = element_text(size=13), strip.text = element_text(size=11), axis.text.y=element_text( size = 10))
dev.off()