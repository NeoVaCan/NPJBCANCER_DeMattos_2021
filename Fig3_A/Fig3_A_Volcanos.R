library("cowplot")
library("dplyr")
library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels

MyData <- read.csv(file="down_2fold.csv", sep = ',', header = TRUE)

volc_DN = ggplot(MyData, aes(foldChange, -log10(pval), color = timepoint)) + #volcanoplot with log2Foldchange versus pvalue
    geom_point() +
    labs(colour = "Time Point") +
    ggtitle("VolcanoPlotDN") + #name of te plot
    xlab("Log(2) FoldChange") #x axis label

MyPlot_DN = volc_DN + geom_text_repel(data=head(MyData, 5), aes(label=gene, size = 3)) + #adding text for the top 22 genes
    scale_color_manual(values=c("#FFCD00B2", 'chartreuse3',"#CC0C00B2")) + #seting color manually
    geom_hline(yintercept=2, linetype="dashed", color = "red") + #dashed horizontal line
    scale_y_continuous(limit = c(-0.05, 5)) #adjust axis scale
    MyData <- read.csv(file="up_2fold.csv", sep = ',', header = TRUE)

volc_UP = ggplot(MyData, aes(foldChange, -log10(pval), color = timepoint)) + #volcanoplot with log2Foldchange versus pvalue
    geom_point() +
    labs(colour = "Time Point") +
    ggtitle("VolcanoPlotUP") + #name of te plot
    xlab("Log(2) FoldChange") #x axis label

MyPlot_UP = volc_UP + geom_text_repel(data=head(MyData, 7), aes(label=gene, size = 3)) + #adding text for the top 22 genes
    scale_color_manual(values=c("#FFCD00B2", 'chartreuse3',"#CC0C00B2")) + #seting color manually
    geom_hline(yintercept=2, linetype="dashed", color = "red") + #dashed horizontal line
    scale_y_continuous(position ="right", limit = c(-0.05, 5)) #adjust axis scale

pdf("VolcanoPlot.pdf", width=10, height=4, useDingbats = FALSE )
plot(test01)
dev.off()