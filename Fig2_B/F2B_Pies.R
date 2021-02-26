library(ggplot2)
library(ggpubr)

#Majority Clonal

DF <- read.csv('ORR__CLONAL.csv', header = TRUE, sep=",")
pdf('Pie_Majority_Clonal_PATIENT.pdf', height = 6, width = 12)
bp<- ggplot(DF, aes(x="", y = Prop, fill=ORR)) +
geom_bar(width = 1, stat = "identity") + theme_bw() + labs(x="", y="Majority Clonal (n = 16)") +
scale_fill_manual("ORR", values = c("CR" = "#30c180", "PD" = "#EB0FE0", "PR" = "#005ECC", "SD" = "#E1D611"))
pie <- bp + coord_polar("y", start=0)
pie
dev.off()

#Majority Subclonal

DF <- read.csv('ORR_SUBCLONAL.csv', header = TRUE, sep=",")
pdf('Pie_Majority_Subclonal_PATIENT.pdf', height = 6, width = 12)
bp<- ggplot(DF, aes(x="", y = Prop, fill=ORR)) +
geom_bar(width = 1, stat = "identity") + theme_bw() + labs(x="", y="Majority Subclonal (n = 11)") +
scale_fill_manual("ORR", values = c("CR" = "#30c180", "PD" = "#EB0FE0", "PR" = "#005ECC", "SD" = "#E1D611"))
pie <- bp + coord_polar("y", start=0)
pie
DF
dev.off()