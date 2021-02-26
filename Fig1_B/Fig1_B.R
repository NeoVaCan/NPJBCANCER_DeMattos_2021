
library(rms)
library(dplyr)
library(magrittr)
library(survival)

PAM50 <- read.csv("PAM50.txt",sep="\t")
outcome <- read.csv("outcome:s.csv",sep=";")

v1 <- subset(PAM50,Visit=="V1")
v1$Subtype.Call <- as.character(v1$Subtype.Call)
v1$Subtype.Call[v1$Subtype.Call==""] <- "NA"
windows(20,20)
pie(round(prop.table(table(v1$Subtype.Call))*100,1),col=c("red","grey",brewer.pal(8,"Blues")[c(7,8)],"grey"))
windows(20,20)
pie(round(prop.table(table(outcome$ORR))*100,1),col=c('darkgreen',"#FFCD00B2", "#CC0C00B2","darkred"))
