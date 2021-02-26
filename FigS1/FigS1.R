clonality <- read.delim("clonality.txt")
outcome <- read.csv("outcome.csv",sep=";")
PAM50 <- read.csv("PAM50.txt",sep="\t")
aux <- read.csv("HR_data_neoerib.csv")
id_general <- read.csv("id_general.csv",sep=";")
drivers<- read.csv("Drivers.csv")
signatures <- read.csv("Signature.csv")

library(rms)
library(dplyr)
library(magrittr)
library(survival)

############################
##### Figure S1.A
############################

drivers$patient <- substr(drivers$Case.ID, 1,4)
drivers$patient <- as.factor(drivers$patient)
results <- data.frame(genes=colnames(drivers)[-c(1,ncol(drivers))])
results$n <- NA

for (j in levels(results$genes)){
  cont <- 0
  for (i in levels(drivers$patient)){
    val <- ifelse(sum(!is.na(drivers[drivers$patient==i,j]))==0,0,1)
    cont <- cont + val
  }
  results$n[results$genes==j] <- cont
}
results2<-results[order(results$n,decreasing=TRUE),]

mypal <-c('chartreuse3',"#FFCD00B2", "#CC0C00B2","#5C88DAB2", 'deeppink3',brewer.pal(n = 10, name = 'Set3')[4:10])
mypal <- mypal[c(1:8,10,9,12,11)]
par(mfrow=c(3,1),mar = c(2, 5, 3, 2))
windows(20,12)
layout(matrix(c(1,1,1,2), nrow = 4, ncol = 1, byrow = TRUE))
barplot(results2$n[1:11],las=2,names=results2$genes[1:11], ylim=c(0,20),cex.axis=2.5,cex=2.5,col=mypal,font.axis = 2, main="Most prevalent genes (n=28)", cex.main=3)

text(0.7,16,"15",font=1,cex=2.2)
text(1.9,12,"11",font=1,cex=2.2)
text(3.1,4,"3",font=1,cex=2.2)
text(4.3,3,"2",font=1,cex=2.2)
text(5.5,3,"2",font=1,cex=2.2)
text(6.7,3,"2",font=1,cex=2.2)
text(7.9,3,"2",font=1,cex=2.2)
text(9.1,3,"2",font=1,cex=2.2)
text(10.3,3,"2",font=1,cex=2.2)
text(11.5,3,"2",font=1,cex=2.2)
text(12.7,3,"2",font=1,cex=2.2)
text(11.5,3,"2",font=1,cex=2.2)

############################
####### Figure S1.B
############################

## Date format
outcome$DOB<-as.Date(outcome$DOB,"%d/%m/%Y")
outcome$BC_Diagnosis<-as.Date(outcome$BC_Diagnosis,"%d/%m/%Y")
outcome$Last_visit<-as.Date(outcome$Last_visit,"%d/%m/%Y")
outcome$Recurrence<-as.Date(outcome$Recurrence,"%d/%m/%Y")
outcome$Death <-as.Date(outcome$Death,"%d/%m/%Y")

## DFS estimation
outcome$dfs_time <- as.numeric(outcome$Recurrence-outcome$BC_Diagnosis)/365.25
outcome$dfs_status <- ifelse(!is.na(outcome$dfs_time),1,0)
outcome$dfs_time[is.na(outcome$dfs_time)]<- as.numeric(outcome$Last_visit[is.na(outcome$dfs_time)]-outcome$BC_Diagnosis[is.na(outcome$dfs_time)])/365.25
dfs <- Surv(outcome$dfs_time , outcome$dfs_status)

## Preprocessing
clonality$pts <- substr(clonality$id, 1, 4)
clonality$Gene<- as.character(clonality$Gene)
clonality$pts<-as.factor(clonality$pts)
clonality$visit <- substr(clonality$id, 6, 7)
clonality$patient <- substr(clonality$id, 1, 4)

#Creamos variable categorica para TMB segun low/high
mutations_MB$visit <- substr(mutations_MB$id_visit, 6, 7)
mutations_MB$patient <- substr(mutations_MB$id_visit, 1, 4)
mutations_V1 <- subset(mutations_MB, visit == "V1")
mutations_MB_V2 <- subset(mutations_MB, visit == "V2")
colnames(PAM50)[8] <- "id_visit"
mutations_MB_Basal<- merge(mutations_MB_V2,PAM50 ,by="id_visit", all.x = TRUE)
mutations_MB_V2_Basal <- subset(mutations_MB_Basal, mutations_MB_Basal$Subtype.Call == "BasalLike")

colnames(mutations_MB_V2_Basal)[7] <- "patient"
outcome_MB <- merge(outcome,mutations_MB_V2_Basal ,by="patient", all.x = TRUE)
mut_v1 <-subset(clonality,visit=="V1" & Gene=="TP53")
final0 <- merge(x = outcome_MB, y = mut_v1, by = "patient", all.x = TRUE)
PAM50_1 <- subset(PAM50, Visit=="V1")
final <- merge(x = final0, y = PAM50_1, by = "patient", all.x = TRUE)
final$TP53 <- ifelse(is.na(final$id),0,1)

### Kaplan-meier
mypal <- c(brewer.pal(5,"Blues")[4],brewer.pal(7,"YlOrRd")[c(4,7)])
surv2 <- survfit(Surv(final$dfs_time,final$dfs_status) ~ final$TP53)


windows(20,16)
plot( 1,1,col="white")
ggsurvplot(surv2, data = final,
           title = "",
           pval = F, pval.method = F,    # Add p-value &  method name
           palette = mypal,                   # Use JCO journal color palette o AAAS palette, D3 etc.
           risk.table = T,                  # Add No at risk table
           cumevents = F,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = F,
           conf.int = F, # Hide tables y axis text
           xlab= "Time (years)",
           ylab="Proportion of event-free (RFS)",
           font.tickslab=18,
           pval.size=4.5,
           risk.table.title="Number at risk",
           risk.table.fontsize=5.5,
           font.y=c(18,"bold"),
           size=1.8,
           font.x=c(18,"bold"),
           linetype=c(1,1),
           legend=c(1.5,0.3),
           legend.title = "",           # Change legend titles
           legend.labs =  c(" ","  "),  # Change legend labels
           font.legend=c(18,"bold"),
           break.time.by=1,
           xlim=c(0,6))

### Kaplan-Meier 2
mut_v1 <-subset(clonality,visit=="V2" & Gene=="PIK3CA")
final <- merge(x = outcome_MB, y = mut_v1, by = "patient", all.x = TRUE)
final$PIK <- ifelse(is.na(final$id),0,1)
surv3 <- survfit(Surv(final$dfs_time,final$dfs_status) ~ final$PIK)

windows(20,16)
plot( 1,1,col="white")
ggsurvplot(surv3, data = final,
           title = "",
           pval = F, pval.method = F,    # Add p-value &  method name
           palette = mypal,                   # Use JCO journal color palette o AAAS palette, D3 etc.
           risk.table = T,                  # Add No at risk table
           cumevents = F,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = F,
           conf.int = F, # Hide tables y axis text
           xlab= "Time (years)",
           ylab="Proportion of event-free (RFS)",
           font.tickslab=18,
           pval.size=4.5,
           risk.table.title="Number at risk",
           risk.table.fontsize=5.5,
           font.y=c(18,"bold"),
           size=1.8,
           font.x=c(18,"bold"),
           linetype=c(1,1),
           legend=c(1.5,0.3),
           legend.title = "",           # Change legend titles
           legend.labs =  c(" ","  "),  # Change legend labels
           font.legend=c(18,"bold"),
           break.time.by=1,
           xlim=c(0,6))


##############################
######## Figure S1.C
##############################

signatures$patient <- substr(signatures$Column1, 1,4)
signatures$APOBEC <- signatures$APOBEC.2. + signatures$Signature.13
names(signatures)

sign_group <- signatures %>%
  group_by(patient) %>%
  summarise(APOBEC = sum(APOBEC),BRCA= sum(BRCA.3.), DNA.repair= sum(Def_DNA_rep.6.), Poly=sum(Polymerase.9.),Pole=sum(POLE.mutations..10.))

final <- merge(outcome,sign_group ,by="patient")

final$APOBEC <- ifelse(final$APOBEC>0,1,0)
table(final$APOBEC, final$Response)  
apobec <- prop.table(table(final$APOBEC, final$Response),1)*100
fisher.test(table(final$APOBEC, final$Response)) 

final$BRCA <- ifelse(final$BRCA>0,1,0)
table(final$BRCA, final$Response) 
brca <-prop.table(table(final$BRCA, final$Response),1)*100
fisher.test(table(final$BRCA, final$Response))  

final$DNA.repair <- ifelse(final$DNA.repair>0,1,0)
table(final$DNA.repair, final$Response)  
dna_reapir <- prop.table(table(final$DNA.repair, final$Response),1)*100
fisher.test(table(final$DNA.repair, final$Response))  

final$Poly <- ifelse(final$Poly>0,1,0)
table(final$Poly, final$Response)  
final$Pole <- ifelse(final$Pole>0,1,0)
table(final$Pole, final$Response) 
pole <- prop.table(table(final$Pole, final$Response),1 )*100

signature <- c(pole[2,2],dna_reapir[2,2],brca[2,2],apobec[2,2])
no_signature <- -c(pole[1,2],dna_reapir[1,2],brca[1,2],apobec[1,2])
windows(18,8)
barplot(signature,horiz = T, col = "darkblue",xlim=c(-100,70), xlab="Proportion of responses (%)",cex.lab=2)
barplot(no_signature,add = TRUE,horiz = T, axes = F, col = "darkorange")
