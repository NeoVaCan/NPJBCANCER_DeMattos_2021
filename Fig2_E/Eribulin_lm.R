
library(ggplot2)
library(ggpubr)
library(dplyr)
library(car)
library(ggpmisc)

# Linear regression Nonsynonymous Mutation load and Mutations generating neoantigens 
mut.load=read.table(file="Sample_TILs_PAM50_ORR_Neoprep_Maf.txt", header = TRUE, sep = "\t", na.strings = "NA", dec =".")


# save predictions of the model in the new data frame 
fit <- lm(maf_nonsyn.mutations~Mutations, data=mut.load)
summary(fit)

pdf("Samples_nonsyn-mutload_mutneo_lm.pdf",useDingbats = FALSE, height = 6, width = 7)
ggplot(fit, aes(maf_nonsyn.mutations, Mutations)) +labs(x="Nonsynonymous Mutations load", y="Mutation generating neoantigens")+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE)+
  theme_classic()+
  theme(text = element_text(size=10), axis.title = element_text(size = 15), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12))+
  annotate("text", x = 50, y = 100, label = "paste(R^2 == 0.9453)", colour="black", size = 5, parse=TRUE)
dev.off()


# Linear regression TILs and Neoantigens per patient
tils=read.table(file="Patient_TILsMeans_PAM50_ORR_HR_status_Neoprep_maf.csv", header = TRUE, sep = ",", na.strings = "NA", dec =".")
tils

# save predictions of the model in the new data frame 
fit <- lm(TILs_Mean~Peptides, data=tils)
summary(fit)


pdf("Patient_TILsMean_Neoantigen_lm.pdf",useDingbats = FALSE, height = 6, width = 7)
ggplot(fit, aes(TILs_Mean, Peptides)) + labs(x="Median TILs per patient",y="Peptides")+
  geom_point() +xlim(0, 40)+ylim(0, 160)+
  geom_smooth(method = "lm", se = TRUE)+
  theme_classic()+
  theme(text = element_text(size=10), axis.title = element_text(size = 15), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12))+
  annotate("text", x = 20, y = 150, label = "paste(R^2 == 0.1605)", colour="black", size = 5, parse=TRUE)
dev.off()