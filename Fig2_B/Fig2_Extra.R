outcome <- read.csv("outcome.csv",sep=";")

library(rms)
library(dplyr)
library(magrittr)
library(survival)

############################
##### Figure 2B.2
############################

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

mutations2 <- clonality %>%
  group_by(id) %>%
  summarise (n = n(),nclonal = sum(Clonal_Subclonal_mutation=="Clonal"),nno = sum(Clonal_Subclonal_mutation!="Clonal")) %>% 
  ungroup();mutations2 

mutations <- as.data.frame(mutations2)
names(mutations)[1]<-"id_visit"
mutations_MB <- merge(id_general,mutations ,by="id_visit", all.x = TRUE)

mutations_MB$visit <- substr(mutations_MB$id_visit, 6, 7)
mutations_MB$patient <- substr(mutations_MB$id_visit, 1, 4)
mutations_V1 <- subset(mutations_MB, visit == "V1")
ff <- merge(outcome, mutations_V1 ,by="patient", all.x = TRUE)
ff$subclonal <- ff$nno/ff$n
mod1 <- glm(Response ~ subclonal, data = ff, family = "binomial");summary(mod1)
dat1<-predict(mod1, ff, se.fit=TRUE)
ff$estimation <- exp(dat1$fit)/(1+exp(dat1$fit))
ff2 <- merge(ff, aux ,by="patient", all.x = TRUE)

library(rms)
ff2 <- subset(ff2, !is.na(estimation))
ff2$HR_status <- droplevels(ff2$HR_status)
model1 <- cph(Surv(dfs_time,dfs_status)~  rcs(subclonal,3) , data=ff2,x=T,y=T,surv=TRUE);model1
model2  <- lrm(Response ~  rcs(subclonal,3), data=ff2,x=T,y=T)

dd <- datadist(ff2);options(datadist="dd") 
windows(20,20)
plot(Predict(model2, subclonal),xlim=c(0,0.6), ylim=c(-3,3))
windows(20,20)
plot(Predict(model1, subclonal),xlim=c(0,0.6), ylim=c(-4,3))

model3 <- cph(Surv(dfs_time,dfs_status)~  rcs(subclonal,3) + HR_status, data=ff2,x=T,y=T,surv=TRUE);model1
model4  <- lrm(Response ~  rcs(subclonal,3)+ HR_status, data=ff2,x=T,y=T)
dd <- datadist(ff2);options(datadist="dd") 
windows(20,20)
plot(Predict(model4, subclonal,HR_status,conf.int=FALSE),xlim=c(0,0.6), ylim=c(-3,3))
windows(20,20)
plot(Predict(model3, subclonal,HR_status,conf.int=FALSE),xlim=c(0,0.6), ylim=c(-4,3))