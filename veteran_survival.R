#load relevant packages

library(survival)
library(survminer)
install.packages("ranger")
library(ranger)
library(ggplot2)
library(dplyr)
install.packages("ggfortify")
library(ggfortify)


#load data

data(veteran)
str(veteran)
glimpse(veteran)
summary(veteran)

#EDA

veteran$trt <- factor(veteran$trt, labels = c("standard","test"))
veteran$prior <- factor(veteran$prior, labels = c("no","yes"))

vet_temp <- veteran
vet_temp$status <- as.factor(vet_temp$status)
summary(vet_temp)
rm(vet_temp)

#data visualisation

ggplot(data=veteran)+
  geom_histogram(binwidth = 100, mapping=aes( x=time, fill=trt))+
  facet_grid(trt ~.)+
  ggtitle("Figure 1")

ggplot(data=veteran)+
  geom_histogram(binwidth=5, mapping=aes(x=karno,fill=celltype))+
  ggtitle("Fig.2")


 fit <- survfit(Surv(time, status) ~celltype, data=veteran)
summary(fit)  
ggsurvplot(fit, pval=TRUE,conf.int=TRUE,xlab="Time in days", ggtheme=theme_classic(),
           censor.mark=TRUE, surv.median.line="hv",cumevents = TRUE, risk.table = TRUE,
           fun = "cumhaz", )
fit1 <- survfit(Surv(time,status)~prior, data=veteran)
summary(fit1)
ggsurvplot(fit1, pval=TRUE,conf.int=TRUE,xlab="Time in days", ggtheme=theme_classic(),
           censor.mark=TRUE, surv.median.line="hv")


fit_cox_uni <- coxph(Surv(time, status) ~ , data = veteran)
summary(fit_cox_uni)
fit_cox_multi <- coxph(Surv(time, status) ~ ., data = veteran)
summary(fit_cox_multi)
