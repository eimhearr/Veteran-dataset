#load relevant packages
library(survival)
install.packages("survminer")
library(survminer)
library(dplyr)
library(tidyverse)
library(ggplot2)
install.packages("plotly")
library(plotly)

#load data into dataframe & examine data
df <- read.table("survival.txt", header=TRUE)
glimpse(df)
summary(df)

#make temp object of survival status as factor type to see number of events vs non events.
temp_surv <- as.factor(df$SurvStatus)
summary(temp_surv)

#Menopause status data type is incorrectly stored as integer data type so I'm
#going to coerce it to categorical type, I'm also going to recode the 1's back
#to 'Yes' and 2's back to 'No' per the source data so its human readable

df$MenopauseStatus <- gsub(df$MenopauseStatus, pattern="1", replacement = "No")
df$MenopauseStatus <- gsub(df$MenopauseStatus, pattern="2", replacement="Yes")
df$MenopauseStatus <- as.factor(df$MenopauseStatus)

#same step as above for Hormone therapy
df$HormoneTherapy <- gsub(df$HormoneTherapy, pattern="1", replacement="No")
df$HormoneTherapy <- gsub(df$HormoneTherapy, pattern="2", replacement ="Yes")
df$HormoneTherapy <- as.factor(df$HormoneTherapy)

#Coercing TumourGrade to ordinal data type
df$TumourGrade <- as.factor(df$TumourGrade)
levels(df$TumourGrade)

#Some exploratory data analysis
summary(df)
range(df$Age)
age_histogram <- ggplot(df, aes(x=Age)) +geom_histogram()
age_histogram
age_histogram_plotly <- ggplotly(p=age_histogram)
age_histogram_plotly
shapiro.test(df$Age)

#Draw survival curves. Use survdiff to examine p-value between groups.
fit1 <- survfit(Surv(SurvTime, SurvStatus) ~TumourGrade, data=df)
summary(fit1)
ggsurvplot(fit1, data=df, pval=TRUE, conf.int=FALSE,xlab="Time in days", ggtheme=theme_classic(),
           risk.table = TRUE, censor.mark=TRUE, title="Overall survival of patients stratified by tumour grade",
           surv.median.line="hv")

survdiff(formula = Surv(SurvTime, SurvStatus) ~ TumourGrade, data=df)
fit2 <- survfit(Surv(SurvTime, SurvStatus) ~MenopauseStatus, data=df)
summary(fit2)
ggsurvplot(fit2, data=df, pval=TRUE, conf.int=FALSE,xlab="Time in days",risk.table = TRUE,
           ggtheme=theme_classic(), title= "Overall survival of patients statified by menopausal status",
           censor.mark=TRUE, surv.median.line="hv")
survdiff(formula = Surv(SurvTime, SurvStatus) ~ MenopauseStatus, data=df)

fit3 <- survfit(Surv(SurvTime, SurvStatus) ~HormoneTherapy, data=df)
summary(fit3)
ggsurvplot(fit3, data=df, pval=TRUE, conf.int=FALSE,xlab="Time in days",title="Overall survival of patients stratified by Hormone therapy",
           ggtheme=theme_classic(),
           censor.mark=TRUE,risk.table = TRUE, surv.median.line="hv")
survdiff(formula = Surv(SurvTime, SurvStatus) ~ HormoneTherapy, data=df)

#Generate univariate CPH models for each factor in turn.
#Then test cox proportionality assumption for each factor using schoenfeld plots
#then test linearity of continuous variables using martingale residuals

grade_cox <- coxph(Surv(SurvTime,SurvStatus)~factor(TumourGrade),method="efron",data=df)
summary(grade_cox)
grade_coxph <- cox.zph(grade_cox, transform='log')
plot(grade_coxph)
plot(grade_coxph[1,])

age_cox <- coxph(Surv(SurvTime,SurvStatus)~Age,method="efron",data=df)
summary(age_cox)
age_coxph <- cox.zph(age_cox, transform='log')
age_coxph
plot(age_coxph[1,])

#alternative survminer method for making schoenfeld plots. More visually appealing.
ggcoxdiagnostics(age_cox,type="schoenfeld",ox.scale="time")
ggcoxdiagnostics(age_cox,type="martingale",ox.scale="linear.predictions")
#martingale residuals using survminer.

size_cox <- coxph(Surv(SurvTime,SurvStatus)~TumourSize,method="efron",data=df)
summary(size_cox)
size_coxph <- cox.zph(size_cox, transform='log')
size_coxph
plot(size_coxph[1,])
ggcoxdiagnostics(size_cox,type="schoenfeld",ox.scale="time")
ggcoxdiagnostics(size_cox,type="martingale",ox.scale="linear.predictions")#martingale residuals

menopause_cox <- coxph(Surv(SurvTime,SurvStatus)~MenopauseStatus,method="efron",data=df)
summary(menopause_cox)
ggcoxdiagnostics(menopause_cox,type="schoenfeld",ox.scale="time")

hormone_cox<- coxph(Surv(SurvTime,SurvStatus)~HormoneTherapy,method="efron",data=df)
summary(hormone_cox)
ggcoxdiagnostics(hormone_cox,type="schoenfeld",ox.scale="time")

#Multivariate model containing all variables
full_model <- coxph(Surv(SurvTime,SurvStatus)~Age+TumourGrade+TumourSize+HormoneTherapy
                    +MenopauseStatus,method="efron",data=df)
summary(full_model)
ggforest(full_model)

#Model selection using AIC
install.packages("MASS")
library(MASS)

all.coxph <-coxph(Surv(SurvTime,SurvStatus)~Age+HormoneTherapy+MenopauseStatus+TumourGrade
        +TumourSize, method="efron", data=df)
stepAIC(all.coxph)

#Final Multivariate Model. AIC suggests a two-factor model: TumourGrade & TumourSize
two_factor.coxph<-coxph(Surv(SurvTime,SurvStatus)~TumourGrade+TumourSize, method="efron",
                        data=df)
summary(two_factor.coxph)
adjusted <-ggadjustedcurves(two_factor.coxph,data=df, variable="TumourGrade", xlab="Time (days)")
ggpar(p = adjusted,legend.title="TumourGrade")

#Create tables and other statistics for report
install.packages("Tableone")
library(tableone)

df1 <- df
df1$SurvStatus <- as.factor(df1$SurvStatus)
MyVars <-c("Age","MenopauseStatus","HormoneTherapy","TumourSize","TumourGrade",
           "SurvTime","SurvStatus")
tab1 <- CreateTableOne(data =df1, vars=MyVars)
print(tab1, showAllLevels = TRUE)
summary(tab1)
tab2 <-CreateTableOne(vars=MyVars, strata="SurvStatus", data=df1)
print(tab2, showAllLevels = TRUE)
summary(tab2)
median.summary <- survfit(Surv(SurvTime,SurvStatus)~1, data=df)
median.summary
