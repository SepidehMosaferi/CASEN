# ============================================================================================
# Name: R code for the CASEN 2009  
# Author: Sepideh Mosaferi
# Date: Spring 2016
# ============================================================================================

require(stats); 
require(graphics)
library(plyr)
library(MASS)
library(lattice)
library(memisc)
library(stargazer)
library(foreign)
library(OpenMx)
library(PracTools)
library(sampling)
library(samplingbook)
library(survey)
library(reshape)
library(nlme)
library(ICC)
require(Rcpp)
library(lme4)
require(faraway)
require(psych)


casendata2009 <- read.spss("Q://Casen 2009//Data 2009//Casen2009Spss.sav")
casendata2009 <- as.data.frame(casendata2009)


# sort dataset in a hierarchical geographical area
SORTDATA <- casendata2009[order(casendata2009$REGION,casendata2009$COMUNA,
                                casendata2009$ESTRATO,casendata2009$SEGMENTO),]

SORTDATA <- as.data.frame(SORTDATA)

#Subsetting the Missing Part of dataset based on the YMONEHAJ
MISS <- subset(SORTDATA,is.na(SORTDATA$YMONEHAJ),select=c(SEGMENTO:ESTRATO,EDAD,SEXO,YMONEHAJ,
                                                          YMONEHAJ,YAIMHAJ))  

length(unique(MISS$FOLIO))

ASSMISS <- sapply(1:nrow(MISS),function(i){SORTDATA$YMONEHAJ[SORTDATA$FOLIO==MISS$FOLIO[i]]})  #Assessing Missing Situation

#Hint: We do not have missing values based on the YMONEHAJ v.

SUBSTITUTE <- sapply(1:length(ASSMISS),function(i){ASSMISS[[i]][is.na(ASSMISS[[i]])] <- ASSMISS[[i]][1]})

SORTDATA$YMONEHAJ[is.na(SORTDATA$YMONEHAJ)] <- SUBSTITUTE

SORTDATA$HUID <- paste(SORTDATA$ESTRATO,SORTDATA$SEGMENTO,SORTDATA$IDVIV,sep="") ##unique housingunit identification

SORDATAnew <- SORTDATA[order(SORTDATA$HUID),]

HHFRAME <- SORDATAnew[!duplicated(SORDATAnew$FOLIO),]

HUFRAME <- HHFRAME[!duplicated(HHFRAME$HUID),]

HUFRAME$Zvariable <- aggregate(HHFRAME$YMONEHAJ~HHFRAME$REGION+HHFRAME$COMUNA+    
                                 HHFRAME$ESTRATO+HHFRAME$SEGMENTO+HHFRAME$HUID,FUN=sum)[,6]

HUFRAME$Hvariable <- aggregate(HHFRAME$FOLIO~HHFRAME$REGION+HHFRAME$COMUNA+     ##household in housing units
                                 HHFRAME$ESTRATO+HHFRAME$SEGMENTO+HHFRAME$HUID,FUN=function(x){length(x)})[,6] 


#investigating the number of PSU in strata:
PSUinSTRAT <- tapply(HUFRAME$SEGMENTO,HUFRAME$ESTRATO,FUN=function(x){length(unique(x))}) 
length(PSUinSTRAT[PSUinSTRAT==1]) #7 strata out of 602 strata contains only 1 PSU

######Variance Estimation and Estimators ######

##Based on Ultimate Cluster variance estimator (Survey Package)

#National Level:
options(survey.lonely.psu="adjust")
onestage<-svydesign(id=~SEGMENTO, strata=~ESTRATO, weight=~EXPR, data=HUFRAME, nest=TRUE)
National <- svyratio(HUFRAME$Zvariable,HUFRAME$Hvariable,onestage)
VarDesignNational <- unlist(National[2])

#ESTRATO level:
Uniq_ESTRATO <- unique(HUFRAME$ESTRATO)
RatioStrata <- VarDesignStrata <- rep(0,602)  

for(i in 1:602){
  SUBSTRAT <- subset(onestage,HUFRAME$ESTRATO==Uniq_ESTRATO[i])
  SUBSTRATRATIO <- svyratio(HUFRAME$Zvariable[HUFRAME$ESTRATO==Uniq_ESTRATO[i]],
                            HUFRAME$Hvariable[HUFRAME$ESTRATO==Uniq_ESTRATO[i]],SUBSTRAT)
  RatioStrata[i] <- SUBSTRATRATIO[1]
  VarDesignStrata[i] <- SUBSTRATRATIO[2]
}
RatioStrata; VarDesignStrata
RatioStrata <- unlist(RatioStrata)
VarDesignStrata <- unlist(VarDesignStrata)


#REGION and ZONA level:
Region <- unique(HUFRAME$REGION)
Zona <- unique(HUFRAME$ZONA)

RatioRegionZona <- VarDesignRegionZona <- matrix(0,15,2)

for(i in 1:15){
  for(j in 1:2){
    SUBRegionZona <- HUFRAME[which(HUFRAME$REGION==Region[i] & HUFRAME$ZONA==Zona[j]),]
    DESIGNRegionZona <- svydesign(id=~SEGMENTO, strata=~ESTRATO, weight=~EXPR, 
                                  data=SUBRegionZona, nest=TRUE)
    
    SVYRatioRegionZona <- svyratio(SUBRegionZona$Zvariable,SUBRegionZona$Hvariable,
                                   DESIGNRegionZona)
    RatioRegionZona[i,j] <- unlist(SVYRatioRegionZona[1])
    VarDesignRegionZona[i,j] <- unlist(SVYRatioRegionZona[2])
  }
}
RatioRegionZona
VarDesignRegionZona


#COMUNA level
COMUNA_NAME <- as.factor(levels(HUFRAME$COMUNA))
NACOMUNA <- tapply(HUFRAME$Zvariable,HUFRAME$COMUNA,FUN=function(x){length(x)})
REMCOMUNA <- names(which(is.na(NACOMUNA)))
FINALCOMUNA <- COMUNA_NAME[! COMUNA_NAME %in% REMCOMUNA]

RatioCOMUNA <- VarDesignCOMUNA <- rep(0,length(FINALCOMUNA))
for(i in 1:length(FINALCOMUNA)){
  options(survey.lonely.psu="adjust")
  DESIGN <- svydesign(id=~SEGMENTO, 
                      strata=~ESTRATO, 
                      weight=~EXPC, 
                      data=subset(HUFRAME,HUFRAME$COMUNA==FINALCOMUNA[i]), nest=TRUE)
  SVY_COMUNA <- svyratio(HUFRAME$Zvariable[HUFRAME$COMUNA==FINALCOMUNA[i]],
                         HUFRAME$Hvariable[HUFRAME$COMUNA==FINALCOMUNA[i]],DESIGN)
  RatioCOMUNA[i] <- unlist(SVY_COMUNA[1])
  VarDesignCOMUNA[i] <- unlist(SVY_COMUNA[2])
}
RatioCOMUNA; VarDesignCOMUNA

# number of PSU in the interaction comuna and zona
PSUINCOMUNAZONA <- aggregate(HUFRAME$SEGMENTO~HUFRAME$COMUNA+HUFRAME$ZONA,
                             FUN=function(x){length(unique(x))})[,3]

#The distribution of PSU in comuna
PSUinCOMUNA <- tapply(HUFRAME$SEGMENTO,HUFRAME$COMUNA,FUN=function(x){length(unique(x))})

cbind(as.vector(COMUNA_NAME),as.vector(PSUinCOMUNA))
boxplot(PSUinCOMUNA,ylim=c(0,25),main="The Distribution of Number of PSUs within
        \n the COMUNA (YMONEHAJ v.)",ylab="Number of SEGMENTO (PSUs)")


#Assessing the Small Area Estimation based on the number of PSUs in COMUNA:
PSUinCOMUNA_woNA <- as.vector(PSUinCOMUNA[!is.na(PSUinCOMUNA)])
STR <- cut(PSUinCOMUNA_woNA,breaks=c(-Inf,quantile(PSUinCOMUNA_woNA,
                                                   probs=c(0.25,0.5,0.75)),Inf)) 

min(PSUinCOMUNA_woNA); max(PSUinCOMUNA_woNA)
levels(STR) <- c("[2,9] \n Q1","(9,11] \n Q2","(11,15] \n Q3","(15,24] \n Q4")
L <- c(levels(STR))
COMUNAFRAME <- data.frame(FINALCOMUNA,VarDesignCOMUNA,PSUinCOMUNA_woNA,STR)
UPCOMUNAFRAME <- COMUNAFRAME[order(STR),]

table(UPCOMUNAFRAME$STR)
boxplot(VarDesignCOMUNA~STR,main="Variation of Variance among COMUNA
        ",xlab="Number of PSUs in COMUNA (YMONEHAJ v.)",
        ylab="Variance per COMUNA",col="gray80",
        font.main=3, cex.main=1.2, font.lab=3)

par(mfrow=c(1,4))
plot(UPCOMUNAFRAME$VarDesignCOMUNA[UPCOMUNAFRAME$STR==L[1]],ylab="",col="blue",pch=17,cex=1,
     ylim=c(0,max(VarDesignCOMUNA)),main="[2,9]")
plot(UPCOMUNAFRAME$VarDesignCOMUNA[UPCOMUNAFRAME$STR==L[2]],ylab="",col="red",pch=17,cex=1,
     ylim=c(0,max(VarDesignCOMUNA)),main="(9,11]")
plot(UPCOMUNAFRAME$VarDesignCOMUNA[UPCOMUNAFRAME$STR==L[3]],ylab="",col="green",pch=17,cex=1,
     ylim=c(0,max(VarDesignCOMUNA)),main="(11,15]")
plot(UPCOMUNAFRAME$VarDesignCOMUNA[UPCOMUNAFRAME$STR==L[4]],ylab="",col="sienna",pch=17,cex=1,
     ylim=c(0,max(VarDesignCOMUNA)),main="(15,24]")

COMPAR_PSUvsVAR <- sapply(1:length(levels(UPCOMUNAFRAME$STR)),
                          function(i){summary(UPCOMUNAFRAME$VarDesignCOMUNA[UPCOMUNAFRAME$STR==levels(UPCOMUNAFRAME$STR)[i]])})

colnames(COMPAR_PSUvsVAR) <- c("GROUP1","GROUP2","GROUP3","GROUP4")

##Formula of variance based on the model National Level

Ybarnation <- sum(HHFRAME$EXPR*HHFRAME$YMONEHAJ)/sum(HHFRAME$EXPR)

r <- length(HHFRAME$FOLIO)

SIGMA2 <- (sum(HHFRAME$EXPR*(HHFRAME$YMONEHAJ-Ybarnation)^2)/(sum(HHFRAME$EXPR)-1))

Wtilde <- HHFRAME$EXPR*(r/sum(HHFRAME$EXPR))

CV2wtilde <- (sum(Wtilde^2)/r)-1
bbar <- length(HHFRAME$FOLIO)/length(unique(HHFRAME$SEGMENTO))
#Intra Class Correlation
MixModel <- lmer(YMONEHAJ ~ 1 + (1 |SEGMENTO), HHFRAME)
MSB <- as.vector(sqrt(unlist(VarCorr(MixModel))))^2
MSE <- attr(VarCorr(MixModel), "sc")^2
Rho <- MSB/(MSE+MSB)
VarModelNational <- (SIGMA2/r)*(1+CV2wtilde)*(1+((bbar-1)*Rho)); 
VarModelNational
VarDesignNational
VarModelNational/VarDesignNational


##Formula of variance based on the model within ESTRATO
Uniq_ESTRATO <- unique(HHFRAME$ESTRATO)
SmallESTRATO <- as.numeric(names(PSUinSTRAT[PSUinSTRAT==1]))
Uniq_ESTRATOBIG <- Uniq_ESTRATO[!is.element(Uniq_ESTRATO, SmallESTRATO)]

VarModelStrata <- YbarModelStrata <-rep(0,length(Uniq_ESTRATOBIG))

for(i in 1:length(Uniq_ESTRATOBIG)){
  SUBESTRATO <- subset(HHFRAME,HHFRAME$ESTRATO==Uniq_ESTRATOBIG[i])
  
  ybar <- sum(SUBESTRATO$EXPR*SUBESTRATO$YMONEHAJ)/sum(SUBESTRATO$EXPR)
  
  r <- length(SUBESTRATO$FOLIO)
  
  SIGMA2 <- (sum(SUBESTRATO$EXPR*(SUBESTRATO$YMONEHAJ-ybar)^2)/(sum(SUBESTRATO$EXPR)-1))
  
  Wtilde <- SUBESTRATO$EXPR*(r/sum(SUBESTRATO$EXPR))
  CV2wtilde <- (sum(Wtilde^2)/r)-1
  bbar <- length(SUBESTRATO$FOLIO)/length(unique(SUBESTRATO$SEGMENTO))
  
  #Intra Class Correlation
  MixModel <- lmer(YMONEHAJ ~ 1 + (1 |SEGMENTO), SUBESTRATO)
  MSB <- as.vector(sqrt(unlist(VarCorr(MixModel))))^2
  MSE <- attr(VarCorr(MixModel), "sc")^2
  Rho <- MSB/(MSE+MSB)
  
  YbarModelStrata[i] <- ybar
  VarModelStrata[i] <- (SIGMA2/r)*(1+CV2wtilde)*(1+((bbar-1)*Rho))
}
YbarModelStrata; VarModelStrata 

##Comparisons of Model based formula and designed based formula
##at the ESTRATO levels (when we have 595 strata as we excluded
##the strata with only 1 PSU)

#ESTRATO with 1 PSU
SmallESTRATO
match(c(32012, 104021, 104022, 104031, 104032, 122011, 122012),Uniq_ESTRATO)

###Hint: ESTRATO "50","51","52","53","72","73","205" contains
###only 1 PSU. 
 
VarDesignStrata_nw <- VarDesignStrata[-c(50,51,52,53,72,73,205)]
VarModelStrata 
which(is.na(VarModelStrata))

##From Vaiance based on the model, we can understand that
##there are not any "NaN" values.

par(mfrow=c(1,2))
VARcompESTRATO <- cbind(VarDesignStrata_nw,VarModelStrata)
colnames(VARcompESTRATO) <- c("Design Based","Model Based")

boxplot(VARcompESTRATO[,2],VARcompESTRATO[,1],names=c("Model Based","Design Based"),
        main="Comparison of Variance Among ESTRATO (YMONEHAJ v.)")

##Ratio Comparison of variances for ESTRATO:
summary(VARcompESTRATO[,2]/VARcompESTRATO[,1])
boxplot(VARcompESTRATO[,2]/VARcompESTRATO[,1],
        main="Ratio comparison of Model variance to Design variance
        \n within ESTRATO (YMONEHAJ v.)")
stem(VARcompESTRATO[,2]/VARcompESTRATO[,1])


##Formula of variance based on the model within REGION & ZONA
Region <- unique(HHFRAME$REGION)
Zona <- unique(HHFRAME$ZONA)

VarModelRegionZona <- YbarModelRegionZona <- matrix(0,15,2)

for(i in 1:15){
  for(j in 1:2){
    SUBRegionZona <- HHFRAME[which(HHFRAME$REGION==Region[i] & HHFRAME$ZONA==Zona[j]),]
    
    ybar <- sum(SUBRegionZona$EXPR*SUBRegionZona$YMONEHAJ)/sum(SUBRegionZona$EXPR)
    
    r <- length(SUBRegionZona$FOLIO)
    
    SIGMA2 <- (sum(SUBRegionZona$EXPR*(SUBRegionZona$YMONEHAJ-ybar)^2)/(sum(SUBRegionZona$EXPR)-1))
    
    Wtilde <- SUBRegionZona$EXPR*(r/sum(SUBRegionZona$EXPR))
    CV2wtilde <- (sum(Wtilde^2)/r)-1
    bbar <- length(SUBRegionZona$FOLIO)/length(unique(SUBRegionZona$SEGMENTO))
    
    #Intra Class Correlation
    MixModel <- lmer(YMONEHAJ ~ 1 + (1 |SEGMENTO), SUBRegionZona)
    MSB <- as.vector(sqrt(unlist(VarCorr(MixModel))))^2
    MSE <- attr(VarCorr(MixModel), "sc")^2
    Rho <- MSB/(MSE+MSB)
    
    YbarModelRegionZona[i,j] <- ybar
    VarModelRegionZona[i,j] <- (SIGMA2/r)*(1+CV2wtilde)*(1+((bbar-1)*Rho))
  }
}
#Row stands for REGION and Col stands for ZONA (col1 for "Urbano" & col2 for "Rural") 
YbarModelRegionZona; VarModelRegionZona

VarModelRegionZona <- c(VarModelRegionZona[,1],VarModelRegionZona[,2])
VarDesignRegionZona <- c(VarDesignRegionZona[,1],VarDesignRegionZona[,2])

YbarModelRegionZona <- c(YbarModelRegionZona[,1],YbarModelRegionZona[,2])
RatioRegionZona <- c(RatioRegionZona[,1],RatioRegionZona[,2])

cbind(VarModelRegionZona,VarDesignRegionZona)
par(mfrow=c(1,2))
boxplot(VarModelRegionZona,VarDesignRegionZona,names=c("Model Based","Design Based"),
        main="Comparison of Variance Among REGION+ZONA (YMONEHAJ v.)",col="lightblue")

##Ratio comparison of variances for Region and Zona (Model) and (Design)
VARcompRegZona <- cbind(VarModelRegionZona,VarDesignRegionZona)
summary(VARcompRegZona[,1]/VARcompRegZona[,2])
boxplot(VARcompRegZona[,1]/VARcompRegZona[,2],
        main="Ratio comparison of Model variance to Design variance
        \n within Region and Zona (YMONEHAJ v.)")

##Formula of variance based on the model within COMUNA
COMUNA_NAME <- as.factor(levels(HHFRAME$COMUNA))
NACOMUNA <- tapply(HHFRAME$YMONEHAJ,HHFRAME$COMUNA,FUN=function(x){length(x)})
REMCOMUNA <- names(which(is.na(NACOMUNA)))
FINALCOMUNA <- COMUNA_NAME[! COMUNA_NAME %in% REMCOMUNA]

VarModelCOMUNA <- YbarModelCOMUNA <- rep(0,length(FINALCOMUNA))

for(i in 1:length(FINALCOMUNA)){
  SUBCOMUNA <- subset(HHFRAME,HHFRAME$COMUNA==FINALCOMUNA[i])
  
  ybar <- sum(SUBCOMUNA$EXPC*SUBCOMUNA$YMONEHAJ)/sum(SUBCOMUNA$EXPC)
  
  r <- length(SUBCOMUNA$FOLIO)
  
  SIGMA2 <- (sum(SUBCOMUNA$EXPC*(SUBCOMUNA$YMONEHAJ-ybar)^2)/(sum(SUBCOMUNA$EXPC)-1))
  
  Wtilde <- SUBCOMUNA$EXPC*(r/sum(SUBCOMUNA$EXPC))
  CV2wtilde <- (sum(Wtilde^2)/r)-1
  bbar <- length(SUBCOMUNA$FOLIO)/length(unique(SUBCOMUNA$SEGMENTO))
  
  #Intra Class Correlation
  MixModel <- lmer(YMONEHAJ ~ 1 + (1 |SEGMENTO), SUBCOMUNA)
  MSB <- as.vector(sqrt(unlist(VarCorr(MixModel))))^2
  MSE <- attr(VarCorr(MixModel), "sc")^2
  Rho <- MSB/(MSE+MSB)
  
  YbarModelCOMUNA[i] <- ybar
  VarModelCOMUNA[i] <- (SIGMA2/r)*(1+CV2wtilde)*(1+((bbar-1)*Rho))
}

YbarModelCOMUNA; VarModelCOMUNA
RatioCOMUNA; VarDesignCOMUNA


cbind(YbarModelCOMUNA,RatioCOMUNA)
cbind(VarModelCOMUNA,VarDesignCOMUNA)
par(mfrow=c(1,2))
boxplot(VarModelCOMUNA,VarDesignCOMUNA,names=c("Model Based","Design Based"),
        main="Comparison of Variance Among COMUNA (YMONEHAJ v.)",col="pink")



##Ratio comparison of variances for Comuna (Model) and (Design)
VarcompCOMUNA <- cbind(VarModelCOMUNA,VarDesignCOMUNA)
summary(VarcompCOMUNA[,1]/VarcompCOMUNA[,2])
boxplot(VarcompCOMUNA[,1]/VarcompCOMUNA[,2],
        main="Ratio comparison of Model variance to Design variance
        within Comuna (YMONEHAJ v.)")

par(mfrow=c(1,1))
YbarcompCOMUNA <- cbind(YbarModelCOMUNA,RatioCOMUNA)
boxplot(YbarcompCOMUNA[,1],YbarcompCOMUNA[,2],main="Point Estimation Comparison among COMUNA")


##########################################################
#####STUDY THE RELATIONSHIP of weights and YMONEHAJ#######
##########################################################

par(mfrow=c(1,1))
hist(HHFRAME$EXPR,main="National Level (YMONEHAJ v.)",
     xlim=c(0,1000),xlab="EXPR weight",nclass=200)

summary(HHFRAME$EXPR)

par(mfrow=c(4,4))
sapply(1:15,function(i){SUBRegionHH <- HHFRAME[which(HHFRAME$REGION==Region[i]),];
hist(SUBRegionHH$EXPR,xlim=c(0,1000),xlab="EXPR weight",main=Region[i],nclass=20)})

sapply(1:15,function(i){SUBRegionHH <- HHFRAME[which(HHFRAME$REGION==Region[i]),];
summary(SUBRegionHH$EXPR)})

par(mfrow=c(1,1))
hist(HHFRAME$EXPC,main="National Level",
     xlim=c(0,1000),xlab="EXPC weight",nclass=200)

summary(HHFRAME$EXPC)

par(mfrow=c(4,4))
sapply(1:15,function(i){SUBRegionHH <- HHFRAME[which(HHFRAME$REGION==Region[i]),];
hist(SUBRegionHH$EXPC,xlim=c(0,1000),xlab="EXPC weight",main=Region[i],nclass=20)})

sapply(1:15,function(i){SUBRegionHH <- HHFRAME[which(HHFRAME$REGION==Region[i]),];
summary(SUBRegionHH$EXPC)})


#########################################################
#####STUDY THE CORRELATION of weights and YMONEHAJ#######
#########################################################

cor.test(HHFRAME$EXPR,HHFRAME$YMONEHAJ) #household (National)

cor(HUFRAME$EXPR,HUFRAME$Zvariable) #housing unit (National)

plot(HHFRAME$EXPR,HHFRAME$YMONEHAJ,xlab="EXPR",ylab="YMONEHAJ", 
     main="Correlation Study (National Level)")

COMUNA_NAME <- as.factor(levels(HHFRAME$COMUNA))
NACOMUNA <- tapply(HHFRAME$YMONEHAJ,HHFRAME$COMUNA,FUN=function(x){length(x)})
REMCOMUNA <- names(which(is.na(NACOMUNA)))
FINALCOMUNA <- COMUNA_NAME[! COMUNA_NAME %in% REMCOMUNA]

sapply(1:length(FINALCOMUNA),function(i){
  SUBCOMUNA <- subset(HHFRAME,HHFRAME$COMUNA==FINALCOMUNA[i]);
  cor(SUBCOMUNA$EXPC,SUBCOMUNA$YMONEHAJ)
})

sapply(1:length(FINALCOMUNA),function(i){
  SUBCOMUNA <- subset(HUFRAME,HUFRAME$COMUNA==FINALCOMUNA[i]);
  cor(SUBCOMUNA$EXPC,SUBCOMUNA$Zvariable)
})


Region <- unique(HHFRAME$REGION)
sapply(1:15,function(i){
  SUBRegion <- HHFRAME[which(HHFRAME$REGION==Region[i]),];
  cor(SUBRegion$EXPR,SUBRegion$YMONEHAJ)})

par(mfrow=c(4,4))
sapply(1:15,function(i){
  SUBRegion <- HHFRAME[which(HHFRAME$REGION==Region[i]),];
  plot(SUBRegion$EXPR,SUBRegion$YMONEHAJ,xlab="EXPR",
       ylab="YMONEHAJ", main=Region[i])})

sapply(1:15,function(i){
  SUBRegion <- HUFRAME[which(HUFRAME$REGION==Region[i]),];
  cor(SUBRegion$EXPR,SUBRegion$Zvariable)})


Uniq_ESTRATO <- unique(HUFRAME$ESTRATO)
SmallESTRATO <- as.numeric(names(PSUinSTRAT[PSUinSTRAT==1]))
Uniq_ESTRATOBIG <- Uniq_ESTRATO[!is.element(Uniq_ESTRATO, SmallESTRATO)]

sapply(1:length(Uniq_ESTRATOBIG),function(i){
  SUBESTRATO <- subset(HHFRAME,HHFRAME$ESTRATO==Uniq_ESTRATOBIG[i]);
  cor(SUBESTRATO$EXPR,SUBESTRATO$YMONEHAJ)})

sapply(1:length(Uniq_ESTRATOBIG),function(i){
  SUBESTRATO <- subset(HUFRAME,HUFRAME$ESTRATO==Uniq_ESTRATOBIG[i]);
  cor(SUBESTRATO$EXPR,SUBESTRATO$Zvariable)})

UniqueSegment <- unique(HHFRAME$SEGMENTO)
sapply(1:length(UniqueSegment),function(i){
  cor(HHFRAME$EXPR[HHFRAME$SEGMENTO==UniqueSegment[i]],
      HHFRAME$YMONEHAJ[HHFRAME$SEGMENTO==UniqueSegment[i]])})


##Study the relationship of transformed YMONEHAJ and weight household (National)

hist(HHFRAME$YMONEHAJ,nclass=20,xlab="YMONEHAJ",main="Histogram of YMONEHAJ across
     households")
eps = 1e-5  ## Define parameters for boxcox
boxcox((HHFRAME$YMONEHAJ+eps)~HHFRAME$EXPR,plotit = T,lambda=seq(0.1,0.6,by=0.05))
transf.y <- HHFRAME$YMONEHAJ^0.26
fitYX <- lm(transf.y~HHFRAME$EXPR-1)
summary(fitYX)
plot(fitYX)
anova(fitYX, test="Chisq")  

##likelihood ratio test
null<-lm(transf.y~1,data=HHFRAME)
analysis <- anova(null,fitYX,test="Chisq")
analysis  #Full Model vs. Null Model


########################################################################
###Simulation Study to investigate the Model based formula of variance
###and Design based formula based on the ultimate cluster estimator 
###comes from the survey package.

###Here we assume that the psuedo population is the whole of the Chile.
###The Design is 2stage where the first stage is selecting SEGMENTO and
###the second stage is selecting household.

###For convinience we consider Household frame as our main frame and the
###observational unit (household) is similar to the second stage sampling
###unit.

###We select 20% of households in the selected PSUs.
###Also we would like to compare 4 situations as follows:
### (SRS/SRS "selected PSUs=40")
### (PPS/SYSTEMATIC "selected PSUs=40")
### (SRS/SRS "selected PSUs=10")
### (PPS/SYSTEMATIC "selected PSUs=10")
#########################################################################



###########################################################################
### Simulation Study ###
### SRS/SRS design (selected PSU=40) ###
###########################################################################


#Finding True variance for the National level based on the (SRS/SRS) design:

M <- length(unique(HHFRAME$SEGMENTO))
m <- 40  #sample size of cluster in each region

#number of elements in PSU i:
Ni <- tapply(HHFRAME$FOLIO,HHFRAME$SEGMENTO,FUN=function(x){length(x)})
ni <- (20/100)*Ni  #number of sample elements in PSU i

ti <- tapply(HHFRAME$YMONEHAJ,HHFRAME$SEGMENTO,sum)
tbarU <- sum(ti)/M

S2U1 <- sum((ti-tbarU)^2)/(M-1)
S2U2i <- tapply(HHFRAME$YMONEHAJ,HHFRAME$SEGMENTO,FUN=function(x){var(x)})

M0HH <- length(HHFRAME$FOLIO)
TRUEVARsrssrs40 <- (M*((M-m)/m)*S2U1+(M/m)*sum(Ni*((Ni-ni)/ni)*S2U2i))/(M0HH^2)
TRUEVARsrssrs40

one.rep.Nationalsrssrs40 <- function(){
  m <- 40  #sample size of cluster in each region
  M <- length(unique(HHFRAME$SEGMENTO)) 
  
  SAM <- cluster(data=HHFRAME,
                 clustername="SEGMENTO",size=m,method="srswor",
                 description=TRUE)
  
  #Extract data from the sampled PSUs
  SAMCLUS <- getdata(HHFRAME,SAM)  
  SAMCLUS <- rename(SAMCLUS,c(Prob="pi1"))
  
  #Extract which PSUs are selected
  selectedSEGMENTO <- names(which(table(SAMCLUS$SEGMENTO)!=0))
  
  #Process of selecting SSUs
  ni <- (20/100)*sapply(1:m,function(i){length(HHFRAME$FOLIO[HHFRAME$SEGMENTO==selectedSEGMENTO[i]])})
  Ni <- sapply(1:m,function(i){length(HHFRAME$FOLIO[HHFRAME$SEGMENTO==selectedSEGMENTO[i]])})
  ni <- ceiling(ni)
  
  #treat sample clusters (SEGMENTO) as strata and select SRS sample from each
  S <- strata(data=as.data.frame(SAMCLUS),stratanames="SEGMENTO",
              size=ni,method="srswor")
  
  #Extract the Observed data
  SAMDATA <- getdata(SAMCLUS,S)
  SAMDATA <- rename(SAMDATA,c(Prob="pi2"))
  
  #Extract which SSUs (or Housing units) are selected
  selectedHH <- names(which(table(SAMDATA$FOLIO)!=0))
  
  SAMDATA$wt <- unlist(1/SAMDATA$pi1/SAMDATA$pi2) #final weight
  
  ####Designed Based formula of variance
  options(survey.lonely.psu="adjust")
  ONESTAGE <- svydesign(id=~SEGMENTO, weight=~wt, data=SAMDATA)  ##survey design option
  SVYMEAN <- svymean(SAMDATA$YMONEHAJ,ONESTAGE)
  MEANDESIGNsim <- as.vector(mean(SVYMEAN))
  VARDESIGNsim <- as.vector(SE(SVYMEAN)^2) #(**)Designed Based Variance Based on the simulation
  
  
  ####Model Based preparation
  ybar <- sum(SAMDATA$wt*SAMDATA$YMONEHAJ)/sum(SAMDATA$wt)
  r <- length(SAMDATA$FOLIO)
  SIGMA2 <- (sum(SAMDATA$wt*(SAMDATA$YMONEHAJ-ybar)^2)/(sum(SAMDATA$wt)-1))
  
  Wtilde <- SAMDATA$wt*(r/sum(SAMDATA$wt))
  CV2wtilde <- 0  ##CV2wtilde is zero (just checking)
  bbar <- length(SAMDATA$FOLIO)/length(unique(SAMDATA$SEGMENTO))
  
  #Intra Class Correlation
  MixModel <- lmer(YMONEHAJ ~ 1 + (1 |SEGMENTO), SAMDATA)
  MSB <- as.vector(sqrt(unlist(VarCorr(MixModel))))^2
  MSE <- attr(VarCorr(MixModel), "sc")^2
  Rho <- MSB/(MSE+MSB)
  VARMODELsim <- (SIGMA2/(r-1))*(1+CV2wtilde)*(1+((bbar-1)*Rho)) 
  
  c(MEANDESIGNsim,VARDESIGNsim,VARMODELsim,Rho)   #Final Output of Loop
  
}

## Replicate 1000 sample:
Rep <- 1000
many.reps.Nationalsrssrs40 <- replicate(n=Rep,one.rep.Nationalsrssrs40())
many.reps.Nationalsrssrs40 <- t(many.reps.Nationalsrssrs40)
colnames(many.reps.Nationalsrssrs40) <- c("MEANEDSIGNsim","Design Based","Model Based","Rho")
many.reps.Nationalsrssrs40
dim(many.reps.Nationalsrssrs40)


write.csv(file="C://Users//smosafer//Desktop//YMONEHAJSRSSRS40.csv",(as.matrix(many.reps.Nationalsrssrs40)))
YMONEHAJSRSSRS40 <- read.csv("YMONEHAJSRSSRS40.csv")


###########################################################################
### Simulation Study ###
### PPS/SYSTEMATIC design (selected PSU=40) ###
###########################################################################

one.rep.Nationalppssystematic40 <- function(){
  m <- 40  #sample size of cluster in each region
  
  Ni <- as.vector(tapply(HHFRAME$FOLIO,HHFRAME$SEGMENTO,FUN=function(x){length(x)}))
  
  probi <- m*Ni/sum(Ni)
  
  SAM <- cluster(data=HHFRAME,
                 clustername="SEGMENTO",size=m,method="systematic",
                 pik=probi,description=TRUE)
  
  #Extract data from the sampled PSUs
  SAMCLUS <- getdata(HHFRAME,SAM)  
  SAMCLUS <- rename(SAMCLUS,c(Prob="pi1"))
  
  #Extract which PSUs are selected
  selectedSEGMENTO <- names(which(table(SAMCLUS$SEGMENTO)!=0))
  
  #Process of selecting SSUs
  ni <- (20/100)*sapply(1:m,function(i){length(HHFRAME$FOLIO[HHFRAME$SEGMENTO==selectedSEGMENTO[i]])})
  Nisam <- sapply(1:m,function(i){length(HHFRAME$FOLIO[HHFRAME$SEGMENTO==selectedSEGMENTO[i]])})
  ni <- ceiling(ni)
  
  probii <- sapply(1:length(ni),function(i){rep(sapply(1:length(ni),
                                                       function(i){ni[i]/Nisam[i]})[i],Nisam[i])})
  probii <- unlist(probii)
  
  #treat sample clusters (SEGMENTO) as strata and select SRS sample from each
  S <- strata(data=as.data.frame(SAMCLUS),stratanames="SEGMENTO",
              size=ni,method="systematic",pik=probii)
  
  #Extract the Observed data
  SAMDATA <- getdata(SAMCLUS,S)
  SAMDATA <- rename(SAMDATA,c(Prob="pi2"))
  
  #Extract which SSUs (or Housing units) are selected
  selectedHH <- names(which(table(SAMDATA$FOLIO)!=0))
  
  SAMDATA$wt <- unlist(1/SAMDATA$pi1/SAMDATA$pi2) #final weight
  
  ####Designed Based formula of variance
  options(survey.lonely.psu="adjust")
  ONESTAGE <- svydesign(id=~SEGMENTO, weight=~wt, data=SAMDATA)  ##survey design option
  SVYMEAN <- svymean(SAMDATA$YMONEHAJ,ONESTAGE)
  MEANDESIGNsim <- as.vector(mean(SVYMEAN))
  VARDESIGNsim <- as.vector(SE(SVYMEAN)^2) #(**)Designed Based Variance Based on the simulation
  
  
  ####Model Based preparation
  ybar <- sum(SAMDATA$wt*SAMDATA$YMONEHAJ)/sum(SAMDATA$wt)
  r <- length(SAMDATA$FOLIO)
  SIGMA2 <- (sum(SAMDATA$wt*(SAMDATA$YMONEHAJ-ybar)^2)/(sum(SAMDATA$wt)-1))
  
  Wtilde <- SAMDATA$wt*(r/sum(SAMDATA$wt))
  CV2wtilde <- 0  ##CV2wtilde is zero (just checking)
  bbar <- length(SAMDATA$FOLIO)/length(unique(SAMDATA$SEGMENTO))
  
  #Intra Class Correlation
  MixModel <- lmer(YMONEHAJ ~ 1 + (1 |SEGMENTO), SAMDATA)
  MSB <- as.vector(sqrt(unlist(VarCorr(MixModel))))^2
  MSE <- attr(VarCorr(MixModel), "sc")^2
  Rho <- MSB/(MSE+MSB)
  VARMODELsim <- (SIGMA2/(r-1))*(1+CV2wtilde)*(1+((bbar-1)*Rho)) 
  
  c(MEANDESIGNsim,VARDESIGNsim,VARMODELsim,Rho)   #Final Output of Loop
  
}

## Replicate 1000 sample:
Rep <- 1000
many.reps.Nationalppssystematic40 <- replicate(n=Rep,one.rep.Nationalppssystematic40())
many.reps.Nationalppssystematic40 <- t(many.reps.Nationalppssystematic40)
colnames(many.reps.Nationalppssystematic40) <- c("MEANDESIGN","Design Based","Model Based","Rho")
many.reps.Nationalppssystematic40
dim(many.reps.Nationalppssystematic40)

##Simulated True Variance (should be estimated based on the design info from simulation)
TRUEVARppssystematic40 <- sum((many.reps.Nationalppssystematic40[,1]-mean(many.reps.Nationalppssystematic40[,1]))^2)/1000
TRUEVARppssystematic40


write.csv(file="C://Users//YMONEHAJPPSSYSTEMATIC40.csv",(as.matrix(many.reps.Nationalppssystematic40)))
YMONEHAJPPSSYSTEMATIC40 <- read.csv("YMONEHAJPPSSYSTEMATIC40.csv")


###########################################################################

#Finding True variance for the National level based on the (SRS/SRS) design:

M <- length(unique(HHFRAME$SEGMENTO))
m <- 10  #sample size of cluster in each region

#number of elements in PSU i:
Ni <- tapply(HHFRAME$FOLIO,HHFRAME$SEGMENTO,FUN=function(x){length(x)})
ni <- (20/100)*Ni  #number of sample elements in PSU i

ti <- tapply(HHFRAME$YMONEHAJ,HHFRAME$SEGMENTO,sum)
tbarU <- sum(ti)/M

S2U1 <- sum((ti-tbarU)^2)/(M-1)
S2U2i <- tapply(HHFRAME$YMONEHAJ,HHFRAME$SEGMENTO,FUN=function(x){var(x)})

M0HH <- length(HHFRAME$FOLIO)
TRUEVARSRSSRS10 <- (M*((M-m)/m)*S2U1+(M/m)*sum(Ni*((Ni-ni)/ni)*S2U2i))/(M0HH^2)
TRUEVARSRSSRS10

one.rep.NationalSRSSRS10 <- function(){
  m <- 10  #sample size of cluster in each region
  M <- length(unique(HHFRAME$SEGMENTO)) 
  
  SAM <- cluster(data=HHFRAME,
                 clustername="SEGMENTO",size=m,method="srswor",
                 description=TRUE)
  
  #Extract data from the sampled PSUs
  SAMCLUS <- getdata(HHFRAME,SAM)  
  SAMCLUS <- rename(SAMCLUS,c(Prob="pi1"))
  
  #Extract which PSUs are selected
  selectedSEGMENTO <- names(which(table(SAMCLUS$SEGMENTO)!=0))
  
  #Process of selecting SSUs
  ni <- (20/100)*sapply(1:m,function(i){length(HHFRAME$FOLIO[HHFRAME$SEGMENTO==selectedSEGMENTO[i]])})
  Ni <- sapply(1:m,function(i){length(HHFRAME$FOLIO[HHFRAME$SEGMENTO==selectedSEGMENTO[i]])})
  ni <- ceiling(ni)
  
  #treat sample clusters (SEGMENTO) as strata and select SRS sample from each
  S <- strata(data=as.data.frame(SAMCLUS),stratanames="SEGMENTO",
              size=ni,method="srswor")
  
  #Extract the Observed data
  SAMDATA <- getdata(SAMCLUS,S)
  SAMDATA <- rename(SAMDATA,c(Prob="pi2"))
  
  #Extract which SSUs (or Housing units) are selected
  selectedHH <- names(which(table(SAMDATA$FOLIO)!=0))
  
  SAMDATA$wt <- unlist(1/SAMDATA$pi1/SAMDATA$pi2) #final weight
  
  ####Designed Based formula of variance
  options(survey.lonely.psu="adjust")
  ONESTAGE <- svydesign(id=~SEGMENTO, weight=~wt, data=SAMDATA)  ##survey design option
  SVYMEAN <- svymean(SAMDATA$YMONEHAJ,ONESTAGE)
  MEANDESIGNsim <- as.vector(mean(SVYMEAN))
  VARDESIGNsim <- as.vector(SE(SVYMEAN)^2) #(**)Designed Based Variance Based on the simulation
  
  
  ####Model Based preparation
  ybar <- sum(SAMDATA$wt*SAMDATA$YMONEHAJ)/sum(SAMDATA$wt)
  r <- length(SAMDATA$FOLIO)
  SIGMA2 <- (sum(SAMDATA$wt*(SAMDATA$YMONEHAJ-ybar)^2)/(sum(SAMDATA$wt)-1))
  
  Wtilde <- SAMDATA$wt*(r/sum(SAMDATA$wt))
  CV2wtilde <- 0  ##CV2wtilde is zero (just checking)
  bbar <- length(SAMDATA$FOLIO)/length(unique(SAMDATA$SEGMENTO))
  
  #Intra Class Correlation
  MixModel <- lmer(YMONEHAJ ~ 1 + (1 |SEGMENTO), SAMDATA)
  MSB <- as.vector(sqrt(unlist(VarCorr(MixModel))))^2
  MSE <- attr(VarCorr(MixModel), "sc")^2
  Rho <- MSB/(MSE+MSB)
  VARMODELsim <- (SIGMA2/(r-1))*(1+CV2wtilde)*(1+((bbar-1)*Rho)) 
  
  c(MEANDESIGNsim,VARDESIGNsim,VARMODELsim,Rho)   #Final Output of Loop
  
}

## Replicate 1000 sample:
Rep <- 1000
many.reps.NationalSRSSRS10 <- replicate(n=Rep,one.rep.NationalSRSSRS10())
many.reps.NationalSRSSRS10 <- t(many.reps.NationalSRSSRS10)
colnames(many.reps.NationalSRSSRS10) <- c("MEANDESIGNsim","Design Based","Model Based","Rho")
many.reps.NationalSRSSRS10
dim(many.reps.NationalSRSSRS10)


write.csv(file="C://Users//YMONEHAJSRSSRS10.csv",(as.matrix(many.reps.NationalSRSSRS10)))
YMONEHAJSRSSRS10 <- read.csv("YMONEHAJSRSSRS10.csv")


###########################################################################
### Simulation Study ###
### PPS/SYSTEMATIC design (selected PSU=40) ###
###########################################################################

one.rep.Nationalppssystematic10 <- function(){
  m <- 10  #sample size of cluster in each region
  
  Ni <- as.vector(tapply(HHFRAME$FOLIO,HHFRAME$SEGMENTO,FUN=function(x){length(x)}))
  
  probi <- m*Ni/sum(Ni)
  
  SAM <- cluster(data=HHFRAME,
                 clustername="SEGMENTO",size=m,method="systematic",
                 pik=probi,description=TRUE)
  
  #Extract data from the sampled PSUs
  SAMCLUS <- getdata(HHFRAME,SAM)  
  SAMCLUS <- rename(SAMCLUS,c(Prob="pi1"))
  
  #Extract which PSUs are selected
  selectedSEGMENTO <- names(which(table(SAMCLUS$SEGMENTO)!=0))
  
  #Process of selecting SSUs
  ni <- (20/100)*sapply(1:m,function(i){length(HHFRAME$FOLIO[HHFRAME$SEGMENTO==selectedSEGMENTO[i]])})
  Nisam <- sapply(1:m,function(i){length(HHFRAME$FOLIO[HHFRAME$SEGMENTO==selectedSEGMENTO[i]])})
  ni <- ceiling(ni)
  
  probii <- sapply(1:length(ni),function(i){rep(sapply(1:length(ni),
                                                       function(i){ni[i]/Nisam[i]})[i],Nisam[i])})
  probii <- unlist(probii)
  
  #treat sample clusters (SEGMENTO) as strata and select SRS sample from each
  S <- strata(data=as.data.frame(SAMCLUS),stratanames="SEGMENTO",
              size=ni,method="systematic",pik=probii)
  
  #Extract the Observed data
  SAMDATA <- getdata(SAMCLUS,S)
  SAMDATA <- rename(SAMDATA,c(Prob="pi2"))
  
  #Extract which SSUs (or Housing units) are selected
  selectedHH <- names(which(table(SAMDATA$FOLIO)!=0))
  
  SAMDATA$wt <- unlist(1/SAMDATA$pi1/SAMDATA$pi2) #final weight
  
  ####Designed Based formula of variance
  options(survey.lonely.psu="adjust")
  ONESTAGE <- svydesign(id=~SEGMENTO, weight=~wt, data=SAMDATA)  ##survey design option
  SVYMEAN <- svymean(SAMDATA$YMONEHAJ,ONESTAGE)
  MEANDESIGNsim <- as.vector(mean(SVYMEAN))
  VARDESIGNsim <- as.vector(SE(SVYMEAN)^2) #(**)Designed Based Variance Based on the simulation
  
  
  ####Model Based preparation
  ybar <- sum(SAMDATA$wt*SAMDATA$YMONEHAJ)/sum(SAMDATA$wt)
  r <- length(SAMDATA$FOLIO)
  SIGMA2 <- (sum(SAMDATA$wt*(SAMDATA$YMONEHAJ-ybar)^2)/(sum(SAMDATA$wt)-1))
  
  Wtilde <- SAMDATA$wt*(r/sum(SAMDATA$wt))
  CV2wtilde <- 0  ##CV2wtilde is zero (just checking)
  bbar <- length(SAMDATA$FOLIO)/length(unique(SAMDATA$SEGMENTO))
  
  #Intra Class Correlation
  MixModel <- lmer(YMONEHAJ ~ 1 + (1 |SEGMENTO), SAMDATA)
  MSB <- as.vector(sqrt(unlist(VarCorr(MixModel))))^2
  MSE <- attr(VarCorr(MixModel), "sc")^2
  Rho <- MSB/(MSE+MSB)
  VARMODELsim <- (SIGMA2/(r-1))*(1+CV2wtilde)*(1+((bbar-1)*Rho)) 
  
  c(MEANDESIGNsim,VARDESIGNsim,VARMODELsim,Rho)   #Final Output of Loop
  
}

## Replicate 1000 sample:
Rep <- 1000
many.reps.Nationalppssystematic10 <- replicate(n=Rep,one.rep.Nationalppssystematic10())
many.reps.Nationalppssystematic10 <- t(many.reps.Nationalppssystematic10)
colnames(many.reps.Nationalppssystematic10) <- c("MEANDESIGN","Design Based","Model Based","Rho")
many.reps.Nationalppssystematic10
dim(many.reps.Nationalppssystematic10)

##Simulated True Variance (should be estimated based on the design info from simulation)
TRUEVARppssystematic10 <- sum((many.reps.Nationalppssystematic10[,1]-mean(many.reps.Nationalppssystematic10[,1]))^2)/1000
TRUEVARppssystematic10


write.csv(file="C://Users//YMONEHAJPPSSYSTEMATIC10.csv",(as.matrix(many.reps.Nationalppssystematic10)))
YMONEHAJPPSSYSTEMATIC10 <- read.csv("YMONEHAJPPSSYSTEMATIC10.csv")


##################################################################
########Calling the Stored Simulated Dataset###################### 
########and conducting the analysis###############################



YMONEHAJSRSSRS10 <- read.csv("Q://Data Simulation Casen2009//YMONEHAJSRSSRS10.csv")

YMONEHAJSRSSRS40 <- read.csv("Q://Data Simulation Casen2009//YMONEHAJSRSSRS40.csv")

YMONEHAJPPSSYSTEMATIC10 <- read.csv("Q://Data Simulation Casen2009//YMONEHAJPPSSYSTEMATIC10.csv")

YMONEHAJPPSSYSTEMATIC40 <- read.csv("Q://Data Simulation Casen2009//YMONEHAJPPSSYSTEMATIC40.csv")


TRUEVARSRSSRS10 <- 24323714449
TRUEVARSRSSRS40 <- 6057963185  
TRUEVARppssystematic10 <- (1/1000)*sum((YMONEHAJPPSSYSTEMATIC10[,2]-mean(YMONEHAJPPSSYSTEMATIC10[,2]))^2)
TRUEVARppssystematic10
TRUEVARppssystematic40 <- (1/1000)*sum((YMONEHAJPPSSYSTEMATIC40[,2]-mean(YMONEHAJPPSSYSTEMATIC40[,2]))^2)
TRUEVARppssystematic40

TRUEMEAN <- mean(HHFRAME$YMONEHAJ)


DeviMean = data.frame(devi = c(YMONEHAJSRSSRS10[,2],YMONEHAJSRSSRS40[,2],
                               YMONEHAJPPSSYSTEMATIC10[,2],YMONEHAJPPSSYSTEMATIC40[,2]),
                      method = rep(c("SRS/SRS(10)","SRS/SRS(40)",
                                     "PPS/SYSTEMATIC(10)","PPS/SYSTEMATIC(40)"), each = 1000))
boxplot(devi~method, data = DeviMean,main = "Deviation from the True Mean (YMONEHAJ v.)")
points(1:4,c(mean(YMONEHAJSRSSRS10[,2]),mean(YMONEHAJSRSSRS40[,2]),
             mean(YMONEHAJPPSSYSTEMATIC10[,2]),mean(YMONEHAJPPSSYSTEMATIC40[,2])),col="green",lwd=2,pch=23)
abline(h=TRUEMEAN,col="red",lwd=2)

DeviVar <- data.frame(devi = c(YMONEHAJSRSSRS10[,3],YMONEHAJSRSSRS10[,4],
                               YMONEHAJSRSSRS40[,3],YMONEHAJSRSSRS40[,4],
                               YMONEHAJPPSSYSTEMATIC10[,3],YMONEHAJPPSSYSTEMATIC10[,4],
                               YMONEHAJPPSSYSTEMATIC40[,3],YMONEHAJPPSSYSTEMATIC40[,4]),
                      method = rep(c("SRS(10)D","SRS(10)M",
                                     "SRS(40)D","SRS(40)M",
                                     "PPS(10)D","PPS(10)M",
                                     "PPS(40)D","PPS(40)M"), 
                                   each = 1000))
boxplot(devi~method, data = DeviVar,las=2,main = "Deviation from the True Variances (YMONEHAJ v.)")


par(mfrow=c(1,4))
DeviVarSRS10 <- data.frame(devi = c(YMONEHAJSRSSRS10[,3],YMONEHAJSRSSRS10[,4]),
                           method = rep(c("SRS(10)D","SRS(10)M"), 
                                        each = 1000))
boxplot(devi~method, data = DeviVarSRS10,main = "Variance SRS/SRS(10)"
        )
points(1:2,c(mean(YMONEHAJSRSSRS10[,3]),mean(YMONEHAJSRSSRS10[,4])),col="green",lwd=2,pch=23)
abline(h=TRUEVARSRSSRS10,col="red",lwd=2)

DeviVarSRS40 <- data.frame(devi = c(YMONEHAJSRSSRS40[,3],YMONEHAJSRSSRS40[,4]),
                           method = rep(c("SRS(40)D","SRS(40)M"), 
                                        each = 1000))
boxplot(devi~method, data = DeviVarSRS40,main = "Variance SRS/SRS(40)"
        )
points(1:2,c(mean(YMONEHAJSRSSRS40[,3]),mean(YMONEHAJSRSSRS40[,4])),col="green",lwd=2,pch=23)
abline(h=TRUEVARSRSSRS40,col="red",lwd=2)

DeviVarPPS10 <- data.frame(devi = c(YMONEHAJPPSSYSTEMATIC10[,3],YMONEHAJPPSSYSTEMATIC10[,4]),
                           method = rep(c("PPS(10)D","PPS(10)M"), 
                                        each = 1000))
boxplot(devi~method, data = DeviVarPPS10,main = "Variance PPS/SYSTEMATIC(10)"
        )
points(1:2,c(mean(YMONEHAJPPSSYSTEMATIC10[,3]),mean(YMONEHAJPPSSYSTEMATIC10[,4])),col="green",lwd=2,pch=23)
abline(h=TRUEVARppssystematic10,col="red",lwd=2)

DeviVarPPS40 <- data.frame(devi = c(YMONEHAJPPSSYSTEMATIC40[,3],YMONEHAJPPSSYSTEMATIC40[,4]),
                           method = rep(c("PPS(40)D","PPS(40)M"), 
                                        each = 1000))
boxplot(devi~method, data = DeviVarPPS40,main = "Variance PPS/SYSTEMATIC(40)"
        )
points(1:2,c(mean(YMONEHAJPPSSYSTEMATIC40[,3]),mean(YMONEHAJPPSSYSTEMATIC40[,4])),col="green",lwd=2,pch=23)
abline(h=TRUEVARppssystematic40,col="red",lwd=2)


table(YMONEHAJSRSSRS10[,5]==0)  #Number of situations getting 0 for Rho
length(YMONEHAJSRSSRS10[,5][is.na(YMONEHAJSRSSRS10[,5])])  ##NA values for Rho

table(YMONEHAJSRSSRS40[,5]==0)  #Number of situations getting 0 for Rho
length(YMONEHAJSRSSRS40[,5][is.na(YMONEHAJSRSSRS40[,5])])  ##NA values for Rho

table(YMONEHAJPPSSYSTEMATIC10[,5]==0)  #Number of situations getting 0 for Rho
length(YMONEHAJPPSSYSTEMATIC10[,5][is.na(YMONEHAJPPSSYSTEMATIC10[,5])])  ##NA values for Rho

table(YMONEHAJPPSSYSTEMATIC40[,5]==0)  #Number of situations getting 0 for Rho
length(YMONEHAJPPSSYSTEMATIC40[,5][is.na(YMONEHAJPPSSYSTEMATIC40[,5])])  ##NA values for Rho


############################################################################
#### Informative and Non-informative Samplind Study
############################################################################

#### 1. Real Data Analysis (National Level)

## Model1) At the National Level 
## yij=alpha_i+eij; alpha_i is SEGMENTO random effect
MixModel1 <- lmer(YMONEHAJ ~ 1 + (1|SEGMENTO), HHFRAME)
COR11 <- cor.test(residuals(MixModel1),HHFRAME$EXPR)
COR12 <- cor.test((residuals(MixModel1)^2),HHFRAME$EXPR)
cor.test(HHFRAME$YMONEHAJ,HHFRAME$EXPR)
COR11$statistic; COR11$estimate; COR11$p.value
COR12$statistic; COR12$estimate; COR12$p.value

## Model2) At the National Level 
## yijk=alpha_i+beta_ij+eijk; beta_ij is nested in alpha_i (both are random)
## where alpha_i is ESTRATO effect & beta_ij is SEGMENTO effect
MixModel2 <- lmer(YMONEHAJ~1+(1|ESTRATO/SEGMENTO),data=HHFRAME)
COR21 <- cor.test(residuals(MixModel2),HHFRAME$EXPR)
COR22 <- cor.test((residuals(MixModel2)^2),HHFRAME$EXPR)
cor.test(HHFRAME$YMONEHAJ,HHFRAME$EXPR)
COR21$statistic; COR21$estimate; COR21$p.value
COR22$statistic; COR22$estimate; COR22$p.value

## Model3) At the National Level
## yijk=alpha_i+beta_j+eijk:
## alpha_i is ESTRATO effect (fixed) & beta_j is SEGMENTO effect (random)
MixModel3 <- lmer(YMONEHAJ~ESTRATO+(1|SEGMENTO),data=HHFRAME)
COR31 <- cor.test(residuals(MixModel3),HHFRAME$EXPR)
COR32 <- cor.test((residuals(MixModel3)^2),HHFRAME$EXPR)
cor.test(HHFRAME$YMONEHAJ,HHFRAME$EXPR)
COR31$statistic; COR31$estimate; COR31$p.value
COR32$statistic; COR32$estimate; COR32$p.value


## Model4) At the National Level
## yijk=alpha_i+beta_j+eijk:
## alpha_i is ESTRATO effect (random) & beta_j is SEGMENTO effect (random)
MixModel4 <- lmer(YMONEHAJ~(1|ESTRATO)+(1|SEGMENTO),data=HHFRAME)
COR41 <- cor.test(residuals(MixModel4),HHFRAME$EXPR)
COR42 <- cor.test((residuals(MixModel4)^2),HHFRAME$EXPR)
cor.test(HHFRAME$YMONEHAJ,HHFRAME$EXPR)
COR41$statistic; COR41$estimate; COR41$p.value
COR42$statistic; COR42$estimate; COR42$p.value


## Model5) At the National Level
## yijk=alpha_i+beta_ij+eijk:
## alpha_i is ESTRATO effect (fixed) & beta_ij is SEGMENTO effect (random & nested in ESTRATO)
## BUT we have a problem here since ESTRATO can be considered in both Random effects and
## Fixed effects for the output.
MixModel5 <- lmer(YMONEHAJ~ESTRATO+(1|ESTRATO/SEGMENTO),data=HHFRAME)
COR51 <- cor.test(residuals(MixModel5),HHFRAME$EXPR)
COR52 <- cor.test((residuals(MixModel5)^2),HHFRAME$EXPR)
cor.test(HHFRAME$YMONEHAJ,HHFRAME$EXPR)
COR51$statistic; COR51$estimate; COR51$p.value
COR52$statistic; COR52$estimate; COR52$p.value


#### 1. Real Data Analysis (COMUNA Level)

## Model1) At the COMUNA Level 
## yij=alpha_i+eij; alpha_i is SEGMENTO random effect
COMUNA_NAME <- as.factor(levels(HHFRAME$COMUNA))
NACOMUNA <- tapply(HHFRAME$YMONEHAJ,HHFRAME$COMUNA,FUN=function(x){length(x)})
REMCOMUNA <- names(which(is.na(NACOMUNA)))
FINALCOMUNA <- COMUNA_NAME[! COMUNA_NAME %in% REMCOMUNA]

COR11statistic <- COR11estimate <- COR11p.value <- 
  COR12statistic <- COR12estimate <- COR12p.value <- rep(0,length(FINALCOMUNA))

for(i in 1:length(FINALCOMUNA)){
  SUBCOMUNA <- subset(HHFRAME,HHFRAME$COMUNA==FINALCOMUNA[i])
  MixModel1 <- lmer(YMONEHAJ ~ 1 + (1|SEGMENTO), SUBCOMUNA)
  COR11 <- cor.test(residuals(MixModel1),SUBCOMUNA$EXPR)
  COR12 <- cor.test((residuals(MixModel1)^2),SUBCOMUNA$EXPR)
  COR11statistic[i] <- as.vector(COR11$statistic) 
  COR11estimate[i] <- as.vector(COR11$estimate)
  COR11p.value[i] <- as.vector(COR11$p.value)
  COR12statistic[i] <- as.vector(COR12$statistic)
  COR12estimate[i] <- as.vector(COR12$estimate)
  COR12p.value[i] <- as.vector(COR12$p.value)
}
COR11statistic ; COR11estimate ; COR11p.value  
COR12statistic ; COR12estimate ; COR12p.value


## Model2) At the COMUNA Level 
## yijk=alpha_i+beta_ij+eijk; beta_ij is nested in alpha_i (both are random)
## where alpha_i is ESTRATO effect & beta_ij is SEGMENTO effect
COR21statistic <- COR21estimate <- COR21p.value <- 
  COR22statistic <- COR22estimate <- COR22p.value <- rep(0,length(FINALCOMUNA))

for(i in 1:length(FINALCOMUNA)){
  SUBCOMUNA <- subset(HHFRAME,HHFRAME$COMUNA==FINALCOMUNA[i])
  MixModel2 <- lmer(YMONEHAJ~1+(1|ESTRATO/SEGMENTO), SUBCOMUNA)
  COR21 <- cor.test(residuals(MixModel2),SUBCOMUNA$EXPR)
  COR22 <- cor.test((residuals(MixModel2)^2),SUBCOMUNA$EXPR)
  COR21statistic[i] <- as.vector(COR21$statistic) 
  COR21estimate[i] <- as.vector(COR21$estimate)
  COR21p.value[i] <- as.vector(COR21$p.value)
  COR22statistic[i] <- as.vector(COR22$statistic)
  COR22estimate[i] <- as.vector(COR22$estimate)
  COR22p.value[i] <- as.vector(COR22$p.value)
}
COR21statistic ; COR21estimate ; COR21p.value  
COR22statistic ; COR22estimate ; COR22p.value

## Model3) At the COMUNA Level
## yijk=alpha_i+beta_j+eijk:
## alpha_i is ESTRATO effect (fixed) & beta_j is SEGMENTO effect (random)
COR31statistic <- COR31estimate <- COR31p.value <- 
  COR32statistic <- COR32estimate <- COR32p.value <- rep(0,length(FINALCOMUNA))

for(i in 1:length(FINALCOMUNA)){
  SUBCOMUNA <- subset(HHFRAME,HHFRAME$COMUNA==FINALCOMUNA[i])
  MixModel3 <- lmer(YMONEHAJ~ESTRATO+(1|SEGMENTO), data=SUBCOMUNA)
  COR31 <- cor.test(residuals(MixModel3),SUBCOMUNA$EXPR)
  COR32 <- cor.test((residuals(MixModel3)^2),SUBCOMUNA$EXPR)
  COR31statistic[i] <- as.vector(COR31$statistic) 
  COR31estimate[i] <- as.vector(COR31$estimate)
  COR31p.value[i] <- as.vector(COR31$p.value)
  COR32statistic[i] <- as.vector(COR32$statistic)
  COR32estimate[i] <- as.vector(COR32$estimate)
  COR32p.value[i] <- as.vector(COR32$p.value)
}
COR31statistic ; COR31estimate ; COR31p.value  
COR32statistic ; COR32estimate ; COR32p.value



## Model4) At the COMUNA Level
## yijk=alpha_i+beta_j+eijk:
## alpha_i is ESTRATO effect (random) & beta_j is SEGMENTO effect (random)
COR41statistic <- COR41estimate <- COR41p.value <- 
  COR42statistic <- COR42estimate <- COR42p.value <- rep(0,length(FINALCOMUNA))

for(i in 1:length(FINALCOMUNA)){
  SUBCOMUNA <- subset(HHFRAME,HHFRAME$COMUNA==FINALCOMUNA[i])
  MixModel4 <- lmer(YMONEHAJ~(1|ESTRATO)+(1|SEGMENTO), data=SUBCOMUNA)
  COR41 <- cor.test(residuals(MixModel4),SUBCOMUNA$EXPR)
  COR42 <- cor.test((residuals(MixModel4)^2),SUBCOMUNA$EXPR)
  COR41statistic[i] <- as.vector(COR41$statistic) 
  COR41estimate[i] <- as.vector(COR41$estimate)
  COR41p.value[i] <- as.vector(COR41$p.value)
  COR42statistic[i] <- as.vector(COR42$statistic)
  COR42estimate[i] <- as.vector(COR42$estimate)
  COR42p.value[i] <- as.vector(COR42$p.value)
}
COR41statistic ; COR41estimate ; COR41p.value  
COR42statistic ; COR42estimate ; COR42p.value


## Model5) At the COMUNA Level
## yijk=alpha_i+beta_ij+eijk:
## alpha_i is ESTRATO effect (fixed) & beta_ij is SEGMENTO effect (random & nested in ESTRATO)
## BUT we have a problem here since ESTRATO can be considered in both Random effects and
## Fixed effects for the output.
COR51statistic <- COR51estimate <- COR51p.value <- 
  COR52statistic <- COR52estimate <- COR52p.value <- rep(0,length(FINALCOMUNA))

for(i in 1:length(FINALCOMUNA)){
  SUBCOMUNA <- subset(HHFRAME,HHFRAME$COMUNA==FINALCOMUNA[i])
  MixModel5 <- lmer(YMONEHAJ~ESTRATO+(1|ESTRATO/SEGMENTO), data=SUBCOMUNA)
  COR51 <- cor.test(residuals(MixModel5),SUBCOMUNA$EXPR)
  COR52 <- cor.test((residuals(MixModel5)^2),SUBCOMUNA$EXPR)
  COR51statistic[i] <- as.vector(COR51$statistic) 
  COR51estimate[i] <- as.vector(COR51$estimate)
  COR51p.value[i] <- as.vector(COR51$p.value)
  COR52statistic[i] <- as.vector(COR52$statistic)
  COR52estimate[i] <- as.vector(COR52$estimate)
  COR52p.value[i] <- as.vector(COR52$p.value)
}
COR51statistic ; COR51estimate ; COR51p.value  
COR52statistic ; COR52estimate ; COR52p.value


## Model6) At the COMUNA Level; SEGMENTO is random effect (Logistic Model)
COR61statistic <- COR61estimate <- COR61p.value <- 
  COR62statistic <- COR62estimate <- COR62p.value <- rep(0,length(FINALCOMUNA))

for(i in 1:length(FINALCOMUNA)){
  SUBCOMUNA <- subset(HHFRAME,HHFRAME$COMUNA==FINALCOMUNA[i])
  MixModel6 <- glmer(YMONEHAJ ~ 1 + (1 |SEGMENTO), SUBCOMUNA,family = binomial)
  COR61 <- cor.test(residuals(MixModel6),SUBCOMUNA$EXPR)
  COR62 <- cor.test((residuals(MixModel6)^2),SUBCOMUNA$EXPR)
  COR61statistic[i] <- as.vector(COR61$statistic) 
  COR61estimate[i] <- as.vector(COR61$estimate)
  COR61p.value[i] <- as.vector(COR61$p.value)
  COR62statistic[i] <- as.vector(COR62$statistic)
  COR62estimate[i] <- as.vector(COR62$estimate)
  COR62p.value[i] <- as.vector(COR62$p.value)
}
COR61statistic ; COR61estimate ; COR61p.value  
COR62statistic ; COR62estimate ; COR62p.value






