# ============================================================================================
# Name: Program for Estimating the 2013 Poverty Rate Variance Estimator from Chilean Dataset 
# Author: Sepideh Mosaferi
# Date: Spring 2016
# ============================================================================================

require(stats) 
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

casendata2009 <- read.spss("Q://Casen Whole Datasets//Casen Whole Datasets//Casen2009.sav")
casendata2009 <- as.data.frame(casendata2009)

casendata2011 <- read.spss("Q://Casen Whole Datasets//Casen Whole Datasets//Casen2011.sav")
casendata2011 <- as.data.frame(casendata2011)

casendata2013 <- read.spss("Q://Casen Whole Datasets//Casen Whole Datasets//Casen2013.sav")
casendata2013 <- as.data.frame(casendata2013)


# Sort datasets
SORTDATA2013 <- casendata2013[order(casendata2013$region,casendata2013$comuna),]
SORTDATA2013 <- as.data.frame(SORTDATA2013)

SORTDATA2011 <- casendata2011[order(casendata2011$region,casendata2011$comuna),]
SORTDATA2011 <- as.data.frame(SORTDATA2011)

SORTDATA2009 <- casendata2009[order(casendata2009$REGION,casendata2009$COMUNA),]
SORTDATA2009 <- as.data.frame(SORTDATA2009)


# Assessing the Missing Part of dataset based on the poor identification
  ## 2013
MISS2013 <- subset(SORTDATA2013,is.na(SORTDATA2013$pobreza_MN),select=c(edad,sexo,pobreza_MN,
                                                         ypchautcor,ypchtot,folio)) 
length(unique(MISS2013$folio))
ASSMISS2013 <- sapply(1:nrow(MISS2013),function(i)
  {SORTDATA2013$pobreza_MN[SORTDATA2013$folio==MISS2013$folio[i]]})  # Assessing Missing Situation
SUBSTITUTE2013 <- sapply(1:length(ASSMISS2013),
                         function(i){ASSMISS2013[[i]][is.na(ASSMISS2013[[i]])] <- ASSMISS2013[[i]][1]})
SORTDATA2013$pobreza_MN[is.na(SORTDATA2013$pobreza_MN)] <- SUBSTITUTE2013

  ## 2011
MISS2011 <- subset(SORTDATA2011,is.na(SORTDATA2011$corte),select=c(edad,sexo,corte,
                                                                     ymonehaj,yaimhaj,folio)) 
length(unique(MISS2011$folio))
ASSMISS2011 <- sapply(1:nrow(MISS2011),function(i)
{SORTDATA2011$corte[SORTDATA2011$folio==MISS2011$folio[i]]})  # Assessing Missing Situation
SUBSTITUTE2011 <- sapply(1:length(ASSMISS2011),
                         function(i){ASSMISS2011[[i]][is.na(ASSMISS2011[[i]])] <- ASSMISS2011[[i]][1]})
SORTDATA2011$corte[is.na(SORTDATA2011$corte)] <- SUBSTITUTE2011

    ## 2009
MISS2009 <- subset(SORTDATA2009,is.na(SORTDATA2009$CORTE),select=c(EDAD,SEXO,CORTE,
                                                       YMONEHAJ,YAIMHAJ,FOLIO))  
length(unique(MISS2009$FOLIO))
ASSMISS2009 <- sapply(1:nrow(MISS2009),function(i)
  {SORTDATA2009$CORTE[SORTDATA2009$FOLIO==MISS2009$FOLIO[i]]})  # Assessing Missing Situation
SUBSTITUTE2009 <- sapply(1:length(ASSMISS2009),
                         function(i){ASSMISS2009[[i]][is.na(ASSMISS2009[[i]])] <- ASSMISS2009[[i]][1]})
SORTDATA2009$CORTE[is.na(SORTDATA2009$CORTE)] <- SUBSTITUTE2009


## Producing the Frame based on the Households
HHFRAME2013 <- SORTDATA2013[!duplicated(SORTDATA2013$folio),]
HHFRAME2011 <- SORTDATA2011[!duplicated(SORTDATA2011$folio),]
HHFRAME2009 <- SORTDATA2009[!duplicated(SORTDATA2009$FOLIO),]


#### Producing psuedo population based on the largest region ####
table(HHFRAME2013$region)  # finding the largest region
max(table(HHFRAME2013$region))  # "Metropolitana": the largest one
# Take Metropolitana as the psuedo population

psuedopop2013 <- data.frame(HHFRAME2013[HHFRAME2013$region=="Metropolitana",])
psuedopop2013$Ipoor <- ifelse(psuedopop2013$pobreza_MN=="No pobres",0,1)

psuedopop2011 <- data.frame(HHFRAME2011[HHFRAME2011$region=="Regi?n Metropolitana",])
psuedopop2011$Ipoor <- ifelse(psuedopop2011$corte=="No pobres",0,1)

psuedopop2009 <- data.frame(HHFRAME2009[HHFRAME2009$REGION=="Regi?n Metropolitana",])
psuedopop2009$Ipoor <- ifelse(psuedopop2009$CORTE=="No pobre",0,1)

#### Sample size allocation
    ## 2013
Nh2013 <- as.vector(table(psuedopop2013$comuna))  # COMUNAs as STRATA
Nh2013 <- Nh2013[Nh2013!=0]
N2013 <- sum(Nh2013)
n2013 <- floor((10/100)*N2013)  # Overall sample size
ALLOC2013 <- strAlloc(n.tot=n2013,Nh=Nh2013,alloc="prop")  # Proportional Allocation
nh2013 <- floor(ALLOC2013$nh+0.5) # final selectd sample size
Wh2013 <- Nh2013/N2013
fh2013 <- nh2013/Nh2013
  ## 2011
Nh2011 <- as.vector(table(psuedopop2011$comuna))  # COMUNAs as STRATA
Nh2011 <- Nh2011[Nh2011!=0]
N2011 <- sum(Nh2011)
n2011 <- floor((10/100)*N2011)  # Overall sample size
ALLOC2011 <- strAlloc(n.tot=n2011,Nh=Nh2011,alloc="prop")  #Proportional Allocation
nh2011 <- floor(ALLOC2011$nh+0.5) #final selectd sample size
Wh2011 <- Nh2011/N2011
fh2011 <- nh2011/Nh2011
  ## 2009
Nh2009 <- as.vector(table(psuedopop2009$COMUNA))  # COMUNAs as STRATA
Nh2009 <- Nh2009[Nh2009!=0]
N2009 <- sum(Nh2009)
n2009 <- floor((10/100)*N2009)  # Overall sample size
ALLOC2009 <- strAlloc(n.tot=n2009,Nh=Nh2009,alloc="prop")  # Proportional Allocation
nh2009 <- floor(ALLOC2009$nh+0.5) # final selectd sample size
Wh2009 <- Nh2009/N2009
fh2009 <- nh2009/Nh2009

unique(psuedopop2013$comuna[psuedopop2013$region=="Metropolitana"])
unique(psuedopop2011$comuna[psuedopop2011$region=="Regi?n Metropolitana"])
unique(psuedopop2009$COMUNA[psuedopop2009$REGION=="Regi?n Metropolitana"])

 # first and last comuna in Metropolitana
which(levels(psuedopop2013$comuna)=="Santiago")
which(levels(psuedopop2013$comuna)=="Pe?aflor")
which(levels(psuedopop2011$comuna)=="Santiago")
which(levels(psuedopop2011$comuna)=="Pe?aflor")
which(levels(psuedopop2009$COMUNA)=="Santiago")
which(levels(psuedopop2009$COMUNA)=="Pe?aflor")

  # comunas nested in Urbano, Rural, and both Zona

which(levels(psuedopop2013$comuna)=="Santiago") # beginning of urbano
which(levels(psuedopop2013$comuna)=="Puente Alto") # end of urbano
which(levels(psuedopop2013$comuna)=="Pirque") # begining of both
which(levels(psuedopop2013$comuna)=="Mar?a Pinto") 
which(levels(psuedopop2013$comuna)=="Talagante")
which(levels(psuedopop2013$comuna)=="Pe?aflor") # end of both
which(levels(psuedopop2013$comuna)=="San Pedro") # just rural
COMUNAinMetro <- levels(psuedopop2013$comuna)[279:330]
table(psuedopop2013$comuna,psuedopop2013$zona)

P.true2013 <- aggregate(psuedopop2013$Ipoor[psuedopop2013$region=="Metropolitana"]
          ~psuedopop2013$comuna[psuedopop2013$region=="Metropolitana"],
          FUN=mean)[,2]
Var.True2013 <- P.true2013*(1-P.true2013)/nh2013  # Theoretical Variance 2013

one.rep <- function(){
  ## Sample selection
      ## 2013
  SMP.IDs2013 <- strata(data=psuedopop2013,
                    stratanames="comuna",size=nh2013,method="srswor")
  SAMPLE2013 <- getdata(psuedopop2013,SMP.IDs2013)  # Output of selected sample
  P.Direct2013 <- tapply(SAMPLE2013$Ipoor,INDEX=SAMPLE2013$comuna,FUN=mean) # direct estimate
  P.Direct2013 <- P.Direct2013[279:330]  # for Metropolitan Region
  Var.Direct2013 <- P.Direct2013*(1-P.Direct2013)/nh2013  # direct variance per stratum
      ## 2011
  SMP.IDs2011 <- strata(data=psuedopop2011,
                        stratanames="comuna",size=nh2011,method="srswor")
  SAMPLE2011 <- getdata(psuedopop2011,SMP.IDs2011)  # Output of selected sample
  P.Direct2011 <- tapply(SAMPLE2011$Ipoor,INDEX=SAMPLE2011$comuna,FUN=mean) # direct estimate
  P.Direct2011 <- P.Direct2011[279:330]  # for Metropolitan Region
  Var.Direct2011 <- P.Direct2011*(1-P.Direct2011)/nh2011  # direct variance per stratum
      ## 2009
  SMP.IDs2009 <- strata(data=psuedopop2009,
                        stratanames="COMUNA",size=nh2009,method="srswor")
  SAMPLE2009 <- getdata(psuedopop2009,SMP.IDs2009)  # Output of selected sample
  P.Direct2009 <- tapply(SAMPLE2009$Ipoor,INDEX=SAMPLE2009$COMUNA,FUN=mean) # direct estimate
  P.Direct2009 <- P.Direct2009[278:329]  # for Metropolitan Region
  Var.Direct2009 <- P.Direct2009*(1-P.Direct2009)/nh2009  # direct variance per stratum

  P.Synth <- sum(SAMPLE2013$Ipoor)/nrow(SAMPLE2013)  # synthetic estimate
  Var.Synth <- P.Synth*(1-P.Synth)/nh2013  # synthetic variance per stratum
  
  ## Synthetic By Zona
  SUMZONA <- aggregate(SAMPLE2013$Ipoor~SAMPLE2013$comuna+SAMPLE2013$zona,FUN=sum)  
  nZONA <- aggregate(SAMPLE2013$Ipoor~SAMPLE2013$comuna+SAMPLE2013$zona,FUN=length)
  FRAMEZona <- data.frame(SUMZONA[,1],SUMZONA[,2],SUMZONA[,3],nZONA[,3])
  colnames(FRAMEZona) <- c("COMUNA","ZONA","sumIpoor","sumn")
  
  P.SynthUrbano <- sum(FRAMEZona$sumIpoor[FRAMEZona$COMUNA %in% COMUNAinMetro[1:33]])/
  sum(FRAMEZona$sumn[FRAMEZona$COMUNA %in% COMUNAinMetro[1:33]])
  
  P.SynthBoth <- sum(FRAMEZona$sumIpoor[FRAMEZona$COMUNA %in% c(COMUNAinMetro[34:46],
                                                                COMUNAinMetro[48:52])])/
    sum(FRAMEZona$sumn[FRAMEZona$COMUNA %in% c(COMUNAinMetro[34:46],
                                               COMUNAinMetro[48:52])])
  
  P.SynthRural <- sum(FRAMEZona$sumIpoor[FRAMEZona$COMUNA %in% COMUNAinMetro[47]])/
    sum(FRAMEZona$sumn[FRAMEZona$COMUNA %in% COMUNAinMetro[47]])
    
  P.SynthZona <- c(rep(P.SynthUrbano,33),rep(P.SynthBoth,13),P.SynthRural,rep(P.SynthBoth,5))
  Var.SynthZona <- P.SynthZona*(1-P.SynthZona)/nh2013  
    
  # time-series estimate
  P.Time <- (nh2013*P.Direct2013+nh2011*P.Direct2011+nh2009*P.Direct2009)/(nh2013+nh2011+nh2009) 
  Var.Time <- P.Time*(1-P.Time)/nh2013  # time-series variance per stratum
  
  cbind(Var.Direct2013,Var.Synth,Var.SynthZona,Var.Time)
}
R <- 1000
many.reps <- replicate(n=R,one.rep()); dim(many.reps)  ## 52*4*1000  each row corresponds to a comuna; 
                                                       ## last row is aggregated stratifies variance


  ## Boxplot per comuna
sapply(1:52,function(i){boxplot(many.reps[i,1,],many.reps[i,2,],many.reps[i,3,],
                                many.reps[i,4,],names=
                                  c("Direct","Synthetic","Synth-Zona","Time-Series"),
main="Boxplot of Variances per Comuna in Metropolitana Region")
points(1:4,c(mean(many.reps[i,1,]),mean(many.reps[i,2,]),
             mean(many.reps[i,3,]),mean(many.reps[i,4,])),col="green",lwd=2,pch=23)
abline(h=Var.True2013[i],col="red",lwd=2)})

  # Strip Chart for the empirical mean
Empmean <- data.frame(sapply(1:52,function(i){mean(many.reps[i,1,])}),
                      sapply(1:52,function(i){mean(many.reps[i,2,])}),
                      sapply(1:52,function(i){mean(many.reps[i,3,])}),
                      sapply(1:52,function(i){mean(many.reps[i,4,])})) 
par(mfrow=c(1,1)) 
colnames(Empmean) <- c("Direct","Synthetic","Synth-Zona","Time-Series")
stripchart(data.frame(scale(Empmean)),method="jitter",
           vertical=TRUE,main="The Empirical Mean of Variances per comuna")

  # Relative Bias
par(mfrow=c(1,2))
RelBias.Direct <- sapply(1:52,function(i){(mean(many.reps[i,1,]-Var.True2013[i])/Var.True2013[i])*100}) 
RelBias.Synth <- sapply(1:52,function(i){(mean(many.reps[i,2,]-Var.True2013[i])/Var.True2013[i])*100}) 
RelBias.SynthZona <- sapply(1:52,function(i){(mean(many.reps[i,3,]-Var.True2013[i])/Var.True2013[i])*100}) 
RelBias.Time <- sapply(1:52,function(i){(mean(many.reps[i,4,]-Var.True2013[i])/Var.True2013[i])*100}) 

RelBias <- cbind(RelBias.Direct,RelBias.Synth,RelBias.SynthZona,RelBias.Time)
RelBias
colnames(RelBias) <- c("Direct","Synthetic","Synth-Zona","Time")
boxplot(RelBias,las=1,main="Relative Bias per Comuna")

  ## True Mean Squared Error Comparison: E((estimator-True)^2)
TrueMSE.Direct <- sapply(1:52,function(i){mean((many.reps[i,1,]-Var.True2013[i])^2)}) 
TrueMSE.Synth <- sapply(1:52,function(i){mean((many.reps[i,2,]-Var.True2013[i])^2)}) 
TrueMSE.SynthZona <- sapply(1:52,function(i){mean((many.reps[i,3,]-Var.True2013[i])^2)}) 
TrueMSE.Time <- sapply(1:52,function(i){mean((many.reps[i,4,]-Var.True2013[i])^2)}) 
TrueMSE <- cbind(TrueMSE.Direct,TrueMSE.Synth,TrueMSE.SynthZona,TrueMSE.Time)
TrueMSE
SQRTMSE <- sqrt(TrueMSE); SQRTMSE
colnames(SQRTMSE) <- c("Direct","Synthetic","Synth-Zona","Time")
boxplot(SQRTMSE,las=1,main="Root of True MSE per Comuna")







