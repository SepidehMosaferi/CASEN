# ============================================================================================
# Name: R code for the CASEN 1990-2013 for estimating the survey weighted poverty rates 
#       at the COMUNAs   
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
require(fda)

casendata1990 <- read.spss("Q://Casen Whole Datasets//Casen1990.sav")
casendata1990 <- as.data.frame(casendata1990)

casendata1992 <- read.spss("Q://Casen Whole Datasets//Casen1992.sav")
casendata1992 <- as.data.frame(casendata1992)

casendata1994 <- read.spss("Q://Casen Whole Datasets//Casen1994.sav")
casendata1994 <- as.data.frame(casendata1994)

casendata1996 <- read.spss("Q://Casen Whole Datasets//Casen1996.sav")
casendata1996 <- as.data.frame(casendata1996)

casendata1998 <- read.spss("Q://Casen Whole Datasets//Casen1998.sav")
casendata1998 <- as.data.frame(casendata1998)

casendata2000 <- read.spss("Q://Casen Whole Datasets//Casen2000.sav")
casendata2000 <- as.data.frame(casendata2000)

casendata2003 <- read.spss("Q://Casen Whole Datasets//Casen2003.sav")
casendata2003 <- as.data.frame(casendata2003)

casendata2006 <- read.spss("Q://Casen Whole Datasets//Casen2006.sav")
casendata2006 <- as.data.frame(casendata2006)

casendata2009 <- read.spss("Q://Casen Whole Datasets//Casen2009.sav")
casendata2009 <- as.data.frame(casendata2009)

casendata2011 <- read.spss("Q://Casen Whole Datasets//Casen2011.sav")
casendata2011 <- as.data.frame(casendata2011)

casendata2013 <- read.spss("Q://Casen Whole Datasets//Casen2013.sav")
casendata2013 <- as.data.frame(casendata2013)


casendata1990$Ipoor1990 <- ifelse(casendata1990$corte=="No pobre",0,1)  
sumnum1990 <- tapply(casendata1990$expc*casendata1990$Ipoor1990,casendata1990$comu,FUN=sum,na.rm=T)
sumdenum1990 <- tapply(casendata1990$expc,casendata1990$comu,FUN=sum,na.rm=T)
Rate1990 <- sumnum1990/sumdenum1990


casendata1992$Ipoor1992 <- ifelse(casendata1992$corte=="No pobre",0,1)  
sumnum1992 <- tapply(casendata1992$expc*casendata1992$Ipoor1992,casendata1992$comu,FUN=sum,na.rm=T)
sumdenum1992 <- tapply(casendata1992$expc,casendata1992$comu,FUN=sum,na.rm=T)
Rate1992 <- sumnum1992/sumdenum1992

casendata1994$Ipoor1994 <- ifelse(casendata1994$corte=="No pobres",0,1)  
sumnum1994 <- tapply(casendata1994$expc*casendata1994$Ipoor1994,casendata1994$comu,FUN=sum,na.rm=T)
sumdenum1994 <- tapply(casendata1994$expc,casendata1994$comu,FUN=sum,na.rm=T)
Rate1994 <- sumnum1994/sumdenum1994

casendata1996$Ipoor1996 <- ifelse(casendata1996$corte=="NO POBRE",0,1)  
sumnum1996 <- tapply(casendata1996$expc*casendata1996$Ipoor1996,casendata1996$comu,FUN=sum,na.rm=T)
sumdenum1996 <- tapply(casendata1996$expc,casendata1996$comu,FUN=sum,na.rm=T)
Rate1996 <- sumnum1996/sumdenum1996

casendata1998$Ipoor1998 <- ifelse(casendata1998$corte=="No Pobre",0,1)  
sumnum1998 <- tapply(casendata1998$expc*casendata1998$Ipoor1998,casendata1998$comu,FUN=sum,na.rm=T)
sumdenum1998 <- tapply(casendata1998$expc,casendata1998$comu,FUN=sum,na.rm=T)
Rate1998 <- sumnum1998/sumdenum1998

casendata2000$Ipoor2000 <- ifelse(casendata2000$corte=="No pobre",0,1)  
sumnum2000 <- tapply(casendata2000$expc*casendata2000$Ipoor2000,casendata2000$comu,FUN=sum,na.rm=T)
sumdenum2000 <- tapply(casendata2000$expc,casendata2000$comu,FUN=sum,na.rm=T)
Rate2000 <- sumnum2000/sumdenum2000

casendata2003$Ipoor2003 <- ifelse(casendata2003$CORTE=="No Pobre",0,1)  
sumnum2003 <- tapply(casendata2003$EXPC*casendata2003$Ipoor2003,casendata2003$COMU,FUN=sum,na.rm=T)
sumdenum2003 <- tapply(casendata2003$EXPC,casendata2003$COMU,FUN=sum,na.rm=T)
Rate2003 <- sumnum2003/sumdenum2003

casendata2006$Ipoor2006 <- ifelse(casendata2006$CORTE=="No pobre",0,1)  
sumnum2006 <- tapply(casendata2006$EXPC*casendata2006$Ipoor2006,casendata2006$COMUNA,FUN=sum,na.rm=T)
sumdenum2006 <- tapply(casendata2006$EXPC,casendata2006$COMUNA,FUN=sum,na.rm=T)
Rate2006 <- sumnum2006/sumdenum2006 

casendata2009$Ipoor2009 <- ifelse(casendata2009$CORTE=="No pobre",0,1)  
sumnum2009 <- tapply(casendata2009$EXPC*casendata2009$Ipoor2009,casendata2009$COMUNA,FUN=sum,na.rm=T)
sumdenum2009 <- tapply(casendata2009$EXPC,casendata2009$COMUNA,FUN=sum,na.rm=T)
Rate2009 <- sumnum2009/sumdenum2009 

casendata2011$Ipoor2011 <- ifelse(casendata2011$corte=="No pobres",0,1)  
sumnum2011 <- tapply(casendata2011$expc_r2*casendata2011$Ipoor2011,casendata2011$comuna,FUN=sum,na.rm=T)
sumdenum2011 <- tapply(casendata2011$expc_r2,casendata2011$comuna,FUN=sum,na.rm=T)
Rate2011 <- sumnum2011/sumdenum2011 

casendata2013$Ipoor2013 <- ifelse(casendata2013$pobreza_MN=="No pobres",0,1)  
sumnum2013 <- tapply(casendata2013$expc*casendata2013$Ipoor2013,casendata2013$comuna,FUN=sum,na.rm=T)
sumdenum2013 <- tapply(casendata2013$expc,casendata2013$comuna,FUN=sum,na.rm=T)
Rate2013 <- sumnum2013/sumdenum2013 
       
par(mfrow=c(2,2))
RateSANTIAGO <- as.vector(c(Rate1990["Santiago "],Rate1992["Santiago "],Rate1994["Santiago "],Rate1996["Santiago "],
Rate1998["Santiago "],Rate2000["SANTIAGO"],Rate2003["Santiago "],Rate2006["Santiago "],
Rate2009["Santiago"],Rate2011["Santiago"],Rate2013["Santiago"]))
plot(RateSANTIAGO,ylab="Poverty Rate",axes=FALSE, ann=FALSE)
axis(1,at=1:11,lab=c(1990,1992,1994,1996,1998,2000,2003,2006,2009,2011,2013))
axis(2,at=seq(0,0.4,0.2),ylab="Poverty Rate")
box()
lines(RateSANTIAGO, type="o", pch=22, lty=1, col="black")
title(main="Poverty Rates for SANTIAGO", col.main="black", font.main=4)
title(xlab="Year")
title(ylab="Poverty Rate")


RateLOLOL <- as.vector(c(Rate1990["Lolol  "],Rate1992["Lolol  "],Rate1994["Lolol  "],
Rate1996["Lolol  "],Rate1998["Lolol  "],Rate2000["LOLOL"],Rate2003["Lolol  "],Rate2006["Lolol "],
Rate2009["Lolol"],Rate2011["Lolol"],Rate2013["Lolol"]))
plot(RateLOLOL,ylab="Poverty Rate",axes=FALSE, ann=FALSE)
axis(1,at=1:11,lab=c(1990,1992,1994,1996,1998,2000,2003,2006,2009,2011,2013))
axis(2,at=seq(0,0.4,0.2),ylab="Poverty Rate")
box()
lines(RateLOLOL, type="o", pch=22, lty=1, col="black")
title(main="Poverty Rates for LOLOL", col.main="black", font.main=4)
title(xlab="Year")
title(ylab="Poverty Rate")


RateHUALPEN <- as.vector(c(Rate1990["Hualpen  "],Rate1992["Hualpen  "],Rate1994["Hualpen  "],
Rate1996["Hualpen  "],Rate1998["Hualpen  "],NA,Rate2003["Hualpen  "],Rate2006["Hualp?n "],
Rate2009["Hualp?n"],Rate2011["Hualp?n"],Rate2013["Hualp?n"]))
plot(RateHUALPEN,ylab="Poverty Rate",axes=FALSE, ann=FALSE)
axis(1,at=1:11,lab=c(1990,1992,1994,1996,1998,2000,2003,2006,2009,2011,2013))
axis(2,at=seq(0,0.4,0.2),ylab="Poverty Rate")
box()
lines(RateHUALPEN, type="o", pch=22, lty=1, col="black")
title(main="Poverty Rates for HUALPEN", col.main="black", font.main=4)
title(xlab="Year")
title(ylab="Poverty Rate")

RateCONCON <- as.vector(c(Rate1990["Conc?n "],Rate1992["Conc?n "],Rate1994["Conc?n "],
Rate1996["Conc?n "],Rate1998["Conc?n "],Rate2000["CON -CON"],Rate2003["Conc?n "],
Rate2006["Conc?n "],Rate2009["Conc?n"],Rate2011["Conc?n"],Rate2013["Conc?n"]))
plot(RateCONCON,ylab="Poverty Rate",axes=FALSE, ann=FALSE)
axis(1,at=1:11,lab=c(1990,1992,1994,1996,1998,2000,2003,2006,2009,2011,2013))
axis(2,at=seq(0,0.4,0.2),ylab="Poverty Rate")
box()
lines(RateCONCON, type="o", pch=22, lty=1, col="black")
title(main="Poverty Rates for CONCON", col.main="black", font.main=4)
title(xlab="Year")
title(ylab="Poverty Rate")


  ##matplot
par(mfrow=c(1,1))
act <- c(2006,2009,2011,2013)
Rates <- c(Rate2006["Santiago "],
           Rate2009["Santiago"],Rate2011["Santiago"],Rate2013["Santiago"],
           Rate2006["Lolol "],
           Rate2009["Lolol"],Rate2011["Lolol"],Rate2013["Lolol"],
           Rate2006["Hualp?n "],
           Rate2009["Hualp?n"],Rate2011["Hualp?n"],Rate2013["Hualp?n"],
           Rate2006["Conc?n "],Rate2009["Conc?n"],Rate2011["Conc?n"],Rate2013["Conc?n"])
rlabels = c("2006", "2009", "2011", "2013")
clabels = c("SANTIAGO", "LOLOL", "HUALPEN", "CONCON")
Res <- matrix(Rates,nrow=4,ncol=4,byrow=T)

matplot(act, t(Res), type ="b", pch=18, lty=c(1,1,1,1),lwd = 2,xlab="Year",ylab="Poverty rate",
        main="Time series direct poverty rate estimates for a four selected comunas",col=c("orange","red","green","blue"))
legend(2006,0.35, lwd = 2, legend = clabels,col=c("orange","red","green","blue"))

  ##Confidence Interval Plot for all COMUNAs
  ##Look at program CORTEv.R (need to run
  ##the program from there from line1--line160)
COMUNAFRAME2 <- data.frame(RatioCOMUNA,sqrt(VarDesignCOMUNA),PSUinCOMUNA_woNA)

LOWER <- COMUNAFRAME2$RatioCOMUNA-
  1.96*COMUNAFRAME2$sqrt.VarDesignCOMUNA.

COMUNAFRAME2$LOWERBOUND <- ifelse(LOWER<0,0,LOWER)

COMUNAFRAME2$UPPERBOUND <- COMUNAFRAME2$RatioCOMUNA+
  1.96*COMUNAFRAME2$sqrt.VarDesignCOMUNA.


UPDATECOMUNA1 <- COMUNAFRAME2[order(COMUNAFRAME2$PSUinCOMUNA_woNA),]

UPDATECOMUNA2 <- COMUNAFRAME2[order(COMUNAFRAME2$sqrt.VarDesignCOMUNA.),]

UPDATECOMUNA3 <- COMUNAFRAME2[order(COMUNAFRAME2$RatioCOMUNA),]

plot(UPDATECOMUNA1$RatioCOMUNA,type="p",ylab="Poverty Rate",
     xlab="COMUNAs",main="Confidence Interval for Poverty Rate of COMUNAs in 
     CASEN 2009 Sorted by the Number of PSUs in COMUNAs",lwd=2.5)
lines(UPDATECOMUNA1$LOWERBOUND,lty=1,col="cadetblue2",cex=5)
lines(UPDATECOMUNA1$UPPERBOUND,lty=1,col="cadetblue2",cex=5)


plot(UPDATECOMUNA2$RatioCOMUNA,type="p",ylab="Poverty rate",
     xlab="Comuna",main="Direct poverty rates and associated direct confidence intervals 
for Chilean comunas; 
     comunas sorted by the direct variance estimates",col="sienna",lwd=4)
lines(UPDATECOMUNA2$LOWERBOUND,lty=1,col="cadetblue2",cex=5)
lines(UPDATECOMUNA2$UPPERBOUND,lty=1,col="cadetblue2",cex=5)


plot(UPDATECOMUNA2$RatioCOMUNA,type="p",ylab="Poverty Rate",
               xlab="COMUNAs",main="Confidence Interval for Poverty Rate of COMUNAs in 
               CASEN 2009 Sorted by the Variance",col="sienna",lwd=2)
dgrid <- seq(1,334)
polygon(x=c(dgrid, rev(dgrid)),
        y=c(UPDATECOMUNA2$UPPERBOUND,
            UPDATECOMUNA2$LOWERBOUND),
        col="lightgray", border=NA)



plot(UPDATECOMUNA3$RatioCOMUNA,type="p",ylab="Poverty Rate",
               xlab="COMUNAs",main="Confidence Interval for Poverty Rate of COMUNAs in 
               CASEN 2009 Sorted by the Poverty Rate",col="sienna",lwd=2)
lines(UPDATECOMUNA3$LOWERBOUND,lty=1,col="cadetblue2",cex=5)
lines(UPDATECOMUNA3$UPPERBOUND,lty=1,col="cadetblue2",cex=5)


scatter.smooth(UPDATECOMUNA3$RatioCOMUNA,type="p",ylab="Poverty Rate",
               xlab="COMUNAs",main="Confidence Interval for Poverty Rate of COMUNAs in 
               CASEN 2009 Sorted by the Poverty Rate",col="sienna",lwd=2)
dgrid <- seq(1,334)
polygon(x=c(dgrid, rev(dgrid)),
        y=c(UPDATECOMUNA3$UPPERBOUND,
            UPDATECOMUNA3$LOWERBOUND),
        col="lightgray", border=NA)





