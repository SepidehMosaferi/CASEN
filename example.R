data("MDarea.pop")
attach(MDarea.pop)
  ######Assessing the Multistage Design Formula############
Ni <- table(MDarea.pop$TRACT)
m <- 20
probi <- m*Ni/sum(Ni)

set.seed(-780087528)
sam <- cluster(data=MDarea.pop,clustername = "TRACT",
               size=m,method="systematic",
               pik=probi,description = TRUE)
samclus <- getdata(MDarea.pop,sam)
samclus <- rename(samclus,c(Prob="pi1"))
s <- strata(data=as.data.frame(samclus),
            stratanames = "TRACT",
            size=rep(50,m),method="srswor")
samdat <- getdata(samclus,s)
samdat <- rename(samdat,c(Prob="pi2"))
pick <- names(Ni) %in% sort(unique(samdat$TRACT))
Ni.sam <- Ni[pick]
d1 <- Ni.sam/sum(Ni)
wt <- 1/samdat$pi1/samdat$pi2
BW2stagePPSe(Ni=Ni.sam,ni=rep(50,20),X=samdat$y1,
             psuID=samdat$TRACT,w=wt,
             m=20,pp=d1)

  #Based on the model formula: [1+(bbar-1)rho]sigma^2/a
bbar <- 1000/20
rho <- 0.0029149
sigma2 <- sd(samdat$y1)^2
a <- length(unique(samdat$TRACT))
(sigma2/a)*(1+((bbar-1)*rho))  ##answer

  #Based on the ultimate cluster estimator (manually computed)
b <- unique(tapply(samdat$y1,samdat$TRACT,FUN=function(x){length(x)})); b  
A <- length(unique(MDarea.pop$TRACT)) ; A 
B <- mean(table(MDarea.pop$TRACT)); B
ybar_alpha <- tapply(samdat$y1,samdat$TRACT,mean); ybar_alpha
ybar <- mean(samdat$y1); ybar
s2a <- sum((ybar_alpha-ybar)^2)/(a-1); s2a  
(1-((a*b)/(A*B)))*(s2a/a) ##answer (**1)
s2a/a

  #Based on the formula of HHT (p. 144)
DATA <- subset(MDarea.pop,MDarea.pop$TRACT %in% samdat$TRACT)
Bi <- tapply(DATA$y1,DATA$TRACT,FUN=function(x){length(x)})
bi <- rep(50,20)
ybari. <- tapply(samdat$y1,samdat$TRACT,mean); ybari.
 
s2b <- (1/(bi-1))*tapply(samdat$y1,samdat$TRACT,sd)
s2b_2 <- s2b^2
(A*(A-a)/a)*s2a+(A/a)*sum((Bi/bi)*(Bi-bi)*s2b)  ##answer

  #surveydesign
options(survey.lonely.psu="adjust")
samdat$vazn <- wt
DESIGN<-svydesign(id=~TRACT, weight=~vazn, data=samdat)
svymean(samdat$y1,DESIGN)  ##answer (**2)

DESIGN

#*************************************************************************
  ##Conclusion:
##The svyesign uses the s2a/a for the ultimate cluster estimator based on 
##the results from the (**1) and (**2).
#*************************************************************************


######Assessing the Stratified Multistage Design Formula############
  #Take a sample separately in each stratum and assume BLKGROUP is
  #the stratum identification

STRATDATA <- MDarea.pop[order(BLKGROUP,TRACT),]
table(BLKGROUP) #we have 6 strata

POPSTRAT1 <- subset(STRATDATA,BLKGROUP==1) #stratum1
POPSTRAT2 <- subset(STRATDATA,BLKGROUP==2) #stratum2
POPSTRAT3 <- subset(STRATDATA,BLKGROUP==3) #stratum3
POPSTRAT4 <- subset(STRATDATA,BLKGROUP==4) #stratum4
POPSTRAT5 <- subset(STRATDATA,BLKGROUP==5) #stratum5
POPSTRAT6 <- subset(STRATDATA,BLKGROUP==6) #stratum6

N1 <- table(POPSTRAT1$TRACT); N2 <- table(POPSTRAT2$TRACT); N3 <- table(POPSTRAT3$TRACT)
N4 <- table(POPSTRAT4$TRACT); N5 <- table(POPSTRAT5$TRACT); N6 <- table(POPSTRAT6$TRACT)
m1 <- m2 <- m3 <- m4 <- m5 <- 10; m6 <- 1
prob1 <- m1*N1/sum(N1); prob2 <- m2*N2/sum(N2); prob3 <- m3*N3/sum(N3)
prob4 <- m4*N4/sum(N4); prob5 <- m5*N5/sum(N5); prob6 <- m6*N6/sum(N6)

  ##since the BLKGROUP6 just contains 1 PSU, we excludes it from
  ##the variance estimation since it is certainty PSU.
set.seed(-780087528)
sam1 <- cluster(data=POPSTRAT1,clustername = "TRACT",
               size=m1,method="systematic",
               pik=prob1,description = TRUE)
samclus1 <- getdata(POPSTRAT1,sam1)
samclus1 <- rename(samclus1,c(Prob="pi1"))
s1 <- strata(data=as.data.frame(samclus1),
            stratanames = "TRACT",
            size=rep(50,m1),method="srswor")
samdat1 <- getdata(samclus1,s1)
samdat1 <- rename(samdat1,c(Prob="pi2"))
pick1 <- names(N1) %in% sort(unique(samdat1$TRACT))
N1.sam <- N1[pick]
d11 <- N1.sam/sum(N1)
wt1 <- 1/samdat1$pi1/samdat1$pi2


set.seed(-780087528)
sam2 <- cluster(data=POPSTRAT2,clustername = "TRACT",
                size=m2,method="systematic",
                pik=prob2,description = TRUE)
samclus2 <- getdata(POPSTRAT2,sam2)
samclus2 <- rename(samclus2,c(Prob="pi1"))
s2 <- strata(data=as.data.frame(samclus2),
             stratanames = "TRACT",
             size=rep(50,m2),method="srswor")
samdat2 <- getdata(samclus2,s2)
samdat2 <- rename(samdat2,c(Prob="pi2"))
pick2 <- names(N2) %in% sort(unique(samdat2$TRACT))
N2.sam <- N2[pick]
d12 <- N2.sam/sum(N2)
wt2 <- 1/samdat2$pi1/samdat2$pi2


set.seed(-780087528)
sam3 <- cluster(data=POPSTRAT3,clustername = "TRACT",
                size=m3,method="systematic",
                pik=prob3,description = TRUE)
samclus3 <- getdata(POPSTRAT3,sam3)
samclus3 <- rename(samclus3,c(Prob="pi1"))
s3 <- strata(data=as.data.frame(samclus3),
             stratanames = "TRACT",
             size=rep(50,m3),method="srswor")
samdat3 <- getdata(samclus3,s3)
samdat3 <- rename(samdat3,c(Prob="pi2"))
pick3 <- names(N3) %in% sort(unique(samdat3$TRACT))
N3.sam <- N3[pick]
d13 <- N3.sam/sum(N3)
wt3 <- 1/samdat3$pi1/samdat3$pi2


set.seed(-780087528)
sam4 <- cluster(data=POPSTRAT4,clustername = "TRACT",
                size=m4,method="systematic",
                pik=prob4,description = TRUE)
samclus4 <- getdata(POPSTRAT4,sam4)
samclus4 <- rename(samclus4,c(Prob="pi1"))
s4 <- strata(data=as.data.frame(samclus4),
             stratanames = "TRACT",
             size=rep(50,m4),method="srswor")
samdat4 <- getdata(samclus4,s4)
samdat4 <- rename(samdat4,c(Prob="pi2"))
pick4 <- names(N4) %in% sort(unique(samdat4$TRACT))
N4.sam <- N4[pick]
d14 <- N4.sam/sum(N4)
wt4 <- 1/samdat4$pi1/samdat4$pi2


set.seed(-780087528)
sam5 <- cluster(data=POPSTRAT5,clustername = "TRACT",
                size=m5,method="systematic",
                pik=prob5,description = TRUE)
samclus5 <- getdata(POPSTRAT5,sam5)
samclus5 <- rename(samclus5,c(Prob="pi1"))
s5 <- strata(data=as.data.frame(samclus5),
             stratanames = "TRACT",
             size=rep(50,m5),method="srswor")
samdat5 <- getdata(samclus5,s5)
samdat5 <- rename(samdat5,c(Prob="pi2"))
pick5 <- names(N5) %in% sort(unique(samdat5$TRACT))
N5.sam <- N5[pick]
d15 <- N5.sam/sum(N5)
wt5 <- 1/samdat5$pi1/samdat5$pi2

  ##applying the svydesign individualy:
options(survey.lonely.psu = "adjust")
samdat1$vazn1 <- wt1
DESIGN1 <- svydesign(id=~TRACT,weight=~vazn1,data=samdat1)
svymean(samdat1$y1,DESIGN1)

options(survey.lonely.psu = "adjust")
samdat2$vazn2 <- wt2
DESIGN2 <- svydesign(id=~TRACT,weight=~vazn2,data=samdat2)
svymean(samdat2$y1,DESIGN2)

options(survey.lonely.psu = "adjust")
samdat3$vazn3 <- wt3
DESIGN3 <- svydesign(id=~TRACT,weight=~vazn3,data=samdat3)
svymean(samdat3$y1,DESIGN3)

options(survey.lonely.psu = "adjust")
samdat4$vazn4 <- wt4
DESIGN4 <- svydesign(id=~TRACT,weight=~vazn4,data=samdat4)
svymean(samdat4$y1,DESIGN4)

options(survey.lonely.psu = "adjust")
samdat5$vazn5 <- wt5
DESIGN5 <- svydesign(id=~TRACT,weight=~vazn5,data=samdat5)
svymean(samdat5$y1,DESIGN5)

  ##It seems that the stratum5 has some problems. So, we just
  ##work on the first 4 strata.
##Final answer (**)

NN <- length(N1)+length(N2)+length(N3)+length(N4)

(length(N1)/(NN))^2*(6.5292^2)+(length(N2)/(NN))^2*(7.2354^2)+
  (length(N3)/(NN))^2*(3.5835^2)+(length(N4)/(NN))^2*(3.1434^2) #(***1)


m1*(6.5292^2)+m2*(7.2354^2)+
  m3*(3.5835^2)+m4*(3.1434^2) #(***2)


  ##Based on just introducing a stratified 2stage design:
vazn <- c(wt1,wt2,wt3,wt4)
y1 <- c(samdat1$y1,samdat2$y1,samdat3$y1,samdat4$y1)
TRACT <- c(samdat1$TRACT,samdat2$TRACT,samdat3$TRACT,samdat4$TRACT)
BLKGROUP <- c(samdat1$BLKGROUP,samdat2$BLKGROUP,samdat3$BLKGROUP,samdat4$BLKGROUP)
samdatFINAL <- data.frame(cbind(y1,vazn,BLKGROUP,TRACT))
samdatFINAL <- samdatFINAL[order(BLKGROUP,TRACT),]

options(survey.lonely.psu="adjust")
DESIGNFINAL <- svydesign(id=~TRACT,strata=BLKGROUP,weight=~vazn,nest=T,data=samdatFINAL)
svymean(samdatFINAL$y1,DESIGNFINAL)
3.3122^2 #(***3)

#****************************************************
#####
###conclusion: the results based on the (***1) and (****3) are the same 
###in the method of stratifies multistage sampling variance
#***************************************************





