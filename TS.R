#Time-Series Modeling and Analysis of 
#Monthly arrivals of tourists coming to the United States of America.
#Term-Project for Time-Series Analysis, Modeling and Control - Spring '17
#Swarup Sahoo, Chih-Hsiang Hsieh, Xintong Zhou 
#Operations Research and Industrial Engineering
#The University of Texas at Austin

library(fpp)
library(readr)
library(MTS)
library(forecast)

#Importing data and stablizing-------------------
UST <- read_csv("UST.csv")

tot.tr<-UST[-c(167:176),1]
tot.te<-UST[c(167:176),1]
tsdisplay(tot.tr)     #exploding variance
tot.tr<-log(tot.tr)   #to compensate for the increasing variance
tsdisplay(tot.tr)        #variance is more in control now

#conversion to time-series from 2002 for training and 2016 for testing---
tot.tr<-ts(tot.tr[[1]], frequency=12, start=c(2002,1))
tot.te<-ts(tot.te[[1]], frequency = 12, start = c(2015,11))

#Incremental model Approach ARMA(2n,2n-1)--------
N<-length(tot.tr)
flag1=0
n=1
model1<-Arima(tot.tr, order = c(2*n,0,2*n-1), method="ML")
RSS.1<-sum((tot.tr-fitted(model1))^2)
while(flag1==0){
  model2<-Arima(tot.tr, order = c(2*(n+1),0,2*(n+1)-1), method = "ML")
  RSS.2<-sum((tot.tr-fitted(model2))^2)
  F.model<-((RSS.1-RSS.2)/4)/(RSS.2/N)
  if(F.model>qf(.95, df1=4, df2=N)){
    print(F.model)
    print(qf(.95, df1=4, df2=N))
    RSS.1<-RSS.2
    model1<-model2
    n=n+1
  }
  else{
    print(F.model)
    print(qf(.95, df1=4, df2=N))
    flag1=1;
  }
}

summary(model1) #We get an ARMA 14,13 model (n=7)
plot(model1)

#Forecast for ARMA(14,13)------------------------
fct.arma<-forecast(model1, h=12, level=95)
fct.arma$mean<-exp(fct.arma$mean)
fct.arma$upper<-exp(fct.arma$upper)
fct.arma$lower<-exp(fct.arma$lower)
fct.arma$x<-exp(fct.arma$x)
plot(fct.arma, lwd=2, ylab="#of Tourists")
points(tot.te)
legend(2005,4.5e+06, c("Forecast", "Actual Data"),pch=c("-","o"),col=c("blue", "black"))

coef.I<-data.frame(coef(model1))
phi<-coef.I[c(1:(2*n)),1]
theta<-coef.I[c((2*n+1):(4*n-1)),1]
ar.roots<-polyroot(c(-coef.I[c((2*n):1),1],1))

#Plot of AR and MA Roots-------------------------
plot(model1)


#Formula to calculate g (Strength) for Greens Function
g=complex(real = c(1:(2*n)), imaginary = c(1:(2*n))) #Strength g values for Greens Function
for (i in 1:(2*n)){
  numerator=ar.roots[i]^(2*n-1);
  for(j in 1:(2*n-1)){
    numerator=numerator-theta[j]*ar.roots[i]^(2*n-1-j);
  }
  denominator=1;
  for(j in 1:(2*n)){
    if(i==j) next
    denominator=denominator*(ar.roots[i]-ar.roots[j]);
  }
  g[i]=numerator/denominator;
}


A=2*abs(g) #Amplitute
RtF=A/max(A) #Ratio to Fundamental

#Trend Calculation-------------------------------
ar.roots
#We have two real roots close to 1, hence linear trend exists
#We will use (1-B)^2 as the trend operator

Xt<-tot.tr
#Removing trend and seasonality------------------
#Yt<-lag(Xt,8) + 0.977*lag(Xt,7) + 1.00169*lag(Xt,6) + 0.0103161*lag(Xt,5) + 0.0298983*lag(Xt,4) + 0.0103161*lag(Xt,3) + 1.00169*lag(Xt,2) + 0.977*lag(Xt,1) + Xt
#Yt<-(-lag(Xt,9)) + 0.023*lag(Xt,8)- 0.024685 *lag(Xt,7) + 0.991369*lag(Xt,6) - 0.0195822*lag(Xt,5) + 0.0195822 *lag(Xt,4) - 0.991369 *lag(Xt,3) + 0.024685 *lag(Xt,2) - 0.023 *lag(Xt,1) + Xt
Yt<-(-lag(Xt,10)) - 0.977*lag(Xt,9) - 0.001685*lag(Xt,8) + 0.966684*lag(Xt,7) + 0.971787*lag(Xt,6) - 0.971787 *lag(Xt,4) - 0.966684 *lag(Xt,3) + 0.001685*lag(Xt,2) + 0.977*lag(Xt,1) + Xt
#Calculating parsimonious ARMA model-------------
Yt.arma<-Arima(Yt,order=c(4,0,13))
RSS.parsi<-sum(residuals(Yt.arma)^2)

coef.Yt<-data.frame(coef(Yt.arma))
phi.Yt<-coef.Yt[c(1:4),1]
ar.roots.Yt<-polyroot(c(-coef.Yt[c(4:1),1],1))


#Non-Stationary Analysis-------------------------
t<-c(1:length(Xt))
fit.lm=lm(Xt~t+sin(2*pi/12*t)+cos(2*pi/12*t))
yt=residuals(fit.lm)
yt=ts(yt,frequency = 12, start=c(2002,1))
summary(fit.lm)

#Incremental model Approach ARMA(2n,2n-1)--------
N<-length(yt)       
flag1=0
n=1
model.yt1<-Arima(yt, order = c(2*n,0,2*n-1))
RSS.yt1<-sum(residuals(model.yt1)^2)
while(flag1==0){
  model.yt2<-Arima(yt, order = c(2*(n+1),0,2*(n+1)-1));
  RSS.yt2<-sum((yt-fitted(model.yt2))^2);
  F.model<-((RSS.yt1-RSS.yt2)/4)/(RSS.yt2/N);
  if(F.model>qf(.95, df1=4, df2=N)){
    print(F.model)
    print(qf(.95, df1=4, df2=N))
    RSS.yt1<-RSS.yt2
    model.yt1<-model.yt2
    n=n+1;
  }
  else{
    print(F.model)
    print(qf(.95, df1=4, df2=N))
    flag1=1;
  }
}

summary(model.yt1)

newdata = data.frame(t=c(167:176))
fct.fit<-predict(fit.lm, newdata, h=12)
fct.ns<-forecast(model.yt1, h=12, level=95)
fct.ns$mean<-exp(fct.ns$mean+fct.fit)
fct.ns$upper<-exp(fct.ns$upper+fct.fit)
fct.ns$lower<-exp(fct.ns$lower+fct.fit)
fct.ns$x<-exp(fct.ns$x+fitted(fit.lm))
plot(fct.ns, n1=c(2014,1),lwd=2, main="Non-Stationary ARMA(14,13)", ylab="#of tourists")
points(tot.te)
legend(2014,4.5e+06, c("Forecast", "Actual Data"),pch=c("-","o"),col=c("blue", "black"))


#ARMAV incremental approach ARMAV(2n,2n-1)-------

k=4 #number of vectors
UST <- read_csv("UST.csv")
UST.L<-log(UST[-c(167:176),])
names(UST.L)[1]<-paste("TOTAL")
UST.V<-UST.L[,c("TOTAL","GERMANY","SOUTH KOREA", "AUSTRALIA")]
UST.V<-ts(UST.V, frequency = 12, start = c(2002,1))
#Lets start with ARMA(2n,2n-1)
flag1=0
n=1
modelV<-VARMA(UST.V,p=2*n,q=2*n-1)
RSS.1V<-sum(residuals(modelV)^2)
while(flag1==0){
  modelV2<-VARMA(UST.V,p=2*(n+1),q=2*(n+1)-1)
  RSS.2V<-sum(residuals(modelV2)^2)
  F.modelV<-((RSS.1V-RSS.2V)/(4*k*k))/(RSS.2V/(N*k))
  if(F.modelV>qf(.95, df1=(4*k*k), df2=N*k)){
    print(F.modelV)
    print(qf(.95, df1=(4*k*k), df2=N*k))
    modelV<-modelV2
    RSS.1V<-RSS.2V
    n=n+1
  }
  else{
    print(F.modelV)
    print(qf(.95, df1=(4*k*k), df2=N*k))
    flag1=1
  }
}

vpred<-VARMApred(modelV, h=12)
vpred$pred<-exp(vpred$pred)
vpred$pred
RSS.armav<-sum(residuals(modelV)^2)
