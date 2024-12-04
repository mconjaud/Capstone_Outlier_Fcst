
###################################################################
## AULA 1: CONCEITOS DE SÉRIES DE TEMPO
## Autor: Paloma Vaissman Uribe
## Material desenvolvido para o curso do PADS - Financial Analytics
## Abr 2022
###################################################################

# Motivação
library(readr)
us <- read_csv("https://raw.githubusercontent.com/mconjaud/Capstone/refs/heads/main/preco_anbima_di_spread.csv?token=GHSAT0AAAAAAC3H2RWQJHTHPEQJJHLEZD74Z2MXBIA")
plot(us$date,us$cases,type='l',ylab='cases',xlab='date')

# correlação entre lag 1 e nível
library(dplyr)
us <- us %>% mutate(us_lag=lag(cases,2))
cor(us$cases,us$us_lag,use='na.or.complete')

# matriz de correlação
library(reshape2)
library(ggplot2)
matrix <- matrix(NA,length(us$cases),61)
matrix[,1] <- us$cases
for (i in 2:61){
  matrix[,i] <- lag(us$cases,i)
}
cormat <- cor(matrix,use='na.or.complete')
melted_cormat <- melt(cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Simulando um processo iid
X <- rnorm(1000,mean=0,sd=2)
par(mfrow=c(1,2))
ts.plot(X)
acf(X,main="")

# Processo AR(1) - estacionário
phi <- 0.5
n <- 1000
e <- rnorm(n,mean=0,sd=2)
y <- rep(NA,n)
y[1] <- e[1]
for (i in 2:n){
  y[i] <- phi*y[i-1] + e[i]
}

# Processo passeio aleatório
phi <- 1
n <- 1000
e <- rnorm(n,mean=0,sd=2)
y <- rep(NA,n)
y[1] <- e[1]
for (i in 2:n){
  y[i] <- phi*y[i-1] + e[i]
}

########## COMPONENTES DE UMA SÉRIE ###################

# Decomposição de séries de tempo usando R basis
data(AirPassengers)
AP <- AirPassengers
decomposeAP <- decompose(AP,"multiplicative")
plot(decomposeAP)

# Outros pacotes: forecast
library(forecast)
library(ggplot2)
data(AirPassengers)
AP <- AirPassengers
fit2 <- ets(AP)
autoplot(fit2)

# Outros pacotes: prophet
library(prophet)
library(xts)
data(AirPassengers)
AP <- AirPassengers
df <- data.frame(ds = index(as.xts(AP)),y = coredata(as.xts(AP)))
m <- prophet(df)
future <- make_future_dataframe(m, periods = 12)
forecast = predict(m,future)
prophet_plot_components(m,forecast)

########## TESTE DE RAIZ UNITARIA ###################

# Teste ADF
library(tseries)
library(readr)
GDP_EUA <- read_csv("Documents/Insper/PADS/FA_2022_aula1/GDP_EUA.csv") #trocar caminho
gdp = GDP_EUA$NA000334Q 
adf.test(gdp, k = trunc((length(gdp)-1)^(1/3))) 
# MAIS INFO: ver https://robjhyndman.com/eindhoven/2-3-Differencing.pdf

# teste KPSS: pacote fable
library(fpp3)
library(fable)
us_employment %>%
  filter(Title=='Total Private') %>%
  features(Employed, unitroot_kpss)

########## MODELAGEM ARIMA e BOX-JENKINS ###################

# Simulando um processo ARIMA
set.seed(123)
y <- arima.sim(n=200,list(order=c(1,1,0),ar=0.8),sd=1)
par(mfrow=c(2,3))
ts.plot(y)
acf(y,main="",60)
acf(y,type="partial",main="")
dy <- diff(y)
ts.plot(dy,main="")
acf(dy,60,main="")
acf(dy,type="partial",main="")

# Simulating a IMA(1,1) process
set.seed(123)
ima1 <- arima.sim(list(order = c(0,1,1), ma = 0.9), n = 100)

# Box-Jenkins: identificação das ordens
par(mfrow=c(2,3))
ts.plot(ima1)
acf(ima1,main="")
acf(ima1,type="partial",main="")
ts.plot(diff(ima1))
acf(diff(ima1),main="")
acf(diff(ima1),type="partial",main="")

# Box-Jenkins: estimação
library(forecast)
Arima(ima1,order=c(0,1,1))

# alternativa para identificação visual: looping para minimizar ICs
p <- 3
q <- 3

resultsAIC <- matrix(0,p+1,q+1)
colnames(resultsAIC) <- c(0:p)
rownames(resultsAIC) <- c(0:q)

for (i in 0:p){
  for(j in 0:q){
    resultsAIC[i+1,j+1] <- Arima(ima1, order = c(i, 1, j))$aicc
  }
}
which(resultsAIC == min(resultsAIC), arr.ind = TRUE)

# Box-Jenkins: diagnóstico
res <- Arima(ima1,order=c(0,1,1))$res
par(mfrow=c(1,1))
Acf(res,main="",xlab="")
Box.test(res, lag=24, fitdf=4, type="Ljung")

# previsão
fit <- Arima(ima1,order=c(0,1,1))
plot(forecast(fit,h=12))

########## MODELAGEM DA SAZONALIDADE ###################

# prophet
library(prophet)
library(xts)
data(AirPassengers)
AP <- AirPassengers
df <- data.frame(ds = zoo::index(as.xts(AP)),y = coredata(as.xts(AP)))
m <- prophet(df, seasonality.mode = "multiplicative") # ALTERNATIVA: sazonalidade multiplicativa (default é aditiva)
future <- make_future_dataframe(m, 12, freq = 'm')
forecast <- predict(m, future)
plot(m, forecast)

# SARIMA
## analise visual
y <- df$y # aproveitando df do prophet
par(mfrow=c(2,2))
acf(y,main="Nível",96)
acf(diff(y),main="Primeira diferença",96)
acf(diff(diff(y),12),96,main='',ylab='ACF diff diff')
pacf(diff(diff(y),12),96,main='',ylab='PACF diff diff')

# testando vários modelos
library(forecast)
fit1 = Arima(y,order=c(1,1,0),seasonal=list(order=c(1,1,0),period=12),method="ML")
fit2 = Arima(y,order=c(1,1,1),seasonal=list(order=c(1,1,0),period=12),method="ML")
fit3 = Arima(y,order=c(0,1,0),seasonal=list(order=c(0,1,0),period=12),method="ML")
fit4 = Arima(y,order=c(0,1,1),seasonal=list(order=c(0,1,0),period=12),method="ML")
fit5 = Arima(y,order=c(1,1,0),seasonal=list(order=c(1,1,1),period=12),method="ML")
fit6 = Arima(y,order=c(1,1,0),seasonal=list(order=c(0,1,0),period=12),method="ML")
aic  = c(fit1$aic,fit2$aic,fit3$aic,fit4$aic,fit5$aic,fit6$aic)
aicc = c(fit1$aicc,fit2$aicc,fit3$aicc,fit4$aicc,fit5$aicc,fit6$aicc)
bic  = c(fit1$bic,fit2$bic,fit3$bic,fit4$bic,fit5$bic,fit6$aicc)
cbind(aic,aicc,bic)

# melhor modelo: ajuste, diagnóstico e previsão
fit = Arima(y,order=c(1,1,0),seasonal=list(order=c(0,1,0),period=12),method="ML")
par(mfrow=c(1,3))
acf(residuals(fit),96,main="")
pacf(residuals(fit),96,main="")
plot(forecast(fit,12))

# usando auto.arima: lembrar de colocar o input como objeto ts
library(forecast)
auto.arima(AP)

# usando fable
library(fable)
airline_fable <- as_tsibble(AP) %>%
  model(stepwise = ARIMA(value))
glance(airline_fable)

# Ajustando um modelo linear aditivo (dummies sazonais + tendência) 
library(lubridate)
months <- as.factor(month(df$ds)) #extracting month using lubridate
dummies <- as.data.frame(model.matrix(~months))
colnames(dummies) <- c("intercept","fev","mar","abr","mai","jun",
                       "jul","ago","set","out","nov","dez")
dummies$trend <- c(1:dim(dummies)[1]) # ading a deterministic trend
data.frame <- cbind(y,dummies)
attach(data.frame)
summary(lm(y~fev+mar+abr+mai+jun+jul+ago+set+out+nov+dez+trend))

# usando um ARIMAX: variáveis exógenas
fit = Arima(y[1:132],order=c(1,1,0),method="ML",xreg=as.matrix(data.frame[1:132,3:14]))
par(mfrow=c(1,1))
plot(forecast(fit,12,xreg=as.matrix(dummies[133:144,2:13])))
