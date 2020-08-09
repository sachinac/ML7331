library(tswge)
library(orcutt)
library(vars)
library(nnfor)
library(forecast)
library(dplyr)

uri.data = 'https://covidtracking.com/api/v1/states/tx/daily.csv'
texas = read.csv(uri.data, header = T)
TX <- texas %>% filter(date > 20200618 & date < 20200707)
TX_Future <- texas %>% filter(date > 20200706) 
TX$Curve_Day <- seq(1, length(TX$date), by= 1)
TX_Future$Curve_Day <- seq(length(TX$date)+1, (length(TX$date) + length(TX_Future$date)), by= 1)

# National to get population data ############
uri.data = 'https://covidtracking.com/api/v1/states/ny/daily.csv'
NY = read.csv(uri.data, header = T)
NY <- NY %>% filter(date > 20200319 & date < 20200506)
NY$Curve_Day <- seq(1, length(NY$date), by= 1)

# National to get population data ############
uri.data = 'https://covidtracking.com/api/v1/states/nj/daily.csv'
NJ = read.csv(uri.data, header = T)

NJ <- NJ %>% filter(date > 20200425 & date < 20200519)
NJ$Curve_Day <- seq(1, length(NJ$date), by= 1)

NY <- NY%>% mutate(positive_perc = (positiveIncrease/totalTestResults)*100)
NJ <- NJ%>% mutate(positive_perc = (positiveIncrease/totalTestResults)*100)

####### Future texas Data  ########
FutureData = read.csv("~/Documents/datascience/DS7331/CoronaCurves-master/Corona_MAE.csv", header = TRUE)
texasPos = FutureData$TX.New.Cases
lenTXFuture = length(texasPos)



######## NY VAR Model 
count = length(NY$positiveIncrease)
CoroVar = VAR(cbind(NY$positiveIncrease * (8/29), NY$negativeIncrease * (8/29)), type = "both", lag.max = 10)

preds = predict(CoroVar,n.ahead = lenTXFuture) 

ASE = mean((texasPos * (8/29) - preds$fcst$y1[,1])^2)
ASE

MAE = sqrt(ASE)

VAR_NY_DF = data.frame(State = "New York",Model = 'VAR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

plot(preds)

######## NJ VAR Model 
count = length(NJ$positiveIncrease)
CoroVar = VAR(cbind(NJ$positiveIncrease * (8/29), NJ$negativeIncrease * (8/29)), type = "both", lag.max = 10)

preds = predict(CoroVar,n.ahead = lenTXFuture) 

ASE = mean((texasPos * (8/29) - preds$fcst$y1[,1])^2)
ASE

MAE = sqrt(ASE)

VAR_NJ_DF = data.frame(State = "New Jersey", Model = 'VAR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

plot(preds)


#### Ensemble/Texas Model

preds_TX = (preds_NY$fcst$y1[,1] + preds_NJ$fcst$y1[,1])/2

ASE = mean((texasPos  * (8/29) - preds_TX)^2)
MAE = sqrt(ASE)

VAR_TX_DF = data.frame(State = "Ensemble/Texas",Model = 'VAR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

df = rbind(VAR_NJ_DF,VAR_NY_DF,VAR_TX_DF)
df
