library(tswge)
library(orcutt)
library(vars)
library(nnfor)
library(forecast)
library(dplyr)

uri.data = 'https://covidtracking.com/api/v1/states/tx/daily.csv'
texas = read.csv(uri.data, header = T)

#####################################################
usCovidClean <- texas %>% replace(is.na(.), 0)
DataUsing <- usCovidClean %>% dplyr::select(date,'positiveIncrease')

# select start date 
DataUsing<- DataUsing %>% filter(date > 20200706)   
DataUsing<- as.data.frame(DataUsing)
DataUsing$date <- as.factor(DataUsing$date)
lenPos = length(DataUsing$positiveIncrease)
texasPos = rev(DataUsing$positiveIncrease)




Class = read.csv("~/Documents/datascience/DS7331/CoronaCurves-master/Corona_Curves_TX_NY_NJ.csv", header = TRUE)
#Class <- Class[rev(order(as.Date(Class$Date, format="%m/%d/%Y"))),]
#Class %>% select(State,Date,Daily_New_Cases)

nj = filter(Class, State == "NJ")
ny = filter(Class, State == "NY")
tx = filter(Class, State == "TX")


####### Future texas Data  ########
FutureData = read.csv("~/Documents/datascience/DS7331/CoronaCurves-master/Corona_MAE.csv", header = TRUE)
texasPos = FutureData$TX.New.Cases
lenTXFuture = length(texasPos)


### New Jersey ########
# VAR Model
onlineData = as.data.frame(texasPos *(8/29))
count = length(nj$Daily_New_Cases)
CoroVar = VAR(cbind(nj$Daily_New_Cases, nj$Pct_Change,nj$Pop_Pct,nj$Three_Day_Avg_Pct_Chg), type = "both", lag.max = 10)

preds_NJ = predict(CoroVar,n.ahead = lenPos)

ASE = mean((onlineData$`texasPos * (8/29)` - preds_NJ$fcst$y1[,1])^2)
ASE

MAE = sqrt(ASE)

VAR_NJ_DF = data.frame(State = "New Jersey", Model = 'VAR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

plot(preds_NJ)
  
#### New york ########
count = length(ny$Daily_New_Cases)
CoroVar = VAR(cbind(ny$Daily_New_Cases, ny$Pct_Change,ny$Pop_Pct,ny$Three_Day_Avg_Pct_Chg), type = "both", lag.max = 10)

preds_NY = predict(CoroVar,n.ahead = lenPos)

ASE = mean((onlineData$`texasPos * (8/29)` - preds_NY$fcst$y1[,1])^2)
ASE

MAE = sqrt(ASE)

VAR_NY_DF = data.frame(State = "New York",Model = 'VAR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

plot(preds_NY)

#### Ensemble

preds_TX = (preds_NY$fcst$y1[,1] + preds_NJ$fcst$y1[,1])/2

ASE = mean((onlineData$`texasPos * (8/29)` - preds_TX)^2)
MAE = sqrt(ASE)

VAR_TX_DF = data.frame(State = "Ensemble/Texas",Model = 'VAR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

df = rbind(VAR_NJ_DF,VAR_NY_DF,VAR_TX_DF)
df



