library(tswge)
library(orcutt)
library(vars)
library(nnfor)
library(forecast)
library(dplyr)

uri.data = 'https://covidtracking.com/api/v1/states/tx/daily.csv'
texas = read.csv(uri.data, header = T)

# National to get population data ############
#uri.data = 'https://covidtracking.com/api/v1/states/daily.csv'
#national = read.csv(uri.data, header = T)

#usCovidCleanNat <- national %>% replace(is.na(.), 0)

#DataUsingNat <- usCovidCleanNat %>% dplyr::select(date,positiveIncrease)

# We ignored using only data quality that scored at C because we have 1100 non labelled which is a lot to ignore. We only group by date
#dataGrouped <-DataUsingNat %>% dplyr::group_by(date) %>% dplyr::group_by(date) %>% summarise_each(tibble::lst(sum))
#ddply(usCovid, .(date), summarise, medABV=sum(positive))
#dataGrouped <- dataGrouped %>% filter(date > 20200706)

#####################################################
usCovidClean <- texas %>% replace(is.na(.), 0)
DataUsing <- usCovidClean %>% dplyr::select(date,'positiveIncrease')

# select start date 
DataUsing<- DataUsing %>% filter(date > 20200706)   
DataUsing<- as.data.frame(DataUsing)
DataUsing$date <- as.factor(DataUsing$date)
#DataUsing$Pop_PositiveIncrease <- dataGrouped$sum
#DataUsing$Pop_Pct <- (DataUsing$positiveIncrease / DataUsing$Pop_PositiveIncrease)

#DataUsing$lagpos = dplyr::lead(DataUsing$positiveIncrease,1)

#DataUsing$Pct_Change <- DataUsing$positiveIncrease - DataUsing$lagpos

####### This was used when I didn't know that there was a columnn named positive increased.
#DataUsing$lagpos = dplyr::lead(DataUsing$positive,1)

#DataUsing$Daily <- DataUsing$positive - DataUsing$lagpos

#dataUsed <- DataUsing %>% dplyr::select(date, Daily)

#dataUsed <-dataUsed %>% filter(!is.na(Daily))

#lenPos = length(dataUsed$Daily)
#texasPos = rev(dataUsed$Daily)

lenPos = length(DataUsing$positiveIncrease)
texasPos = rev(DataUsing$positiveIncrease)

############################## Class Data ###################################
Class = read.csv("~/Documents/datascience/DS7331/CoronaCurves-master/Corona_Curves_TX_NY_NJ.csv", header = TRUE)
#Class <- Class[rev(order(as.Date(Class$Date, format="%m/%d/%Y"))),]
#Class %>% select(State,Date,Daily_New_Cases)

nj = filter(Class, State == "NJ")
ny = filter(Class, State == "NY")
tx = filter(Class, State == "TX")

#NewClass <- rbind(ny, nj, tx)
Class <- rbind(nj, ny, tx)

########################################### Combine latest Texas data ##########################

onlineData = as.data.frame(texasPos)
classData <- Class$Daily_New_Cases

#names(classData)[1] <- "Daily_New_Cases"
names(onlineData)[1] <- "Daily_New_Cases"

foreAll = rbind(classData,onlineData)

#CoroNewDaily = rev(Coro$Daily_New_Cases)
#CoroNewDaily = Coro$Daily_New_Cases

plotts.sample.wge(NewCoro$Daily_New_Cases)
plotts.sample.wge(Class$Daily_New_Cases)


factor.wge(NewCoro$Daily_New_Cases)
factor.wge(Class$Daily_New_Cases)

acf(Class$Daily_New_Cases)
pacf(Class$Daily_New_Cases)

### ARMA Forcast #####
y.mle = est.arma.wge(Class$Daily_New_Cases, p = 1, q = 1)
foreARMA = fore.arma.wge(Class$Daily_New_Cases, phi=y.mle$phi, n.ahead = 10, lastn = T, limits = F)
foreARMA$f

######    Plots
plot(seq(1,count,1), Class$Daily_New_Cases[1:count], type = "l",xlim = c(0,count), ylab = "Texas Daily Corona Virus", main = "5 Days Forecast")
lines(seq(count-9,count,1), preds$fcst$y1[,1], type = "l", col = "red")

######### ARUMA Forecast ####
x = Class$Daily_New_Cases
y = artrans.wge(x, phi.tr = c(0,0,0,0,0,0,0,1))

pacf(y)

aic5.wge(y, p=0:10, q = 0:2, type = 'bic')
aic5.wge(Class$Daily_New_Cases, p=0:10, q = 0:2, type = 'bic')

y.mle = est.arma.wge(y, p = 1, q = 1)

sa = fore.aruma.wge(Class$Daily_New_Cases, phi=y.mle$phi,  s = 7, n.ahead = 10, lastn = T, limits = F)
sa


plotts.sample.wge(y.mle$res) # looks like white noise.
lj = ljung.wge(y.mle$res, p = 1)
lj$pval # 0.52 so we fail to reject Ho. There are not enough evidence to suggest that the residuals are serially correlated. we think it is a white noise residual

lj = ljung.wge(y.mle$res, p = 1, K = 48)
lj$pval # 0.52 so we fail to reject the null. There are not enough evidence to suggest that the residuals are serially correlated. This maybe white noise


####################### Predict Class explanatory variables #########################
## Pct_change
plotts.sample.wge(Class$Pct_Change)

aic = aic.wge(Class$Pct_Change,p = 0:10, q = 0:2, type = 'bic')
Pct_Change = est.arma.wge(Class$Pct_Change,p = aic$p, q = aic$q)
acf(Pct_Change$res)
ljung.wge(Pct_Change$res, p = aic$p, q = aic$q)
ljung.wge(Pct_Change$res, p = aic$p, q = aic$q, K = 48)

Pct_ChangeFore = fore.arma.wge(Class$Pct_Change, phi = Pct_Change$phi, theta = Pct_Change$theta, n.ahead = lenTXFuture)

Pct_ChangeFore = Pct_ChangeFore$f
######################################################
# Curve Day
plotts.sample.wge(Class$Curve_Day)
aic = aic.wge(Class$Curve_Day,p = 0:10, q = 0:2, type = 'bic')
Curve_Day = est.arma.wge(Class$Curve_Day,p = aic$p, q = aic$q)
acf(Curve_Day$res)
ljung.wge(Curve_Day$res, p = aic$p, q = aic$q)
ljung.wge(Curve_Day$res, p = aic$p, q = aic$q, K = 48)

Curve_DayFore = fore.arma.wge(Class$Curve_Day, phi = Curve_Day$phi, theta = Curve_Day$theta, n.ahead = lenTXFuture)

Curve_DayFore = Curve_DayFore$f

###########################
## Pop_ct
plotts.sample.wge(Class$Pop_Pct)
aic = aic.wge(Class$Pop_Pct,p = 0:10, q = 0:2, type = 'bic')
Pop_Pct = est.arma.wge(Class$Pop_Pct,p = aic$p, q = aic$q)
acf(Pop_Pct$res)
ljung.wge(Pop_Pct$res, p = aic$p, q = aic$q)
ljung.wge(Pop_Pct$res, p = aic$p, q = aic$q, K = 48)

Pop_PctFore = fore.arma.wge(Class$Pop_Pct, phi = Pop_Pct$phi, theta = Pop_Pct$theta, n.ahead = lenTXFuture)

Pop_PctFore = Pop_PctFore$f

##################################


# VAR Model Initial
count = length(Class$Daily_New_Cases)
#CoroVar = VAR(cbind(Class$Daily_New_Cases,Class$Three_Day_Avg_Pct_Chg,Class$Curve_Day, Class$Pop_Pct, Class$Pct_Change), type = "both", lag.max = 10)
CoroVar = VAR(cbind(Class$Daily_New_Cases, Class$Pct_Change,Class$Curve_Day,Class$Pop_Pct), type = "both", lag.max = 10)

preds = predict(CoroVar,n.ahead = 10, lastn = T) 

ASE = mean((Class$Daily_New_Cases[count-9:count] - preds$fcst$y1[,1])^2)
ASE

plot(preds)
dev.off()
plot(seq(1,count,1), Class$Daily_New_Cases[1:count], type = "l",xlim = c(0,count), ylab = "Texas Daily Corona Virus", main = "5 Days Forecast")
lines(seq(count-9,count,1), preds$fcst$y1[,1], type = "l", col = "red")

##################################
#### Final VAR model
# VAR Model
count = length(Class$Daily_New_Cases)
#CoroVar = VAR(cbind(Class$Daily_New_Cases,Class$Three_Day_Avg_Pct_Chg,Class$Curve_Day, Class$Pop_Pct, Class$Pct_Change), type = "both", lag.max = 10)
CoroVar = VAR(cbind(Class$Daily_New_Cases, Class$Pct_Change,Class$Curve_Day), type = "both", lag.max = 10)

preds = predict(CoroVar,n.ahead = lenPos)

ASE = mean((texasPos - preds$fcst$y1[,1])^2)
ASE

MAE = sqrt(ASE)

VAR_DF = data.frame(Model = 'VAR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

plot(preds)
dev.off()
plot(seq(1,count+lenPos,1), foreAll$Daily_New_Cases, type = "l",xlim = c(0,count+lenPos), ylab = "Texas Daily Corona Virus Forecast", main = paste(as.character(lenPos),"Days Forecast"))
lines(seq(count+ 1,count+lenPos,1), preds$fcst$y1[,1], type = "l", col = "red")

#################################################################################

##############   MLP    ####################
CoroDF = data.frame(popct = Class$Pop_Pct,pctchange =Class$Pct_Change, cv = Class$Curve_Day)

foreDF = data.frame(cv = as.double( Curve_DayFore), popct = as.double(Pop_PctFore),pctchange = Pct_ChangeFore)

foreAllDF = rbind(CoroDF,foreDF)

CoroDF = data.frame(popct = ts(Class$Pop_Pct),pctchange = ts(Class$Pct_Change), cv = ts(Class$Curve_Day))

fit.mlp = mlp(ts(Class$Daily_New_Cases),reps = 20,comb = "mean",xreg = CoroDF)
fit.mlp
plot(fit.mlp)
fore.mlp = forecast(fit.mlp, h = lenPos,xreg = foreAllDF)
plot(fore.mlp) 
ASE = mean((texasPos - fore.mlp$mean)^2)
ASE

MAE = sqrt(ASE)

MLP_DF = data.frame(Model = 'MLP', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

################################ MLR ########################################
ccf(Class$Pop_Pct,Class$Daily_New_Cases)
ccf(Class$Pct_Change,Class$Daily_New_Cases)
ccf(Class$Three_Day_Avg_Pct_Chg,Class$Daily_New_Cases)
ccf(Class$Curve_Day,Class$Daily_New_Cases)
Class$lagCurve = dplyr::lag(Class$Curve_Day,3)
#Class$lagPct_Change_M1 = dplyr::lag(Class$lagPct_Change_M1,4)

#Class$lagPct_Change_M2 = dplyr::lag(Class$lagPct_Change_M2,6)
#Class$lagPct_Change_M3 = dplyr::lag(Class$lagPct_Change_M3,11)

ccf(Class$Pct_Change_M1 ,Class$Daily_New_Cases)
ccf(Class$Pct_Change_M2,Class$Daily_New_Cases)
ccf(Class$Pct_Change_M3,Class$Daily_New_Cases)
ccf(Class$Pct_Change_M4,Class$Daily_New_Cases)
ccf(Class$Pct_Change_M5,Class$Daily_New_Cases)
ccf(Class$Pct_Change_M6,Class$Daily_New_Cases)
ccf(Class$Pct_Change_M7,Class$Daily_New_Cases)
ccf(Class$Three_Day_Avg_Pct_Chg,Class$Daily_New_Cases)
ccf(Class$Seven_Day_Avg_Pct_Chg,Class$Daily_New_Cases)
ccf(Class$Pop_Pct_M1 ,Class$Daily_New_Cases)
ccf(Class$Pop_Pct_M2 ,Class$Daily_New_Cases)
ccf(Class$Pop_Pct_M3 ,Class$Daily_New_Cases)
ccf(Class$Pop_Pct_M4 ,Class$Daily_New_Cases)
ccf(Class$Pop_Pct_M5 ,Class$Daily_New_Cases)
ccf(Class$Pop_Pct_M6 ,Class$Daily_New_Cases)
ccf(Class$Pop_Pct_M7 ,Class$Daily_New_Cases)
ccf(Class$Three_Day_Avg_Pop_Pct,Class$Daily_New_Cases)
ccf(Class$Seven_Day_Avg_Pop_Pct,Class$Daily_New_Cases)

######################### Begin ##########################
Corofit1 = lm(Daily_New_Cases~as.double( Class$Pop_Pct) + as.double( Class$Pct_Change) + as.double( Class$lagCurve), data = Class)
phi = aic.wge(Corofit1$residuals)

acf(Corofit1$residuals)
ljung.wge(Corofit1$residuals) # pval = .066
ljung.wge(Corofit1$residuals, K = 48) # pval = .0058

resids = fore.arma.wge(Corofit1$residuals,phi = phi$phi,n.ahead = lenPos)

foreDF = data.frame(lagCurve = as.double( Curve_DayFore), Pop_Pct = as.double(Pop_PctFore),Pct_Change = Pct_ChangeFore)
#predict trend manually
preds = predict(Corofit1, newxreg  = foreDF)

predsFinal = preds + resids$f

predsFinal <- as.numeric(predsFinal)

ASE = mean((texasPos - predsFinal)^2)
ASE

MAE = sqrt(ASE)

MLR_DF = data.frame(Model = 'MLR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

df = rbind(VAR_DF,MLP_DF,MLR_DF)
df