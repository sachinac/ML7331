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

#####################################################
Class <- rbind(NJ, NY, TX)
foreAll <- rbind(NJ, NY, TX,TX_Future)
lenForeAll <- length(foreAll$date)
lenTXFuture <- length(TX_Future$date)

colSums(is.na(Class))

########################################### Combine latest Texas data ##########################
plotts.sample.wge(Class$positiveIncrease)

acf(Class$positiveIncrease)
pacf(Class$positiveIncrease)

################################## forecast explanatory Variable. ##############
Class <- Class%>% mutate(positive_perc = (positiveIncrease/totalTestResults)*100)

ggpairs(Class %>% dplyr::select(positiveIncrease,negativeIncrease, deathIncrease,totalTestResults,positive_perc,hospitalizedCurrently,Curve_Day))

##################################################
################## Hospitalized Currrently ######################################################
plotts.sample.wge(Class$hospitalizedCurrently)

aic = aic.wge(Class$hospitalizedCurrently,p = 0:10, q = 0:2, type = 'bic')
hospital = est.arma.wge(Class$hospitalizedCurrently,p = aic$p, q = aic$q)
acf(hospital$res)
ljung.wge(hospital$res, p = aic$p, q = aic$q)
ljung.wge(hospital$res, p = aic$p, q = aic$q, K = 48)

hospitalFore = fore.arma.wge(Class$hospitalizedCurrently, phi = hospital$phi, theta = hospital$theta, n.ahead = lenTXFuture)

hospitalFore = hospitalFore$f

######################## Negative Increase Daily ###############################################
plotts.sample.wge(Class$deathIncrease)
aic = aic.wge(Class$deathIncrease,p = 0:10, q = 0:2, type = 'bic')
death = est.arma.wge(Class$deathIncrease,p = aic$p, q = aic$q)

acf(death$res)
ljung.wge(death$res, p = aic$p, q = aic$q)
ljung.wge(death$res, p = aic$p, q = aic$q, K = 48)

deathFore = fore.aruma.wge(Class$deathIncrease, phi = perc$phi, d = 1, theta = perc$theta, n.ahead = lenTXFuture)

deathFore = deathFore$f

############################# Curve Day ###########################################################
plotts.sample.wge(Class$Curve_Day)
aic = aic.wge(Class$Curve_Day,p = 0:10, q = 0:2, type = 'bic')
death = est.arma.wge(Class$Curve_Day,p = aic$p, q = aic$q)

acf(death$res)
ljung.wge(death$res, p = aic$p, q = aic$q)
ljung.wge(death$res, p = aic$p, q = aic$q, K = 48)

deathFore = fore.aruma.wge(Class$Curve_Day, phi = perc$phi, d = 1, theta = perc$theta, n.ahead = lenTXFuture)

deathFore = deathFore$f

###################################################################################################################

# VAR Model 
count = length(Class$positiveIncrease)
CoroVar = VAR(cbind(Class$positiveIncrease, Class$negativeIncrease, Class$Curve_Day), type = "both", lag.max = 10)

preds = predict(CoroVar,n.ahead = lenTXFuture) 

ASE = mean((TX_Future$positiveIncrease - preds$fcst$y1[,1])^2)
ASE

MAE = sqrt(ASE)

VAR_DF = data.frame(Model = 'VAR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

plot(preds)

dev.off()
plot(seq(1,lenForeAll,1), foreAll$positiveIncrease, type = "l",xlim = c(0,lenForeAll), ylab = "Cardiac Mortality", main = "20 Week Cardiac Mortality Forecast")
lines(seq(count+1,lenForeAll,1), preds$fcst$y1[,1], type = "l", col = "red")

#################################################################################
# MLP Model
ClassDF = data.frame(negInc = ts(Class$negativeIncrease),cv = ts(Class$Curve_Day), test =  ts(Class$totalTestResults),hospCur = ts(Class$hospitalizedCurrently),dInc = ts(Class$deathIncrease))
fit.mlp = mlp(ts(Class$positiveIncrease),reps = 20,comb = "mode",xreg = ClassDF)
fit.mlp
plot(fit.mlp)
foreAllDF = data.frame(negInc = ts(foreAll$negativeIncrease),cv = ts(foreAll$Curve_Day), test =  ts(foreAll$totalTestResults),hospCur = ts(foreAll$hospitalizedCurrently),dInc = ts(foreAll$deathIncrease))
fore.mlp = forecast(fit.mlp, h = lenTXFuture, xreg = foreAllDF)
plot(fore.mlp) 
ASE = mean((TX_Future$positiveIncrease - fore.mlp$mean)^2)
ASE

MAE = sqrt(ASE)

MLP_DF = data.frame(Model = 'MLP', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

dev.off()
plot(seq(1,lenForeAll,1), foreAll$positiveIncrease, type = "l",xlim = c(0,lenForeAll), ylab = "Cardiac Mortality", main = "20 Week Cardiac Mortality Forecast")
lines(seq(count+1,lenForeAll,1), fore.mlp$mean, type = "l", col = "red")

################################ MLR ########################################
ccf(Class$negativeIncrease,Class$positiveIncrease)
ccf(Class$Curve_Day,Class$positiveIncrease)
ccf(Class$totalTestResults,Class$positiveIncrease)
ccf(Class$hospitalizedCurrently,Class$positiveIncrease)
ccf(Class$deathIncrease,Class$positiveIncrease)
Class$leadCurve = dplyr::lead(Class$Curve_Day,5)
Class$lagtest = dplyr::lag(Class$totalTestResults,3)
Class$leadhosp = dplyr::lead(Class$hospitalizedCurrently,2)
Class$leaddeath = dplyr::lead(Class$deathIncrease,4)

plotts.parzen.wge(Class$negativeIncrease)
Bpart1High <- butterworth.wge(Class$positiveIncrease,order = 1, type = 'high', cutoff = 0.03)
plotts.sample.wge(Bpart1High$x.filt)


######################### Begin ##########################

Corofit = lm(positiveIncrease~negativeIncrease + leadCurve + leadhosp + lagtest + leaddeath, data = Class)
phi = aic.wge(Corofit$residuals)
ljung.wge(Corofit$residuals) # pval = .066
ljung.wge(Corofit$residuals, K = 48) # pval = .0058

#TX_FutureDF = TX_Future %>% dplyr::select(negativeIncrease,Curve_Day,totalTestResults,hospitalizedCurrently)
TX_FutureDF = data.frame(leadCurve = TX_Future$Curve_Day,lagtest= TX_Future$totalTestResults,leadhosp = TX_Future$hospitalizedCurrently, negativeIncrease = TX_Future$negativeIncrease, leaddeath = TX_Future$deathIncrease)


resids = fore.arma.wge(Corofit$residuals,phi = phi$phi,n.ahead = lenTXFuture)
#predict trend manually
#preds = predict(Corofit, newdata = TX_FutureDF)
preds = predict(Corofit)

predsFinal = preds + resids$f

predsFinal <- as.numeric(predsFinal)

ASE = mean((TX_Future$positiveIncrease - predsFinal)^2)
ASE

MAE = sqrt(ASE)

MLR_DF = data.frame(Model = 'MLR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)


df = rbind(VAR_DF,MLP_DF,MLR_DF)
df
