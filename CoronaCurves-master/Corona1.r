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
#foreAll <- rbind(NJ, NY, TX,TX_Future)
#lenForeAll <- length(foreAll$date)
#lenTXFuture <- length(TX_Future$date)

colSums(is.na(Class))

####### Future texas Data  ########
FutureData = read.csv("~/Documents/datascience/DS7331/CoronaCurves-master/Corona_MAE.csv", header = TRUE)
texasPos = FutureData$TX.New.Cases
lenTXFuture = length(texasPos)
foreAll <- rbind(NJ, NY, TX,texasPos)


########################################### Combine latest Texas data ##########################
plotts.sample.wge(Class$positiveIncrease)

acf(Class$positiveIncrease)
pacf(Class$positiveIncrease)

################################## forecast explanatory Variable. ##############
Class <- Class%>% mutate(positive_perc = (positiveIncrease/totalTestResults)*100)

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

######################## Death Increase Daily ###############################################
plotts.sample.wge(Class$deathIncrease)
aic = aic.wge(Class$deathIncrease,p = 0:10, q = 0:2, type = 'bic')
death = est.arma.wge(Class$deathIncrease,p = aic$p, q = aic$q)

acf(death$res)
ljung.wge(death$res, p = aic$p, q = aic$q)
ljung.wge(death$res, p = aic$p, q = aic$q, K = 48)

deathFore = fore.aruma.wge(Class$deathIncrease, phi = perc$phi, d = 1, theta = perc$theta, n.ahead = lenTXFuture)

deathFore = deathFore$f

###################### Negative Daily ###############################
negativeInc_01= artrans.wge(Class$negativeIncrease, c(rep(0,0),1))
aic = aic.wge(negativeInc_01,p = 0:10, q = 0:2, type = 'bic')
negative = est.arma.wge(negativeInc_01,p = aic$p, q = aic$q)

acf(negative$res)
ljung.wge(negative$res, p = aic$p, q = aic$q)
ljung.wge(negative$res, p = aic$p, q = aic$q, K = 48)

negativeFore = fore.aruma.wge(dataGrouped$negativeIncrease, phi = perc$phi, d = 1, theta = perc$theta, n.ahead = lenTXFuture)
negativeFore = negativeFore$f
######################### Total Test Forecast #########################
library(orcutt)
x = Class$totalTestResults
n = length(x)
t = 1:n
d = lm(x~t)
x.z = x - d$coefficients[1] - d$coefficients[2]*t


ar.z = aic.wge(x.z, p = 0:15, q = 0:2) # ar.z$p is 10

# Tranform the dependent varible daily rates 
y.trans = artrans.wge(x, phi.tr = ar.z$phi)

# Transform the independent variable t (trend)
t.trans = artrans.wge(t, phi.tr = ar.z$phi)
plotts.sample.wge(y.trans)    # looks like noise
plotts.sample.wge(t.trans)


# regress y hat t on T hat t using OLS
fitTotal = lm(y.trans~t.trans)

# evaluate the residuals( after Cochrane-Orcutt)

plotts.wge(fitTotal$residuals)
acf(fitTotal$residuals)
ljung.wge(fitTotal$residuals) # There is evidence that this is white noise. Pval = 0.99. Fail to reject.
ljung.wge(fitTotal$residuals, K = 48)

# phis and white noise variance
x.z.est = est.arma.wge(x.z, p = ar.z$p, q = ar.z$q) # sig plus noise phis

x.z.est$avar  # the variance

# ASE for signal plus noise
totalFore = fore.sigplusnoise.wge(dataGrouped$totalTestResults, max.p = 10, n.ahead = lenTXFuture)
totalFore = totalFore$f

##################### Actual Curve Day ###########################
# Actual Curve day should be used since we know it would be continuation of the last curve day from previous data
Curve_DayFore = FutureData$Curve_Day

################### Positive percent day #############################################################
perc_07= artrans.wge(dataGrouped$positive_perc, c(rep(0,6),1))
aic = aic.wge(perc_07,p = 0:10, q = 0:2, type = 'bic')
#aic = aic.wge(dataGrouped$positive_perc,p = 0:10, q = 0:2, type = 'bic')
#perc = est.arma.wge(dataGrouped$positive_perc,p = aic$p, q = aic$q)
perc = est.arma.wge(perc_07,p = aic$p, q = aic$q)

acf(perc$res) 
ljung.wge(perc$res, p = aic$p, q = aic$q)
ljung.wge(perc$res, p = aic$p, q = aic$q, K = 48)

percFore = fore.aruma.wge(dataGrouped$positive_perc, phi = perc$phi, s = 7, theta = perc$theta, n.ahead = lenTXFuture)

percFore = percFore$f

# VAR Model 
count = length(texasPos)
CoroVar = VAR(cbind(Class$positiveIncrease, Class$negativeIncrease, Class$Curve_Day), type = "both", lag.max = 10)

preds = predict(CoroVar,n.ahead = lenTXFuture) 

ASE = mean((texasPos - preds$fcst$y1[,1])^2)
ASE

MAE = sqrt(ASE)

VAR_DF = data.frame(Model = 'VAR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

plot(preds)


#################################################################################
# MLP Model
ClassDF = data.frame(negInc = ts(Class$negativeIncrease),cv = ts(Class$Curve_Day), test =  ts(Class$totalTestResults),hospCur = ts(Class$hospitalizedCurrently),dInc = ts(Class$deathIncrease))
fit.mlp = mlp(ts(Class$positiveIncrease),reps = 20,comb = "mode",xreg = ClassDF)
fit.mlp
plot(fit.mlp)
foreAllDF = data.frame(negInc = ts(foreAll$negativeIncrease),cv = ts(foreAll$Curve_Day), test =  ts(foreAll$totalTestResults),hospCur = ts(foreAll$hospitalizedCurrently),dInc = ts(foreAll$deathIncrease))


CoroDF = data.frame(negInc = Class$negativeIncrease,cv = Class$Curve_Day, test = Class$totalTestResults,hospCur = Class$hospitalizedCurrently,dInc = Class$deathIncrease)

foreDF = data.frame(cv = Curve_DayFore, negInc = negativeFore, test = totalFore,hospCur = hospitalFore,dInc = deathFore)

foreAllDF = rbind(CoroDF,foreDF)

fore.mlp = forecast(fit.mlp, h = lenTXFuture, xreg = foreAllDF)
plot(fore.mlp) 
ASE = mean((texasPos - fore.mlp$mean)^2)
ASE

MAE = sqrt(ASE)

MLP_DF = data.frame(Model = 'MLP', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)

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

Corofit = lm(positiveIncrease~negativeIncrease + leadCurve + leadhosp + lagtest, data = Class)
phi = aic.wge(Corofit$residuals)
ljung.wge(Corofit$residuals) # pval = .066
ljung.wge(Corofit$residuals, K = 48) # pval = .0058

#TX_FutureDF = TX_Future %>% dplyr::select(negativeIncrease,Curve_Day,totalTestResults,hospitalizedCurrently)
TX_FutureDF = data.frame(leadCurve = TX_Future$Curve_Day,lagtest= TX_Future$totalTestResults,leadhosp = TX_Future$hospitalizedCurrently, negativeIncrease = TX_Future$negativeIncrease)
foreDF = data.frame(negativeIncrease = as.double( negativeFore), lagtest = as.double(totalFore),leadhosp = hospitalFore,leadPositive_perc = percFore)


resids = fore.arma.wge(Corofit$residuals,phi = phi$phi,n.ahead = lenTXFuture)
#predict trend manually
preds = predict(Corofit, newdata = TX_FutureDF)


predsFinal = preds + resids$f

predsFinal <- as.numeric(predsFinal)

ASE = mean((texasPos - predsFinal)^2)
ASE

MAE = sqrt(ASE)

MLR_DF = data.frame(Model = 'MLR', Mean_ABS_Error = MAE, Mean_Square_Error = ASE)


df = rbind(VAR_DF,MLP_DF,MLR_DF)
df
