library(tseries)
library(forecast)
library(cluster)
### In case it is your first time in R, install them as: ###
#install.packages('tseries')
#install.packages('forecast')
#install.packages('cluster')

path <- 'C:\\Users\\~\\Parameters\\'

PassengersData <- na.omit(read.csv(paste0(path,'passengers.csv'),sep=','))$Passengers.thousands.

DataSet <- PassengersData

# Use the last year as test set (validation set) but remember that you need to redict for the following year, so you wont have the real data to compare. 

smp_size <- length(DataSet)-12
train_ind <- c(1:smp_size)
train <- DataSet[train_ind]
test <- DataSet[-train_ind]

plot(train,col='blue',pch=20,main = '',xlab = 'Time [months]',ylab='# thousands of passengers',xlim=c(1,length(train)),ylim = c(min(DataSet),max(DataSet)))
lines(train,col='blue',lwd=1)

#Data for number of passengers during the past 11 years

plot(train,col='blue',pch=20,main = '',xlab = 'Time [months]',ylab='# thousands of passengers',xlim=c(1,length(DataSet)),ylim = c(min(DataSet),max(DataSet)))
lines(train,col='blue',lwd=1)
points(c((length(train)+1):(length(DataSet))),test,pch=20,col='black')
lines(c((length(train)+1):(length(DataSet))),test,lwd=1,col='black')
abline(v=length(train),col='grey')
legend('topleft',legend=c('Training Set','Test Set'),col=c('blue','black'),lwd=c(2,2))

##############################
### TIME SERIES ANLYSIS ######
##############################

# Plot ACF and PACF to check stationarity and correlations:
acf(train,main='ACF for Data',length(train)/2)
pacf(train,main='PACF for Data',length(train)/2)

# Apply transformations if required:
#log(train)
#diff(train,1)
#diff(log(train),1)
#diff(log(train),S)

###########################################################################################
# choose model Parameters
# Redifine train set
train <- 
###########################################################################################
# Create ARIMA order=(p,d,q) or SARIMA order=(p,d,q)x(P,D,Q)S:
  p = # AR(p)
  d = # differentiate order
  q = # MA(q)
  # Seasonal 
  S = # Seasonality  
  P = # AR(P) Seasonal
  D = # differentiate order Seasonal
  Q = # MA(Q) Seasonal
###########################################################################################

TS_Model <- arima(train,order = c(p,d,q), seasonal = list(order = c(P,D,Q), period = S))
# Check residuals of the model:
hist(TS_Model$residuals,prob = T,breaks = 20,col='deepskyblue1',main='Histogram residuals')
curve(dnorm(x, mean(TS_Model$residuals), sd(TS_Model$residuals)), add=TRUE, col="red", lwd=2)

qqnorm(TS_Model$residuals,main='Q-Q plot residuals')
qqline(TS_Model$residuals)

plot(c(fitted(TS_Model)),c(TS_Model$residuals),pch=20,col='red',xlab = 'Fitted Values',ylab='Residuals',main='Residual vs Fitted residuals')
abline(h=0)

acf(TS_Model$residuals,length(train)/2,main='ACF for residuals')

##############################################
### LETS MAKE PREDICTIONS AND PLOT THEM ######
##############################################

# Predict n.ahead steps:
Predictions <- predict(TS_Model, n.ahead = length(test))
# Plot predictions:
plot(train,type = 'l',col='blue',lwd=2,main = 'Fill Data',xlab = 'Time',ylab='Data',xlim=c(1,length(DataSet)),ylim = c(0,max(DataSet)))
lines(c((length(train)+1):(length(DataSet))),test,lwd=2,col='gray60')
lines(c((length(train)+1):(length(DataSet))),(Predictions$pred),lwd=1,col='red')
lines(c((length(train)+1):(length(DataSet))),(Predictions$pred+Predictions$se),lwd=1,lty=2,col='red')
lines(c((length(train)+1):(length(DataSet))),(Predictions$pred-Predictions$se),lwd=1,lty=2,col='red')
abline(v=length(train),col='grey')
legend('topleft',legend=c('Past observations','Real Observations','Mean Forecast','95% Conf. Intervals Forecast'),col=c('blue','gray60','red','red'),lwd=c(2,2,1,1),lty=c(1,1,1,2))

###################################
### WE CREATE SCENARIOS NOW  ######
###################################

###########################################################################################
## Choose these values according to what is asked in the assingment 
###########################################################################################
Steps <-  # predict the next Horizon steps
Scen <-     # initial number of scenarios
RedScen <-   # Reduce number of scenarios
###########################################################################################

# Data structure to store the scenarios
Scenarios <- matrix(NA,nrow =Steps ,ncol =Scen )
# Loop over scenarios and simulate with the arima model a prediction for the next 24 hours
for(w in 1:Scen){
  Scenarios[,w]  <- simulate(TS_Model, nsim=Steps, future=TRUE, seed=w)
}
# Remember to retransform the scenarios in case you applied log transfrom
TransformScenarios <- exp(Scenarios)

# Plot them
plot(train,type = 'l',col='blue',lwd=1,main = '',xlab = 'Time [months]',ylab='Number of passengers (thousands)',xlim=c(1,length(DataSet)),ylim = c(min(DataSet),max(TransformScenarios)))
points(train,col='blue',pch=20)
for(w in 1:Scen){lines(c((length(train)+1):(length(DataSet))),TransformScenarios[,w],col='grey',lwd=0.5)}
lines(c((length(train)+1):(length(DataSet))),test,lwd=1,col='black')
points(c((length(train)+1):(length(DataSet))),test,pch=20,col='black')
abline(v=length(train),col='grey')
legend('topleft',legend=c('Past observations','Actual Observations','Scenarios'),col=c('blue','black','grey'),lwd=c(1,1,1),lty=c(1,1,1))

# Reduce scenarios using Partitions around medoids
Scenarios_pam <- pam(t(TransformScenarios),RedScen)
Scenarios_pam$medoids

# Plot reduced scenarios
plot(train,type = 'l',col='blue',lwd=1,main = '',xlab = 'Time [months]',ylab='Number of passengers (thousands)',xlim=c(1,length(DataSet)),ylim = c(min(DataSet),max(TransformScenarios)))
points(train,col='blue',pch=20)
for(w in 1:RedScen){lines(c((length(train)+1):(length(DataSet))),Scenarios_pam$medoids[w,],col='grey',lwd=0.5)}
lines(c((length(train)+1):(length(DataSet))),test,lwd=1,col='black')
points(c((length(train)+1):(length(DataSet))),test,pch=20,col='black')
abline(v=length(train),col='grey')
legend('topleft',legend=c('Past observations','Actual Observations','Scenarios'),col=c('blue','black','grey'),lwd=c(1,1,1),lty=c(1,1,1))


# Probability Distribution for reduced number of scenarios
Probability <- 1/length(1:Scen)
Reduce_prob <- c() 
for (i in 1:RedScen) {
  Reduce_prob[i] <- Probability *length( Scenarios_pam$clustering[ Scenarios_pam$clustering==i])
}

# Store these two values for later 
Passengers_Scenarios  <- round(t(Scenarios_pam$medoids)*1000,0) # I multipled by 1000 cause data in the Time Series model are given in thousands
Probability_Scenarios <- Reduce_prob

# Transfrom the values to CSV files that GAMS can directly digest
Passengers_CSV <- Passengers_Scenarios
colnames(Passengers_CSV) <- c(',s1',paste0('s',2:RedScen))
rownames(Passengers_CSV) <- paste0('t',1:Steps)

Probability_CSV <- matrix(Probability_Scenarios,RedScen,1)
rownames(Probability_CSV) <- paste0('s',1:RedScen)

write.table(Passengers_CSV,file =  paste0(path,'Scenarios_Passengers.csv'), sep=",",  dec=".",row.names=TRUE,col.names=TRUE,quote=FALSE)
write.table(Probability_CSV,file =  paste0(path,'Scenarios_Probability.csv'), sep=",",  dec=".",row.names=TRUE,col.names=FALSE,quote=FALSE)

