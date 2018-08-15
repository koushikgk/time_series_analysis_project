library(forecast)
library(tseries)
library(itsmr)
freq = 288

#---------READING THE DATA INTO DATAFRAME---------------------------------------
mydata = read.csv("data.csv", header = TRUE)
mydata$X <- NULL 
dates <- as.POSIXct(mydata[,1],format="%m/%d/%Y %H:%M:%S",tz=Sys.timezone())
dat <- ts(tsclean(mydata[,2]))
png('dat.png')
plot(dates,dat, type = 'l', xlab = "Date", ylab = "Wind speed in MPH", main = "5-minute average Wind Speeds at Tillamook");grid(lwd = 2)
dev.off()
png('acf,pacf.png')
plota(dat,h = 'full')
dev.off()
sink("adf.txt")
adf.test(dat, k = 288)
sink()

#---------SPLITTING DATA INTO TRAINING AND TESTING------------------------------
smp_size <- floor(6/7 * length(dat))
train_dat <- window(dat, end = smp_size)
test_dat <- window(dat, start = smp_size+1)
train_dates <- window(dates, end = smp_size)
test_dates <- window(dates, start = smp_size+1)
adf.test(train_dat, k = 288)

#----------MODEL FITTING USING TRAINING DATA-------------------------------------
l <- length(test_dates)
detach(package:itsmr)
train_dat <- ts(train_dat, frequency = freq)
fit2 <- HoltWinters(train_dat, start.periods = 2)
f2 <- forecast(fit2, h = l)
fit1 <- stlm(train_dat, s.window = 2, method = "arima", lambda = 'auto', biasadj = TRUE)
f1 <- forecast(fit1, h = l)
sink("fit1.txt")
fit1$model
sink()
sink("fit2.txt")
fit2
sink()

pred1 <- ts(f1$mean, frequency = 1, start = smp_size+1, end = smp_size+288)
pred_er1 <- test_dat-pred1
pred2 <- ts(f2$mean, frequency = 1, start = smp_size+1, end = smp_size+288)
pred_er2 <- test_dat-pred2
#----------------PLOTTING THE FIT AND FORECASTS------------------------------------
png('stlm.png')
par(mfrow=c(2,1))
plot(dates, c(train_dat,ts(NA, start = 1, end = 288)), lwd = 2, main = "Forecasts from STL + ARIMA", type = 'l', ylab = "Wind speed in MPH", xlab = "Dates")
grid(lwd = 2)
polygon(c(test_dates,rev(test_dates)),c(f1$lower[,"95%"],rev(f1$upper[,"95%"])),col="grey", border = NA)
lines(train_dates,fit1$fitted, col = 'red')
lines(test_dates,f1$mean, col = 'blue', lwd = 1)
lines(test_dates,test_dat, col = 'black', lwd =1)
legend("bottomleft", legend=c("Data", "Fitted", "Predicted"), col=c("black", "red", "blue"),lty = 1:1,lwd = 2:2, cex=0.8)
plot(test_dates,test_dat, type = "l", ylab = "Wind speed in MPH", xlab = "Time, 8 Jan 2013");grid(lwd = 2)
polygon(c(test_dates,rev(test_dates)),c(f1$lower[,"95%"],rev(f1$upper[,"95%"])),col="grey", border = NA)
lines(test_dates,pred1, col = "blue")
lines(test_dates,test_dat, col = "black")
legend("topleft", legend=c("Data", "Predicted"), col=c("black", "blue"),lty = 1:1,lwd = 2:2, cex=0.8)
dev.off()

png('hw.png')
par(mfrow=c(2,1))
plot(dates, ylim = c(-30,40) , c(train_dat,ts(NA, start = 1, end = 288)), lwd = 2, main = "Forecasts from HoltWinters", type = 'l', ylab = "Wind speed in MPH", xlab = "Dates")
grid(lwd = 2)
polygon(c(test_dates,rev(test_dates)),c(f2$lower[,"95%"],rev(f2$upper[,"95%"])),col="grey", border = NA)
lines(train_dates,f2$fitted, col = 'red')
lines(test_dates,f2$mean, col = 'blue', lwd = 1)
lines(test_dates,test_dat, col = 'black', lwd =1)
legend("bottomleft", legend=c("Data", "Fitted", "Predicted"), col=c("black", "red", "blue"),lty = 1:1,lwd = 2:2, cex=0.8)
plot(test_dates,test_dat, type = "l", ylab = "Wind speed in MPH", xlab = "Time, 8 Jan 2013", ylim = c(-30,40));grid(lwd = 2)
polygon(c(test_dates,rev(test_dates)),c(f2$lower[,"95%"],rev(f2$upper[,"95%"])),col="grey", border = NA)
lines(test_dates,pred2, col = "blue")
lines(test_dates,test_dat, col = "black")
legend("topleft", legend=c("Data", "Predicted"), col=c("black", "blue"),lty = 1:1,lwd = 2:2, cex=0.8)
dev.off()

#----------RESIDUE EVALUATION------------------------------------------------------
res1 <- ts(fit1$residuals, frequency = 1)
res2 <- ts(train_dat-fit2$fitted[,"xhat"], frequency = 1)
library(itsmr)

attach(mtcars)
png('res1.png')
par(mfrow=c(3,1))
spec.pgram(res1, k=kernel("daniell", 40),taper=0,log="no",ylim=c(0,0.1))
acf(res1, lag.max = 'full', ylim=range(-1,1))
pacf(res1, lag.max = 'full', ylim=range(-1,1))
dev.off()
sink("adfres1.txt")
adf.test(res1, k = 288)
sink()

png('res2.png')
par(mfrow=c(3,1))
spec.pgram(res2, k=kernel("daniell", 40),taper=0,log="no")
acf(res2, lag.max = 'full', ylim=range(-1,1))
pacf(res2, lag.max = 'full', ylim=range(-1,1))
dev.off()
sink("adfres2.txt")
adf.test(res2, k = 288)
sink()

#----------FORECAST EVALUATION--------------------------------

sink("pred1.txt")
cat("Mean of the difference between Prediciton and Test set: ", mean(pred_er1), "\n")
cat("Variance of the difference between Prediciton and Test set: ", var(pred_er1), "\n")
cat("Correlation between Prediciton and Test set: ", cor(pred1, test_dat), "\n")
sink()

sink("pred2.txt")
cat("Mean of the difference between Prediciton and Test set: ", mean(pred_er2), "\n")
cat("Variance of the difference between Prediciton and Test set: ", var(pred_er2), "\n")
cat("Correlation between Prediciton and Test set: ", cor(pred2, test_dat), "\n")
sink()

#-------------FORECAST FOR JAN 9 2013---------------------
smp_size <- floor(7/7 * length(dat))
train_dat <- window(dat, end = smp_size)
train_dates <- window(dates, end = smp_size)

test_dates <- seq.POSIXt(as.POSIXct(c("2013-01-09")), as.POSIXct(c("2013-01-10")), by = "5 min")
test_dates <- test_dates[1:288]
l <- length(test_dates)
detach(package:itsmr)
train_dat <- ts(train_dat, frequency = freq)
fit2 <- HoltWinters(train_dat, start.periods = 2)
f2 <- forecast(fit2, h = l)
fit1 <- stlm(train_dat, s.window = 2, method = "arima", lambda = 'auto', biasadj = TRUE)
f1 <- forecast(fit1, h = l)

pred1 <- ts(f1$mean, frequency = 1, start = smp_size+1, end = smp_size+288)
pred2 <- ts(f2$mean, frequency = 1, start = smp_size+1, end = smp_size+288)
#----------------PLOTTING THE FIT AND FORECASTS------------------------------------
png('stlm_pred.png')
par(mfrow=c(2,1))
plot(c(dates,test_dates), c(train_dat,ts(NA, start = 1, end = 288)), lwd = 2, main = "Forecasts from STL + ARIMA", type = 'l', ylab = "Wind speed in MPH", xlab = "Dates")
grid(lwd = 2)
polygon(c(test_dates,rev(test_dates)),c(f1$lower[,"95%"],rev(f1$upper[,"95%"])),col="grey", border = NA)
lines(train_dates,fit1$fitted, col = 'red')
lines(test_dates,f1$mean, col = 'blue', lwd = 1)
legend("bottomleft", legend=c("Data", "Fitted", "Forecast"), col=c("black", "red", "blue"),lty = 1:1,lwd = 2:2, cex=0.8)
plot(test_dates,pred1, col = "blue", type = "l", ylab = "Wind speed in MPH", xlab = "Time, 9 Jan 2013", ylim = c(1,19));grid(lwd = 2)
polygon(c(test_dates,rev(test_dates)),c(f1$lower[,"95%"],rev(f1$upper[,"95%"])),col="grey", border = NA)
lines(test_dates,pred1, col = "blue", type = "l", lwd = 2);grid(lwd = 2)
dev.off()

png('hw_pred.png')
par(mfrow=c(2,1))
plot(c(dates,test_dates), c(train_dat,ts(NA, start = 1, end = 288)), ylim =c(-30,40) , lwd = 2, main = "Forecasts from Holt Winters", type = 'l', ylab = "Wind speed in MPH", xlab = "Dates")
grid(lwd = 2)
polygon(c(test_dates,rev(test_dates)),c(f2$lower[,"95%"],rev(f2$upper[,"95%"])),col="grey", border = NA)
lines(train_dates,f2$fitted, col = 'red')
lines(test_dates,f2$mean, col = 'blue', lwd = 1)
legend("bottomleft", legend=c("Data", "Fitted", "Forecast"), col=c("black", "red", "blue"),lty = 1:1,lwd = 2:2, cex=0.8)
plot(test_dates,pred2, col = "blue", type = "l", ylab = "Wind speed in MPH", xlab = "Time, 9 Jan 2013", ylim = c(-30,40));grid(lwd = 2)
polygon(c(test_dates,rev(test_dates)),c(f2$lower[,"95%"],rev(f2$upper[,"95%"])),col="grey", border = NA)
lines(test_dates,pred2, col = "blue", type = "l", lwd = 2);grid(lwd = 2)
dev.off()
closeAllConnections()

#----------------------------------------------------------

