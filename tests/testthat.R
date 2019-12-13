# install.packages("testthat")
library(testthat)
# install.packages("forecast")
# library(forecast)
# install.packages("stats")
# library(stats)
# install.packages("future.apply")
# library(future.apply)
library(dhReg)
Data1 <- runif(runif(200,100,1000)) #To generate random number for example
Data_ts <- ts(Data1)

test_that("summary of Dynamic harmonic regression model", {
  M <- dhr(Data=Data_ts,XREG=NULL,Range=list(1:2,1),Frequency=c(24,168),Criteria="aicc")
})

M <- dhr(Data=Data_ts,XREG=NULL,Range=list(1:2,1),Frequency=c(24,168),Criteria="aicc")
test_that("function to get best value of K used in dhr function", {
  fourier_K(M)
})

test_that("forecasting the time series data using Dynamic Harmonic Regression", {
  Fcast <- fc(Frequency = c(24,168), XREG_test = NULL, h = 10, Fit = M, Data = Data_ts)
})
