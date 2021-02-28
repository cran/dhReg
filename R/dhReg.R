#' @title Dynamic Harmonic Regression
#'
#' @description Building model for time series data with multiple seasonality using Dynamic Harmonic Regression
#'
#' @param Data a time series data
#' @param Range Range of k in fourier series
#' @param XREG independent variable if any
#' @param Frequency seasonal frequency(can be multiple)
#' @param Criteria can be "aicc", "aic", "bic"
#' @param maxp maximum value of Auto regressive term in auto.arima
#' @param maxq maximum value of Moving average term in auto.arima
#' @param maxd maximum value of integrated term in auto.arima
#'
#' @return summary of Dynamic harmonic regression model
#'
#' @importFrom forecast auto.arima
#' @importFrom forecast fourier
#' @importFrom forecast msts
#' @importFrom future.apply future_sapply
#' @importFrom stats ts
#' @import future
#' @import testthat
#' @examples
#' \donttest{
#' Data1 <- runif(runif(200,100,1000)) #To generate random number for example
#' Data_ts <- ts(Data1)
#' M <- dhr(Data=Data_ts,XREG=NULL,Range=list(1:2,1),Frequency=c(24,168),Criteria="aicc")
#' }
#' @export dhr

dhr <- function(Data, Range, XREG = NULL, Frequency, Criteria = "aicc", maxp = 5, maxq = 5, maxd = 5){#, p1 = NULL, q1 = NULL, d1 = NULL){
  if (future::supportsMulticore()) future::plan(future::multiprocess)
  if (future::supportsMulticore() == FALSE) future::plan(future::multisession)
  if(length(Frequency) == 5){
    bb <- future_sapply(Range[1][[1]], function(i){
      future_sapply(Range[2][[1]], function(j){
        future_sapply(Range[3][[1]], function(k){
          future_sapply(Range[4][[1]], function(l){
            future_sapply(Range[5][[1]], function(m){
              if(is.null(XREG)){
                fit <- auto.arima(Data, max.p = maxp, max.q = maxq, max.d = maxd, #p = p1, q = q1, d = d1,
                                  xreg = fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3], Frequency[4], Frequency[5])), K = c(i, j, k, l, m)),
                                  seasonal = F, lambda = 0)
              }else{
                fit <- auto.arima(Data, max.p = maxp, max.q = maxq, max.d = maxd, #p = p1, q = q1, d = d1,
                                  xreg = cbind(fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3], Frequency[4], Frequency[5])), K = c(i, j, k, l, m)), XREG),
                                  seasonal = F, lambda = 0)
              }
              c(fit = list(fit))
            })
          })
        })
      })
    })
    ##because in bb values are stored in 1st column and all rows then 2nd col all rows
    ##so first column is fixed by considering it as i
    if(Criteria == "aicc"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1 & length(Range[[3]]) == 1 & length(Range[[4]]) == 1 & length(Range[[5]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) + length(Range[[3]]) + length(Range[[4]]) + length(Range[[5]]) <= 6)){
        cri1 <- bb[[1]]$aicc
        cri <- append(cri1, bb[[2]]$aicc)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$aicc
          })
        })}
    }else if(Criteria == "aic"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1 & length(Range[[3]]) == 1 & length(Range[[4]]) == 1 & length(Range[[5]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) + length(Range[[3]]) + length(Range[[4]]) + length(Range[[5]]) <= 6)){
        cri1 <- bb[[1]]$aic
        cri <- append(cri1, bb[[2]]$aic)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$aic
          })
        })}
    }else if(Criteria == "bic"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1 & length(Range[[3]]) == 1 & length(Range[[4]]) == 1 & length(Range[[5]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) + length(Range[[3]]) + length(Range[[4]]) + length(Range[[5]]) <= 6)){
        cri1 <- bb[[1]]$bic
        cri <- append(cri1, bb[[2]]$bic)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$bic
          })
        })}}
      idx <- which.min(cri)
      best_fit <- bb[idx][[1]]
  }else if(length(Frequency) == 4){
    bb <- future_sapply(Range[1][[1]], function(i){
      future_sapply(Range[2][[1]], function(j){
        future_sapply(Range[3][[1]], function(k){
          future_sapply(Range[4][[1]], function(l){
            if(is.null(XREG)){
              fit <- auto.arima(Data, max.p = maxp, max.q = maxq, max.d = maxd, #p = p1, q = q1, d = d1,
                                xreg = fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3], Frequency[4])), K = c(i, j, k, l)),
                                seasonal = F, lambda = 0)
            }else{
              fit <- auto.arima(Data, max.p = maxp, max.q = maxq, max.d = maxd, #p = p1, q = q1, d = d1,
                                xreg = cbind(fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3], Frequency[4])), K = c(i, j, k, l)), XREG),
                                seasonal = F, lambda = 0)
            }
            c(fit = list(fit))
          })
        })
      })
    })
    ##because in bb values are stored in 1st column and all rows then 2nd col all rows
    ##so first column is fixed by considering it as i
    if(Criteria == "aicc"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1 & length(Range[[3]]) == 1 & length(Range[[4]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) + length(Range[[3]]) + length(Range[[4]]) <= 5)){
        cri1 <- bb[[1]]$aicc
        cri <- append(cri1, bb[[2]]$aicc)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$aicc
          })
        })}
    }else if(Criteria == "aic"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1 & length(Range[[3]]) == 1 & length(Range[[4]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) + length(Range[[3]]) + length(Range[[4]]) <= 5)){
        cri1 <- bb[[1]]$aic
        cri <- append(cri1, bb[[2]]$aic)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$aic
          })
        })}
    }else if(Criteria == "bic"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1 & length(Range[[3]]) == 1 & length(Range[[4]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) + length(Range[[3]]) + length(Range[[4]]) <= 5)){
        cri1 <- bb[[1]]$bic
        cri <- append(cri1, bb[[2]]$bic)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$bic
          })
        })}}
      idx <- which.min(cri)
      best_fit <- bb[idx][[1]]
  }else if(length(Frequency) == 3){
    bb <- future_sapply(Range[1][[1]], function(i){
      future_sapply(Range[2][[1]], function(j){
        future_sapply(Range[3][[1]], function(k){
          if(is.null(XREG)){
            fit <- auto.arima(Data, max.p = maxp, max.q = maxq, max.d = maxd, #p = p1, q = q1, d = d1,
                              xreg = fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3])), K = c(i, j, k)),
                              seasonal = F, lambda = 0)
          }else{
            fit <- auto.arima(Data, max.p = maxp, max.q = maxq, max.d = maxd, #p = p1, q = q1, d = d1,
                              xreg = cbind(fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3])), K = c(i, j, k)), XREG),
                              seasonal = F, lambda = 0)
          }
          c(fit = list(fit))
        })
      })
    })
    ##because in bb values are stored in 1st column and all rows then 2nd col all rows
    ##so first column is fixed by considering it as i
    if(Criteria == "aicc"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1 & length(Range[[3]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) + length(Range[[3]]) <= 4)){
        cri1 <- bb[[1]]$aicc
        cri <- append(cri1, bb[[2]]$aicc)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$aicc
          })
        })}
    }else if(Criteria == "aic"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1 & length(Range[[3]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) + length(Range[[3]]) <= 4)){
        cri1 <- bb[[1]]$aic
        cri <- append(cri1, bb[[2]]$aic)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$aic
          })
        })}
    }else if(Criteria == "bic"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1 & length(Range[[3]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) + length(Range[[3]]) <= 4)){
        cri1 <- bb[[1]]$bic
        cri <- append(cri1, bb[[2]]$bic)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$bic
          })
        })}}
      idx <- which.min(cri)
      best_fit <- bb[idx][[1]]
  }else if(length(Frequency) == 2){
    bb <- future_sapply(Range[1][[1]], function(i){
      future_sapply(Range[2][[1]], function(j){
        if(is.null(XREG)){
          fit <- auto.arima(Data, max.p = maxp, max.q = maxq, max.d = maxd, #p = p1, q = q1, d = d1,
                            xreg = fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2])), K = c(i, j)),
                            seasonal = F, lambda = 0)
        }else{
          fit <- auto.arima(Data, max.p = maxp, max.q = maxq, max.d = maxd, #p = p1, q = q1, d = d1,
                            xreg = cbind(fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2])), K = c(i, j)), XREG),
                            seasonal = F, lambda = 0)
        }
        c(fit = list(fit))
      })
    })
    ##because in bb values are stored in 1st column and all rows then 2nd col all rows
    ##so first column is fixed by considering it as i
    if(Criteria == "aicc"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) <= 3)){
        cri1 <- bb[[1]]$aicc
        cri <- append(cri1, bb[[2]]$aicc)
      }else{
      cri <- future_sapply(1:ncol(bb), function(i){
        future_sapply(1:nrow(bb), function(j){
          bb[j,i][[1]]$aicc
        })
      })}
    }else if(Criteria == "aic"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) <= 3)){
        cri1 <- bb[[1]]$aic
        cri <- append(cri1, bb[[2]]$aic)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$aic
          })
        })}
    }else if(Criteria == "bic"){
      if((length(Range[[1]]) == 1 & length(Range[[2]]) == 1) | (length(Range[[1]]) + length(Range[[2]]) <= 3)){
        cri1 <- bb[[1]]$bic
        cri <- append(cri1, bb[[2]]$bic)
      }else{
        cri <- future_sapply(1:ncol(bb), function(i){
          future_sapply(1:nrow(bb), function(j){
            bb[j,i][[1]]$bic
          })
        })}}
      idx <- which.min(cri)
      best_fit <- bb[idx][[1]]
  }else if(length(Frequency) == 1){
    bb <- future_sapply(Range[1][[1]], function(i){
      if(is.null(XREG)){
        fit <- auto.arima(Data, max.p = maxp, max.q = maxq, max.d = maxd, #p = p1, q = q1, d = d1,
                          xreg = fourier(msts(Data, seasonal.periods = c(Frequency[1])), K = c(i)),
                          seasonal = F, lambda = 0)
      }else{
        fit <- auto.arima(Data, max.p = maxp, max.q = maxq, max.d = maxd, #p = p1, q = q1, d = d1,
                          xreg = cbind(fourier(msts(Data, seasonal.periods = c(Frequency[1])), K = c(i)), XREG),
                          seasonal = F, lambda = 0)
      }
      c(fit = list(fit))
    })
    ##because in bb values are stored in 1st column and all rows then 2nd col all rows
    ##so first column is fixed by considering it as i
    if(Criteria == "aicc"){
      cri <- future_sapply(Range[1][[1]], function(i){
        bb[i]$fit$aicc
      })
      idx <- which.min(cri)
      best_fit <- bb[idx]$fit
    }else if(Criteria == "aic"){
      cri <- future_sapply(Range[1][[1]], function(i){
        bb[i]$fit$aicc
      })
      idx <- which.min(cri)
      best_fit <- bb[idx]$fit
    }else if(Criteria == "bic"){
      cri <- future_sapply(Range[1][[1]], function(i){
        bb[i]$fit$aicc
      })
      idx <- which.min(cri)
      best_fit <- bb[idx]$fit
    }
  }
}
#' @title Fourier K
#'
#' @description function to get best value of K used in dhr function
#'
#' @param Fit Model built using dhr function
#'
#' @return optimal value of K used in dhr function
#'
#' @importFrom forecast auto.arima
#' @importFrom forecast fourier
#' @importFrom forecast msts
#' @importFrom future.apply future_sapply
#' @importFrom stats ts
#' @import future
#' @import testthat
#' @examples
#' \donttest{
#' Data1 <- runif(runif(200,100,1000))#To generate random number for example
#' Data_ts <- ts(Data1)
#' M <- dhr(Data=Data_ts,XREG=NULL,Range=list(1:2,1),Frequency=c(24,168),Criteria="aicc")
#' fourier_K(M)
#' }
#' @export fourier_K
fourier_K <- function(Fit){
  names(Fit$coef)
  coeff <- names(Fit$coef)[grepl(pattern = "^S[0-9]?[0-9]?[0-9]-", names(Fit$coef)) | grepl(pattern = "^C[0-9]?[0-9]?[0-9]-", names(Fit$coef))]
  q <- strsplit(coeff, "-")
  K <- c()
  for(i in 1:length(q)){
    K[i] <- q[[i]][2]
  }
  K <- table(K)/2
  return(K)
}

#' @title forecast using Dynamic Harmonic Regression
#'
#' @description forecasting the time series data using Dynamic Harmonic Regression
#'
#' @param Frequency seasonal frequency(can be multiple frequency)
#' @param XREG_test independent variable of test data, if any
#' @param h how much further to forecast
#' @param Fit Model fitted using dhr function
#' @param Data a time series data used while building a model
#'
#' @return forecasted values
#' @importFrom forecast auto.arima
#' @importFrom forecast fourier
#' @importFrom forecast msts
#' @importFrom future.apply future_sapply
#' @importFrom stats ts
#' @importFrom forecast forecast
#' @import future
#' @import testthat
#' @examples
#' \donttest{
#' Data1 <- runif(runif(200,100,1000))#To generate random number for example
#' Data_ts <- ts(Data1)
#' M <- dhr(Data=Data_ts,XREG=NULL,Range=list(1:2,1),Frequency=c(24,168),Criteria="aicc")
#' Fcast <- fc(Frequency = c(24,168), XREG_test = NULL, h = 10, Fit = M, Data = Data_ts)
#' plot(Fcast)
#' }
#' @export fc
fc <- function(Frequency, XREG_test = NULL, h, Fit, Data){
  coeff <- names(Fit$coef)[grepl(pattern = "^S[0-9]?[0-9]?[0-9]-", names(Fit$coef)) | grepl(pattern = "^C[0-9]?[0-9]?[0-9]-", names(Fit$coef))]
  q <- strsplit(coeff, "-")
  K <- c()
  for(i in 1:length(q)){
    K[i] <- q[[i]][2]
  }
  K <- table(K)/2
  if(length(Frequency) == 5){
    if(!is.null(XREG_test)){
      forecast(Fit,
               xreg=cbind(fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3], Frequency[4], Frequency[5])), K = c(K[as.character(Frequency[1])], K[as.character(Frequency[2])], K[as.character(Frequency[3])], K[as.character(Frequency[4])], K[as.character(Frequency[5])]), h = h), XREG_test)
               , h = h)
    }else{
      forecast(Fit,
               xreg = fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3], Frequency[4], Frequency[5])), K = c(K[as.character(Frequency[1])], K[as.character(Frequency[2])], K[as.character(Frequency[3])], K[as.character(Frequency[4])], K[as.character(Frequency[5])]), h = h)
      , h = h)
    }
  }else if(length(Frequency) == 4){
    if(!is.null(XREG_test)){
      forecast(Fit,
               xreg = cbind(fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3], Frequency[4])), K = c(K[as.character(Frequency[1])], K[as.character(Frequency[2])], K[as.character(Frequency[3])], K[as.character(Frequency[4])]), h = h), XREG_test)
               , h = h)
    }else{
      forecast(Fit,
               xreg = fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3], Frequency[4])), K = c(K[as.character(Frequency[1])], K[as.character(Frequency[2])], K[as.character(Frequency[3])], K[as.character(Frequency[4])]), h = h)
               , h = h)
    }
  }else if(length(Frequency) == 3){
    if(!is.null(XREG_test)){
      forecast(Fit,
               xreg=cbind(fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3])), K = c(K[as.character(Frequency[1])], K[as.character(Frequency[2])], K[as.character(Frequency[3])]), h = h), XREG_test)
               , h = h)
    }else{
      forecast(Fit,
               xreg = fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2], Frequency[3])), K = c(K[as.character(Frequency[1])], K[as.character(Frequency[2])], K[as.character(Frequency[3])]), h = h)
      , h = h)
    }
  }else if(length(Frequency) == 2){
    if(!is.null(XREG_test)){
      forecast(Fit,
               xreg=cbind(fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2])), K = c(K[as.character(Frequency[1])], K[as.character(Frequency[2])]), h = h), XREG_test)
               , h = h)
    }else{
      forecast(Fit,
               xreg = fourier(msts(Data, seasonal.periods = c(Frequency[1], Frequency[2])), K = c(K[as.character(Frequency[1])], K[as.character(Frequency[2])]), h = h)
      , h = h)
    }
  }else if(length(Frequency) == 1){
    if(!is.null(XREG_test)){
      forecast(Fit,
               xreg=cbind(fourier(msts(Data, seasonal.periods = c(Frequency[1])), K = c(K[1]), h = h), XREG_test)
               , h = h)
    }else{
      forecast(Fit,
               xreg = fourier(msts(Data, seasonal.periods = c(Frequency[1])), K = c(K[1]), h = h)
      , h = h)
    }
  }
}
