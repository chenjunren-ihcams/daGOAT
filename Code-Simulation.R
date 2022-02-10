rm(list = ls())


# Load libraries - when fail to import, please run:
# install.packages('library-name')
library(pROC)
library(PRROC)
library(Matrix)
library(xgboost)



# Define functions
mutual.information <- function(x, z) {
  p00 <- sum(x == 0 & z == 0)/length(z)
  p01 <- sum(x == 0 & z == 1)/length(z)
  p10 <- sum(x == 1 & z == 0)/length(z)
  p11 <- sum(x == 1 & z == 1)/length(z)
  p1. <- sum(x == 1)/length(z)
  p0. <- sum(x == 0)/length(z)
  p.1 <- sum(z == 1)/length(z)
  p.0 <- sum(z == 0)/length(z)
  out <- 0
  if (p00 > 0)
    out <- out + p00 * log(p00/p0./p.0)
  if (p01 > 0)
    out <- out + p01 * log(p01/p0./p.1)
  if (p10 > 0)
    out <- out + p10 * log(p10/p1./p.0)
  if (p11 > 0)
    out <- out + p11 * log(p11/p1./p.1)
  out
}


composite.xvalid <- function(keys, train, test, start.t = start.t, end.t = end.t,
                             engraftment.aware = engraftment.aware, use.smoothing = TRUE, use.laplace = TRUE,
                             use.stationary = TRUE) {
  
  if (FALSE) {
    "
        yy & y would be seen in pair in this fuction.
        yy stands for negative case 
        y stands for positive case
        "
  }
  
  
  ## Data manipulation
  all.patients <- complete.patients
  all.patient.data <- complete.patient.data
  all.severe <- complete.severe
  
  patients <- all.patients[train]
  severe <- intersect(all.severe, all.patients[train])
  
  
  ## Build Model - thresholds
  for (i in 1:length(train)) {
    m <- train[i]
    if (df$Severe.aGVHD.within.100.days[m] == 1) {
      onset.day <- df$Severe.aGVHD.onset.day[m]
      all.patient.data[[m]][1:min(onset.day - pre.window.size - 1,
                                  end.t) + offset.t, ] <- NA
    }
  }
  
  
  thresholds <- array(NA, dim = c(day, length(keys)))
  for (k in 1:length(keys)) {
    key <- keys[k]
    
    yy <- c()
    for (i in 1:length(patients)) {
      patient <- patients[i]
      if (!(patient %in% severe)) {
        m <- (1:length(all.patients))[all.patients == patient]
        patient.data <- all.patient.data[[m]]
        j <- (1:ncol(patient.data))[names(patient.data) == key]
        yy <- cbind(yy, patient.data[, j])
      }
    }
    
    y <- c()
    for (i in 1:length(severe)) {
      patient <- severe[i]
      m <- (1:length(all.patients))[all.patients == patient]
      patient.data <- all.patient.data[[m]]
      j <- (1:ncol(patient.data))[names(patient.data) == key]
      y <- cbind(y, patient.data[, j])
    }
    
    
    for (t in start.t:end.t + offset.t) {
      if (sum(!is.na(y[t, ])) > 1 & sum(!is.na(yy[t, ])) > 5) {
        min.val <- min(y[t, ], yy[t, ], na.rm = T)
        max.val <- max(y[t, ], yy[t, ], na.rm = T)
        
        vals <- seq(from = min.val, to = max.val, length.out = 100)
        shannon <- array(NA, length(vals))
        for (iiii in 1:length(vals)) {
          thres <- vals[iiii]
          x <- as.numeric(c(y[t, ], yy[t, ]) > thres)
          z <- c(array(1, ncol(y)), array(0, ncol(yy)))
          I <- !is.na(x)
          x <- x[I]
          z <- z[I]
          shannon[iiii] <- mutual.information(x, z)
        }
        
        best.iiii <- (1:length(shannon))[shannon == max(shannon)][1]
        
        thres <- vals[best.iiii]
        p.y.over <- sum(y[t, ] > thres, na.rm = T)/sum(!is.na(y[t,
        ]))
        p.yy.over <- sum(yy[t, ] > thres, na.rm = T)/sum(!is.na(yy[t,
        ]))
        
        thresholds[t, k] <- thres
      }
    }
  }
  
  
  thresholds.high <- array(NA, dim = c(day, length(keys)))
  thresholds.low <- array(NA, dim = c(day, length(keys)))
  for (k in 1:length(keys)) {
    thresholds.high[, k] <- quantile(thresholds[, k], 3/4, na.rm = T)
    thresholds.low[, k] <- quantile(thresholds[, k], 1/4, na.rm = T)
  }
  
  
  # Build Model - odds ratio
  p.y.high <- array(NA, dim = c(day, length(keys)))
  p.y.low <- array(NA, dim = c(day, length(keys)))
  p.yy.high <- array(NA, dim = c(day, length(keys)))
  p.yy.low <- array(NA, dim = c(day, length(keys)))
  for (k in 1:length(keys)) {
    key <- keys[k]
    
    yy <- c()
    for (i in 1:length(patients)) {
      patient <- patients[i]
      if (!(patient %in% severe)) {
        m <- (1:length(all.patients))[all.patients == patient]
        patient.data <- all.patient.data[[m]]
        j <- (1:ncol(patient.data))[names(patient.data) == key]
        yy <- cbind(yy, patient.data[, j])
      }
    }
    
    y <- c()
    for (i in 1:length(severe)) {
      patient <- severe[i]
      m <- (1:length(all.patients))[all.patients == patient]
      patient.data <- all.patient.data[[m]]
      j <- (1:ncol(patient.data))[names(patient.data) == key]
      y <- cbind(y, patient.data[, j])
    }
    
    for (t in start.t:end.t + offset.t) {
      if (!is.na(thresholds.high[t, k]) & sum(!is.na(y[t, ])) > 1 &
          sum(!is.na(yy[t, ])) > 5) {
        thres <- thresholds.high[t, k]
        p.y.high[t, k] <- sum(y[t, ] > thres, na.rm = T)/sum(!is.na(y[t,
        ]))
        p.yy.high[t, k] <- sum(yy[t, ] > thres, na.rm = T)/sum(!is.na(yy[t,
        ]))
      }
      if (!is.na(thresholds.low[t, k]) & sum(!is.na(y[t, ])) > 1 &
          sum(!is.na(yy[t, ])) > 5) {
        thres <- thresholds.low[t, k]
        p.y.low[t, k] <- sum(y[t, ] < thres, na.rm = T)/sum(!is.na(y[t,
        ]))
        p.yy.low[t, k] <- sum(yy[t, ] < thres, na.rm = T)/sum(!is.na(yy[t,
        ]))
      }
    }
  }
  
  
  odds.ratio.high <- array(NA, dim = c(day, length(keys)))
  odds.ratio.medium <- array(NA, dim = c(day, length(keys)))
  odds.ratio.low <- array(NA, dim = c(day, length(keys)))
  for (k in 1:length(keys)) {
    
    if (use.smoothing) {
      x <- start.t:end.t + offset.t
      y <- p.y.high[start.t:end.t + offset.t, k]
      I <- !is.na(y)
      if (sum(I) > 4) {
        model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
        p.y.high[start.t:end.t + offset.t, k] <- predict(model,
                                                         start.t:end.t + offset.t)$y
        p.y.high[start.t:end.t + offset.t, k][p.y.high[start.t:end.t +
                                                         offset.t, k] < 0] <- 0
      } else {
        p.y.high[start.t:end.t + offset.t, k] <- 0
      }
      
      x <- start.t:end.t + offset.t
      y <- p.yy.high[start.t:end.t + offset.t, k]
      I <- !is.na(y)
      if (sum(I) > 4) {
        model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
        p.yy.high[start.t:end.t + offset.t, k] <- predict(model,
                                                          start.t:end.t + offset.t)$y
        p.yy.high[start.t:end.t + offset.t, k][p.yy.high[start.t:end.t +
                                                           offset.t, k] < 0] <- 0
      } else {
        p.yy.high[start.t:end.t + offset.t, k] <- 0
      }
      
      x <- start.t:end.t + offset.t
      y <- p.y.low[start.t:end.t + offset.t, k]
      I <- !is.na(y)
      if (sum(I) > 4) {
        model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
        p.y.low[start.t:end.t + offset.t, k] <- predict(model,
                                                        start.t:end.t + offset.t)$y
        p.y.low[start.t:end.t + offset.t, k][p.y.low[start.t:end.t +
                                                       offset.t, k] < 0] <- 0
      } else {
        p.y.low[start.t:end.t + offset.t, k] <- 0
      }
      
      x <- start.t:end.t + offset.t
      y <- p.yy.low[start.t:end.t + offset.t, k]
      I <- !is.na(y)
      if (sum(I) > 4) {
        model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
        p.yy.low[start.t:end.t + offset.t, k] <- predict(model,
                                                         start.t:end.t + offset.t)$y
        p.yy.low[start.t:end.t + offset.t, k][p.yy.low[start.t:end.t +
                                                         offset.t, k] < 0] <- 0
      } else {
        p.yy.low[start.t:end.t + offset.t, k] <- 0
      }
    }
    
    if (use.laplace) {
      odds.ratio.high[, k] <- (p.y.high[, k] + 0.1)/(p.yy.high[,
                                                               k] + 0.1)
      odds.ratio.low[, k] <- (p.y.low[, k] + 0.1)/(p.yy.low[, k] +
                                                     0.1)
      
      y.high.and.low <- p.y.high[, k] + p.y.low[, k]
      y.high.and.low[y.high.and.low > 1] <- 1
      
      yy.high.and.low <- p.yy.high[, k] + p.yy.low[, k]
      yy.high.and.low[yy.high.and.low > 1] <- 1
      
      odds.ratio.medium[, k] <- (1 - y.high.and.low + 0.1)/(1 - yy.high.and.low +
                                                              0.1)
    } else {
      odds.ratio.high[, k] <- (p.y.high[, k] + 0)/(p.yy.high[, k] +
                                                     0)
      odds.ratio.low[, k] <- (p.y.low[, k] + 0)/(p.yy.low[, k] +
                                                   0)
      
      odds.ratio.medium[, k] <- (1 - p.y.high[, k] - p.y.low[, k] +
                                   0)/(1 - p.yy.high[, k] - p.yy.low[, k] + 0)
    }
    
    
    odds.ratio.high[start.t:end.t + offset.t, ][is.na(odds.ratio.high[start.t:end.t +
                                                                        offset.t, ])] <- 1
    odds.ratio.low[start.t:end.t + offset.t, ][is.na(odds.ratio.low[start.t:end.t +
                                                                      offset.t, ])] <- 1
    odds.ratio.medium[start.t:end.t + offset.t, ][is.na(odds.ratio.medium[start.t:end.t +
                                                                            offset.t, ])] <- 1
    
  }
  
  
  
  ## Test Model
  trials <- test
  
  out.odds <- array(NA, dim = c(day, length(trials)))
  is.severe <- array(NA, length(trials))
  data.points <- array(NA, length(trials))
  real.st <- array(NA, length(trials))
  
  
  for (i in 1:length(trials)) {
    trial <- trials[i]
    
    patient <- all.patients[trial]
    m <- (1:length(all.patients))[all.patients == patient]
    patient.data <- all.patient.data[[m]]
    
    
    if (engraftment.aware) {
      real.start.t <- max(df$Neutrophil.engraftment.days[df$ID ==
                                                           patient], start.t)
    } else {
      real.start.t <- start.t
    }
    
    real.st[i] <- real.start.t
    
    
    odds <- array(NA, day)
    if (sum(!is.na(patient.data)) > 50) {
      for (t in real.start.t:end.t) {
        odds[t + offset.t] <- 0
        for (k in 1:length(keys)) {
          key <- keys[k]
          val <- patient.data[t + offset.t, names(patient.data) ==
                                key]
          if (!is.na(val) & !is.na(thresholds.high[t + offset.t,
                                                   k])) {
            if (val > thresholds.high[t + offset.t, k]) {
              odds[t + offset.t] <- odds[t + offset.t] + log(odds.ratio.high[t +
                                                                               offset.t, k])
            }
          }
          if (!is.na(val) & !is.na(thresholds.low[t + offset.t,
                                                  k])) {
            if (val < thresholds.low[t + offset.t, k]) {
              odds[t + offset.t] <- odds[t + offset.t] + log(odds.ratio.low[t +
                                                                              offset.t, k])
            }
          }
          if (!is.na(val) & !is.na(thresholds.low[t + offset.t,
                                                  k]) & !is.na(thresholds.high[t + offset.t, k])) {
            if (val >= thresholds.low[t + offset.t, k] & val <=
                thresholds.high[t + offset.t, k]) {
              odds[t + offset.t] <- odds[t + offset.t] + log(odds.ratio.medium[t +
                                                                                 offset.t, k])
            }
          }
        }
      }
      
      odds.old <- odds
      for (t in 1:(end.t + offset.t)) {
        if (sum(!is.na(odds.old[max(1, t - pre.window.size + 1):t])) ==
            0) {
          odds[t] <- sum(odds.old[max(1, t - pre.window.size +
                                        1):t], na.rm = T)
        }
      }
      
      
      
      out.odds[, i] <- odds
      
      if (patient %in% all.severe) {
        is.severe[i] <- 1
      } else {
        is.severe[i] <- 0
      }
      
      data.points[i] <- data.density[all.patients == patient]
    }
  }
  
  
  
  ## Turn on use.stationary mode
  if (use.stationary) {
    
    severe.I.raw <- intersect(big.severe.I, train)
    mild.I.raw <- intersect(big.mild.I, train)
    
    
    for (t in start.t:end.t) {
      
      mild.I <- mild.I.raw
      mild.I <- setdiff(mild.I, out.pt[[t]])
      
      severe.I <- severe.I.raw
      severe.I <- setdiff(severe.I, out.pt[[t]])
      
      train.data <- subset(df, ID %in% complete.patients[c(mild.I,
                                                           severe.I)])
      
      model.nb <- naiveBayes(x = train.data[, c(4:18)], y = train.data$Severe.aGVHD.within.100.days,
                             laplace = 50, )
      
      
      for (i in 1:length(trials)) {
        
        trial <- trials[i]
        
        patient <- all.patients[trial]
        m <- (1:length(all.patients))[all.patients == patient]
        patient.data <- df[m, ]
        
        prob.nb <- predict(model.nb, patient.data[, c(4:18)], type = c("raw"),
        )
        
        out.odds[t + offset.t, i] <- sum(c(out.odds[t + offset.t,
                                                    i], log((prob.nb[2] + 0.01)/(prob.nb[1] + 0.01))), na.rm = T)
        
        
      }
    }
    
  }
  
  
  
  # Output
  out <- list(is.severe, out.odds, data.points, odds.ratio.high, odds.ratio.medium,
              odds.ratio.low)
  
  out
  
}


make.risk.model <- function(variable.number, lead.time.steps, smooth, percentage.effectors,
                            keys) {
  
  risk.model <- array(0.5, dim = c(lead.time.steps, variable.number))
  
  for (k in 1:variable.number) {
    if (k <= variable.number * percentage.effectors) {
      if (smooth) {
        temp <- array(1, lead.time.steps)
        for (t in 2:lead.time.steps) {
          temp[t] <- temp[t - 1] * rnorm(1, mean = 1, sd = 0.05)
        }
      } else {
        temp <- runif(lead.time.steps)
      }
      
      risk.model[, k] <- (temp - min(temp))/(max(temp) - min(temp)) *
        0.4 + 0.3
    }
  }
  
  risk.model <- data.frame(log(risk.model/(1 - risk.model)))
  names(risk.model) <- keys
  
  risk.model
  
}



generate.data <- function(risk.model, population.size, keys, missing.rate) {
  
  complete.patient.data <- list()
  outcome <- array(NA, population.size)
  prob.collect <- array(NA, population.size)
  target.cases <- round(percentage.cases * population.size)
  
  i <- 0
  while (i < population.size) {
    
    data <- array(NA, dim = c(lead.time.steps, variable.number))
    
    for (k in 1:variable.number) {
      temp <- array(NA, lead.time.steps)
      temp[1] <- sample(c(-1, 1), 1)
      for (t in 2:lead.time.steps) {
        if (runif(1) < 0.7) {
          temp[t] <- temp[t - 1]
        } else {
          temp[t] <- -temp[t - 1]
        }
      }
      data[, k] <- temp
    }
    data <- data.frame(data)
    names(data) <- keys
    
    total.odds <- sum(data * risk.model)
    prob <- exp(total.odds)/(1 + exp(total.odds))
    
    if (runif(1) < prob) {
      if (sum(outcome == 1, na.rm = T) < target.cases) {
        i <- i + 1
        outcome[i] <- 1
        
        if (missing.rate > 0) {
          rand.matrix <- as.array(Matrix(runif(nrow(data) * ncol(data)),
                                         nrow(data), ncol(data)))
          data[rand.matrix < missing.rate] <- NA
        }
        complete.patient.data[[i]] <- data
        prob.collect[i] <- prob
      }
    } else {
      if (sum(outcome == 0, na.rm = T) < population.size - target.cases) {
        i <- i + 1
        outcome[i] <- 0
        
        if (missing.rate > 0) {
          rand.matrix <- as.array(Matrix(runif(nrow(data) * ncol(data)),
                                         nrow(data), ncol(data)))
          data[rand.matrix < missing.rate] <- NA
        }
        complete.patient.data[[i]] <- data
        prob.collect[i] <- prob
      }
    }
  }
  
  list(complete.patient.data, outcome, prob.collect)
}



data.flatten <- function(data.list) {
  
  data <- c()
  
  for (i in 1:length(data.list)) {
    temp <- data.list[[i]]
    
    datum <- c()
    for (j in 1:nrow(temp)) {
      datum <- c(datum, as.numeric(temp[j, ]))
    }
    
    data <- rbind(data, datum)
  }
  
  Matrix(data, sparse = T)
  
}



# Define hyper parameters
variable.number <- 50
lead.time.steps <- 14
percentage.cases <- 0.15

keys <- c()
for (i in 1:variable.number) {
  keys[i] <- paste("X", i, sep = "")
}

smooth.list <- c(TRUE, FALSE)
percentage.effectors.list <- c(0.1, 0.9)
population.size.list <- c(1000)
missing.rate.list <- c(0, 0.6)

offset.t <- 15
simul.results <- list()




# Build Model
for (ttt in (1 + length(simul.results)):200) {
  
  
  simul.results[[ttt]] <- data.frame()
  
  for (m1 in 1:length(smooth.list)) {
    smooth <- smooth.list[m1]
    for (m2 in 1:length(percentage.effectors.list)) {
      percentage.effectors <- percentage.effectors.list[m2]
      for (m3 in 1:length(population.size.list)) {
        population.size <- population.size.list[m3]
        for (m4 in 1:length(missing.rate.list)) {
          missing.rate <- missing.rate.list[m4]
          
          risk.model <- make.risk.model(variable.number, lead.time.steps,
                                        smooth, percentage.effectors, keys)
          
          temp <- generate.data(risk.model, population.size, keys,
                                missing.rate)
          
          complete.patient.data <- temp[[1]]
          outcome <- temp[[2]]
          
          complete.severe <- (1:length(complete.patient.data))[outcome ==
                                                                 1]
          complete.patients <- 1:length(complete.patient.data)
          
          Severe.aGVHD.onset.day <- array(lead.time.steps + 1,
                                          length(complete.patient.data))
          Severe.aGVHD.within.100.days <- array(0, length(complete.patient.data))
          Severe.aGVHD.within.100.days[complete.severe] <- 1
          df <- data.frame(Severe.aGVHD.onset.day, Severe.aGVHD.within.100.days)
          
          pre.window.size <- lead.time.steps
          
          train <- sample(length(complete.patient.data), round(length(complete.patient.data) *
                                                                 0.8))
          test <- setdiff(complete.patients, train)
          
          day <- lead.time.steps
          data.density <- array(NA, length(complete.patients))
          
          out <- composite.xvalid(keys, train, test, start.t = 1 -
                                    offset.t, end.t = lead.time.steps - offset.t, engraftment.aware = FALSE,
                                  use.smoothing = TRUE, use.laplace = TRUE, use.stationary = FALSE)
          
          auroc <- array(NA, lead.time.steps)
          auprc <- array(NA, lead.time.steps)
          for (t in 1:lead.time.steps) {
            temp <- out[[2]][t, ]
            I <- !is.na(temp)
            auroc[t] <- roc(as.factor(test %in% complete.severe),
                            temp, direction = "<")$auc
            auprc[t] <- pr.curve(scores.class0 = as.numeric(temp[(test %in% complete.severe) & I]), 
                                 scores.class1 = as.numeric(temp[!(test %in% complete.severe) & I]), 
                                 curve = FALSE)$auc.integral
          }
          
          
          
          data.flat <- data.flatten(complete.patient.data)
          
          auroc.xgb <- array(NA, lead.time.steps)
          auprc.xgb <- array(NA, lead.time.steps)
          for (t in 1:lead.time.steps) {
            dtrain <- xgb.DMatrix(data = data.flat[train, 1:(variable.number *
                                                               t)], label = outcome[train] == 1)
            dtest <- xgb.DMatrix(data = data.flat[test, 1:(variable.number *
                                                             t)], label = outcome[test] == 1)
            xgb <- xgboost(data = dtrain, objective = "binary:logistic",
                           nround = 10, verbose = 0)
            prob.xgb <- predict(xgb, newdata = dtest)
            
            auroc.xgb[t] <- roc(outcome[test], prob.xgb, direction = "<")$auc
            auprc.xgb[t] <- pr.curve(scores.class0 = prob.xgb[test %in% complete.severe], 
                                     scores.class1 = prob.xgb[!(test %in% complete.severe)], 
                                     curve = TRUE, max.compute = TRUE, min.compute = TRUE)$auc.integral
          }
          
          simul.results[[ttt]] <- rbind(simul.results[[ttt]], c(smooth,
                                                                percentage.effectors, population.size, missing.rate,
                                                                auroc, auroc.xgb, auprc, auprc.xgb))
          cat(ttt, "\n", c(smooth, percentage.effectors, population.size,
                           missing.rate), "\n", auroc, "\n", auroc.xgb, "\n",
              auprc, "\n", auprc.xgb, "\n")
        }
      }
    }
  }
  
  save(simul.results, file = "simul-results.RData")
  
  
}




# Evaluate Model

temp <- simul.results[[1]][, -c(1:4)]
temp2 <- simul.results[[1]][, -c(1:4)]^2

for (i in 2:150) {
  temp <- temp + simul.results[[i]][, -c(1:4)]
  temp2 <- temp2 + simul.results[[i]][, -c(1:4)]^2
}

temp <- temp/length(simul.results)
temp2 <- temp2/length(simul.results)

simul.results.summary <- simul.results[[1]]
simul.results.summary[, -c(1:4)] <- temp
simul.results.summary <- data.frame(simul.results.summary, sqrt((temp2 -
                                                                   temp^2)/(length(simul.results) - 1)))

colnames(simul.results.summary) <- c("smooth.world", "percentage.effectors",
                                     "population.size", "missing.rate", 
                                     paste("goat.roc.mean.t", 1:lead.time.steps, sep = ""), 
                                     paste("xgb.roc.mean.t", 1:lead.time.steps, sep = ""),
                                     paste("goat.prc.mean.t", 1:lead.time.steps, sep = ""), 
                                     paste("xgb.prc.mean.t", 1:lead.time.steps, sep = ""), 
                                     paste("goat.roc.se.t", 1:lead.time.steps, sep = ""), 
                                     paste("xgb.roc.se.t", 1:lead.time.steps, sep = ""),
                                     paste("goat.prc.se.t", 1:lead.time.steps, sep = ""), 
                                     paste("xgb.prc.se.t", 1:lead.time.steps, sep = ""))

write.csv(simul.results.summary, file='simul-results.csv')



# Statistical test

effector <- 0.7

smooth <- 1
missing <- 0.3

roc.pvals <- c()
prc.pvals <- c()

for (t in 1:lead.time.steps) {
  
  temp0 <- c()
  temp1 <- c()
  
  for (j in 1:length(simul.results)) {
    result <- simul.results[[j]]
    temp0[j] <- result[result[, 1] == smooth & result[, 2] == effector &
                         result[, 4] == missing, 4 + t]
    temp1[j] <- result[result[, 1] == smooth & result[, 2] == effector &
                         result[, 4] == missing, 4 + lead.time.steps + t]
  }
  
  roc.pvals[t] <- t.test(temp0, temp1, alternative = "two.sided", paired = TRUE)$p.value
  
  
  temp0 <- c()
  temp1 <- c()
  for (j in 1:length(simul.results)) {
    result <- simul.results[[j]]
    temp0[j] <- result[result[, 1] == smooth & result[, 2] == effector &
                         result[, 4] == missing, 4 + lead.time.steps * 2 + t]
    temp1[j] <- result[result[, 1] == smooth & result[, 2] == effector &
                         result[, 4] == missing, 4 + lead.time.steps * 3 + t]
  }
  prc.pvals[t] <- t.test(temp0, temp1, alternative = "two.sided", paired = TRUE)$p.value
}


roc.pvals < 0.05
prc.pvals < 0.05

