rm(list = ls())



# Load libraries - when fail to import, please run:
# install.packages('library-name')
library(readxl)
library(e1071)
library(randomForestSRC)
library(randomForest)
library(dplyr)
library(pROC)
library(PRROC)
library(modEvA)
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
      real.start.t <- max(df$Neutrophil.engraftment.days[df$ID == patient], start.t)
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
      
      model.nb <- naiveBayes(x = train.data[, c(4:18)], 
                             y = train.data$Severe.aGVHD.within.100.days,
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



# Load data
root <- "/Users/cuimengxuan/Desktop/Smartphone/smartphone/"
keys <- paste("x", 1:561, sep = "")

status1 <- 4  # sitting
status2 <- 8  # from sitting to standing
lead.time.steps <- 7  # look ahead by how much

# a 'patient' here equals to a 12.8-sec segment of smartphone data
complete.patients <- c()  
# a 'severe' case equals to a 12.8-sec segment of smartphone data that is followed by postural transition
complete.severe <- c()  
# a 'mild' case equals to a 12.8-sec segment of smartphone data that is not followed by postural transition
complete.mild <- c() 
complete.patient.data <- list()
train <- c()
test <- c()


# collect training data
raw.X <- read.csv(paste(root, "Train/X_train.txt", sep = ""),
                  sep = " ", header = FALSE)
colnames(raw.X) <- keys
raw.Y <- read.csv(paste(root, "Train/Y_train.txt", sep = ""),
                  sep = " ", header = FALSE)[, 1]
subject.IDs <- read.csv(paste(root, "Train/subject_id_train.txt",
                              sep = ""), sep = " ", header = FALSE)[, 1]
for (i in 1:(nrow(raw.X) - lead.time.steps * 2)) {
  if (length(unique(subject.IDs[i + 0:(lead.time.steps + 3)])) ==
      1 & sum(raw.Y[i + 0:(lead.time.steps - 1)] == status1) ==
      lead.time.steps) {
    if (sum(raw.Y[i + lead.time.steps + 0:3] == status2) >
        0) {
      patient.id <- length(complete.patients) + 1
      complete.patients <- c(complete.patients, patient.id)
      train <- c(train, patient.id)
      complete.patient.data[[patient.id]] <- raw.X[i +
                                                     0:(lead.time.steps + 3), ]
      complete.severe <- c(complete.severe, patient.id)
    }
    if (sum(raw.Y[i + lead.time.steps + 0:3] == status2) ==
        0) {
      patient.id <- length(complete.patients) + 1
      complete.patients <- c(complete.patients, patient.id)
      train <- c(train, patient.id)
      complete.patient.data[[patient.id]] <- raw.X[i +
                                                     0:(lead.time.steps - 1), ]
      complete.mild <- c(complete.mild, patient.id)
    }
  }
}


# collect test data
raw.X <- read.csv(paste(root, "Test/X_test.txt", sep = ""), sep = " ",
                  header = FALSE)
colnames(raw.X) <- keys
raw.Y <- read.csv(paste(root, "Test/Y_test.txt", sep = ""), sep = " ",
                  header = FALSE)[, 1]
subject.IDs <- read.csv(paste(root, "Test/subject_id_test.txt",
                              sep = ""), sep = " ", header = FALSE)[, 1]
for (i in 1:(nrow(raw.X) - lead.time.steps * 2)) {
  if (length(unique(subject.IDs[i + 0:(lead.time.steps + 3)])) ==
      1 & sum(raw.Y[i + 0:(lead.time.steps - 1)] == status1) ==
      lead.time.steps) {
    if (sum(raw.Y[i + lead.time.steps + 0:3] == status2) >
        0) {
      patient.id <- length(complete.patients) + 1
      complete.patients <- c(complete.patients, patient.id)
      test <- c(test, patient.id)
      complete.patient.data[[patient.id]] <- raw.X[i +
                                                     0:(lead.time.steps + 3), ]
      complete.severe <- c(complete.severe, patient.id)
    }
    if (sum(raw.Y[i + lead.time.steps + 0:3] == status2) ==
        0) {
      patient.id <- length(complete.patients) + 1
      complete.patients <- c(complete.patients, patient.id)
      test <- c(test, patient.id)
      complete.patient.data[[patient.id]] <- raw.X[i +
                                                     0:(lead.time.steps - 1), ]
      complete.mild <- c(complete.mild, patient.id)
    }
  }
}




# Build daGOAT
Severe.aGVHD.onset.day <- array(lead.time.steps + 1,
                                length(complete.patient.data))
Severe.aGVHD.within.100.days <- array(0,
                                      length(complete.patient.data))
Severe.aGVHD.within.100.days[complete.severe] <- 1
df <- data.frame(Severe.aGVHD.onset.day, Severe.aGVHD.within.100.days)

offset.t <- 15
pre.window.size <- lead.time.steps

day <- lead.time.steps
data.density <- array(NA, length(complete.patients))
out <- composite.xvalid(keys, train, test, start.t = (1 - offset.t),
                        end.t = (lead.time.steps - offset.t), engraftment.aware = FALSE,
                        use.smoothing = TRUE, use.laplace = TRUE, use.stationary = FALSE)

auroc <- array(NA, lead.time.steps)
auprc <- array(NA, lead.time.steps)
for (t in 1:lead.time.steps) {
  
  temp <- out[[2]][t, ]
  auroc[t] <- roc(as.factor(test %in% complete.severe), temp,
                  direction = "<")$auc
  auprc[t] <- pr.curve(scores.class0 = as.numeric(temp[test %in% complete.severe]), 
                       scores.class1 = as.numeric(temp[!(test %in% complete.severe)]), curve = FALSE
                       )$auc.integral
  
}




# Build Other Models
nb.auroc <- c()
rf.auroc <- c()
xgb.auroc <- c()

nb.auprc <- c()
rf.auprc <- c()
xgb.auprc <- c()

for (t in 1:lead.time.steps) {
  train.data <- c()
  for (i in 1:length(train)) {
    train.data <- rbind(train.data, as.numeric(as.matrix(complete.patient.data[[train[i]]][1:t,
    ])))
  }
  test.data <- c()
  for (i in 1:length(test)) {
    test.data <- rbind(test.data, as.numeric(as.matrix(complete.patient.data[[test[i]]][1:t,
    ])))
  }
  
  model.nb <- naiveBayes(x = train.data, y = train %in% complete.severe)
  prob.nb <- predict(model.nb, test.data, type = c("raw"))[, 2]
  nbroc <- roc(as.factor(test %in% complete.severe), prob.nb, direction = "<")
  nb.auroc[t] <- nbroc$auc
  nb.auprc[t] <- pr.curve(scores.class0 = prob.nb[test %in% complete.severe],
                          scores.class1 = prob.nb[!(test %in% complete.severe)], curve = TRUE,
                          max.compute = TRUE, min.compute = TRUE)$auc.integral
  
  model.rf <- randomForest(x = train.data, y = as.factor(train %in% complete.severe),
                           ntree = 5, importance = TRUE, proximity = TRUE)
  prob.rf <- predict(model.rf, test.data, type = "prob")[, 2]
  rfroc <- roc(as.factor(test %in% complete.severe), prob.rf, direction = "<")
  rf.auroc[t] <- rfroc$auc
  rf.auprc[t] <- pr.curve(scores.class0 = prob.rf[test %in% complete.severe],
                          scores.class1 = prob.rf[!(test %in% complete.severe)], curve = TRUE,
                          max.compute = TRUE, min.compute = TRUE)$auc.integral
  
  
  xgbtrain <- xgb.DMatrix(data = data.matrix(train.data), 
                          label = data.matrix(as.numeric(train %in% complete.severe)))
  xgbtest <- xgb.DMatrix(data = data.matrix(test.data), 
                         label = data.matrix(as.numeric(test %in% complete.severe)))
  
  watchlist = list(train = xgbtrain, test = xgbtest)
  model <- xgb.train(data = xgbtrain, watchlist = watchlist, objective = "binary:logistic",
                     nrounds = 10, verbose = FALSE)
  
  prob.xgb <- predict(model, data.matrix(test.data))
  xgbroc <- roc(as.factor(test %in% complete.severe), prob.xgb, direction = "<")
  xgb.auroc[t] <- xgbroc$auc
  xgb.auprc[t] <- pr.curve(scores.class0 = prob.xgb[test %in% complete.severe],
                           scores.class1 = prob.xgb[!(test %in% complete.severe)], curve = TRUE,
                           max.compute = TRUE, min.compute = TRUE)$auc.integral
}

write.csv(data.frame(time = (-lead.time.steps:-1) * 1.28,
                     goat.auroc = auroc, xgb.auroc, nb.auroc, rf.auroc,
                     goat.auprc = auprc, xgb.auprc, nb.auprc, rf.auprc),
          file = "test-outcome.csv")

