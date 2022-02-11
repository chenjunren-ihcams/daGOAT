rm(list=ls())




# Load libraries - when fail to import, please run: install.packages("library-name")
library(readxl)
library(dplyr)
library(e1071)
library(randomForest)
library(xgboost)
library(pROC)
library(PRROC)




# Define hyper parameters
population <- 0
engraftment.aware <- 0

start.t <- 1
end.t <- 100 
pre.window.size <- 14
post.window.size <- 14
offset.t <- 15
day <- 100 + offset.t

CV <- 0
CV.Training.Size <- 0
CV.Key.Importance <- 0
CV.Top.Keys <- 0

HO <- 1
HO.Top.Keys <- 0

Save.Files <- 1

if (population){
    peak.t <- 23
    root <- 'pediatric-demo/'
    ntree <- 20
    mtry <- 10
    keys.file <- paste(root, 'keys.csv', sep = "")
    keys <- read.csv(keys.file)$keys
    nfold <- 8
} else{
    peak.t <- 10
    root <- 'pediatric-demo/'
    ntree <- 10
    mtry <- 10
    keys.file <- paste(root, 'keys.csv', sep = "")
    keys <- read.csv(keys.file)$keys
    nfold <- 3
}

factor.col <- c("Severe.aGVHD.within.100.days", 
                "Sex", "Primary.disease", "Conditioning.regimen", 
                "Sex.match", "ABO.match", "Stem.cell.source", "Donor.type",
                "Year", "GVHD.prophylaxis", "Antithymocyte.globulin.in.conditioning", 
                "Death")




# Define functions
colApply <- function(dat, cols = colnames(dat), func = as.factor) {
    dat[cols] <- lapply(dat[cols], func)
    return(dat)
}



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



composite.xvalid <- function(keys, 
                             train, 
                             test, 
                             start.t=start.t, 
                             end.t=end.t,  
                             engraftment.aware=engraftment.aware, 
                             use.smoothing=TRUE, 
                             use.laplace=TRUE,
                             use.stationary=TRUE
                             ) {
    
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
      if (df$Severe.aGVHD.within.100.days[m]==1) {
        onset.day <- df$Severe.aGVHD.onset.day[m]
        all.patient.data[[m]][1 : min(onset.day - pre.window.size - 1, end.t) + offset.t, ] <- NA 
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
                p.y.over <- sum(y[t, ] > thres, na.rm = T)/sum(!is.na(y[t, ]))
                p.yy.over <- sum(yy[t, ] > thres, na.rm = T)/sum(!is.na(yy[t, ]))
                
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
    
    
    ## Build Model - odds ratio
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
            if (!is.na(thresholds.high[t, k]) & sum(!is.na(y[t, ])) > 1 & sum(!is.na(yy[t, ])) > 5) {
                thres <- thresholds.high[t, k]
                p.y.high[t, k] <- sum(y[t, ] > thres, na.rm = T)/sum(!is.na(y[t, ]))
                p.yy.high[t, k] <- sum(yy[t, ] > thres, na.rm = T)/sum(!is.na(yy[t, ]))
            }
            if (!is.na(thresholds.low[t, k]) & sum(!is.na(y[t, ])) > 1 & sum(!is.na(yy[t, ])) > 5) {
                thres <- thresholds.low[t, k]
                p.y.low[t, k] <- sum(y[t, ] < thres, na.rm = T)/sum(!is.na(y[t, ]))
                p.yy.low[t, k] <- sum(yy[t, ] < thres, na.rm = T)/sum(!is.na(yy[t, ]))
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
                p.y.high[start.t:end.t + offset.t, k] <- predict(model, start.t:end.t + offset.t)$y
                p.y.high[start.t:end.t + offset.t, k][p.y.high[start.t:end.t + offset.t, k] < 0] <- 0
            } else {
                p.y.high[start.t:end.t + offset.t, k] <- 0
            }
            
            x <- start.t:end.t + offset.t
            y <- p.yy.high[start.t:end.t + offset.t, k]
            I <- !is.na(y)
            if (sum(I) > 4) {
                model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
                p.yy.high[start.t:end.t + offset.t, k] <- predict(model, start.t:end.t + offset.t)$y
                p.yy.high[start.t:end.t + offset.t, k][p.yy.high[start.t:end.t + offset.t, k] < 0] <- 0
            } else {
                p.yy.high[start.t:end.t + offset.t, k] <- 0
            }
            
            x <- start.t:end.t + offset.t
            y <- p.y.low[start.t:end.t + offset.t, k]
            I <- !is.na(y)
            if (sum(I) > 4) {
                model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
                p.y.low[start.t:end.t + offset.t, k] <- predict(model, start.t:end.t + offset.t)$y
                p.y.low[start.t:end.t + offset.t, k][p.y.low[start.t:end.t + offset.t, k] < 0] <- 0
            } else {
                p.y.low[start.t:end.t + offset.t, k] <- 0
            }
            
            x <- start.t:end.t + offset.t
            y <- p.yy.low[start.t:end.t + offset.t, k]
            I <- !is.na(y)
            if (sum(I) > 4) {
                model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
                p.yy.low[start.t:end.t + offset.t, k] <- predict(model, start.t:end.t + offset.t)$y
                p.yy.low[start.t:end.t + offset.t, k][p.yy.low[start.t:end.t + offset.t, k] < 0] <- 0
            } else {
                p.yy.low[start.t:end.t + offset.t, k] <- 0
            }
        }
        
        if (use.laplace) {
            odds.ratio.high[, k] <- (p.y.high[, k] + 0.1)/(p.yy.high[, k] + 0.1)
            odds.ratio.low[, k] <- (p.y.low[, k] + 0.1)/(p.yy.low[, k] + 0.1)
            
            y.high.and.low <- p.y.high[, k] + p.y.low[, k]
            y.high.and.low[y.high.and.low > 1] <- 1
            
            yy.high.and.low <- p.yy.high[, k] + p.yy.low[, k]
            yy.high.and.low[yy.high.and.low > 1] <- 1
            
            odds.ratio.medium[, k] <- (1 - y.high.and.low + 0.1)/(1 - yy.high.and.low + 0.1)
        } else {
            odds.ratio.high[, k] <- (p.y.high[, k] + 0)/(p.yy.high[, k] + 0)
            odds.ratio.low[, k] <- (p.y.low[, k] + 0)/(p.yy.low[, k] + 0)
            
            odds.ratio.medium[, k] <- (1 - p.y.high[, k] - p.y.low[, k] + 0)/(1 - p.yy.high[, k] - p.yy.low[, k] + 0)
        }
        
        
        odds.ratio.high[start.t:end.t + offset.t, ][is.na(odds.ratio.high[start.t:end.t + offset.t, ])] <- 1
        odds.ratio.low[start.t:end.t + offset.t, ][is.na(odds.ratio.low[start.t:end.t + offset.t, ])] <- 1
        odds.ratio.medium[start.t:end.t + offset.t, ][is.na(odds.ratio.medium[start.t:end.t + offset.t, ])] <- 1
        
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
        real.start.t <- max(df$`Neutrophil.engraftment.days`[df$ID == patient], start.t)
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
            val <- patient.data[t + offset.t, names(patient.data) == key]
            if (!is.na(val) & !is.na(thresholds.high[t + offset.t, k])) {
              if (val > thresholds.high[t + offset.t, k]) {
                odds[t + offset.t] <- odds[t + offset.t] + log(odds.ratio.high[t + offset.t, k])
              }
            }
            if (!is.na(val) & !is.na(thresholds.low[t + offset.t, k])) {
              if (val < thresholds.low[t + offset.t, k]) {
                odds[t + offset.t] <- odds[t + offset.t] + log(odds.ratio.low[t + offset.t, k])
              }
            }
            if (!is.na(val) & !is.na(thresholds.low[t + offset.t, k]) & !is.na(thresholds.high[t + offset.t, k])) {
              if (val >= thresholds.low[t + offset.t, k] & val <= thresholds.high[t + offset.t, k]) {
                odds[t + offset.t] <- odds[t + offset.t] + log(odds.ratio.medium[t + offset.t, k])
              }
            }
          }
        }
        
        odds.old <- odds 
        for (t in 1:(end.t + offset.t)) {
          if (sum(!is.na(odds.old[max(1, t - pre.window.size + 1) : t]))==0) {
            odds[t] <- sum(odds.old[max(1, t - pre.window.size + 1) : t], na.rm=T) 
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
        
        train.data <- subset(df, ID %in% complete.patients[c(mild.I, severe.I)])
        
        model.nb <- naiveBayes(x=train.data[, c(4:18)], 
                               y=train.data$Severe.aGVHD.within.100.days, 
                               laplace=10)
        
        
        for (i in 1:length(trials)) {
          
          trial <- trials[i]
          
          patient <- all.patients[trial]
          m <- (1:length(all.patients))[all.patients == patient]
          patient.data <- df[m, ]
          
          prob.nb <- predict(model.nb, patient.data[, c(4:18)], type = c("raw"),)
          
          out.odds[t + offset.t, i] <- sum(c(out.odds[t + offset.t, i],
                                             log( (prob.nb[2]+0.01) / (prob.nb[1]+0.01) )
                                             ), na.rm=T)  
          
        }
        
      }
      
    }
    
    
    
    ## Output
    out <- list(is.severe, out.odds, data.points, 
                odds.ratio.high, odds.ratio.medium, odds.ratio.low)
    
    out

}




# Internal Validation 
if(CV){
  
  
  
  # Load and wrangle data.
  encoded.file <- paste(root,'encoded.csv', sep = "")
  df <- read.csv(encoded.file)
  df <- df[df$fold<=nfold, ]
  df <- df[, c(-2, -3)]
  df <- colApply(df, factor.col)
  
  complete.patients <- df$ID
  death.or.not <- df$Death
  death.date <- df$Last.followup.time.for.death
  fold <- c(df$fold)
  
  complete.severe <- df$ID[df$Severe.aGVHD.within.100.days == 1]
  complete.mild <- df$ID[df$Severe.aGVHD.within.100.days == 0]
  
  
  data.density <- c()
  complete.patient.data <- list()
  for (i in 1:length(complete.patients)) {
    
    patient <- complete.patients[i]
    file.name <- paste(root, "processed/pt_", patient, ".csv", sep = "")
    patient.data <- read.csv(file.name, fileEncoding = "gbk")
    patient.data <- patient.data[1:day, 2:ncol(patient.data)]
    
    # time-limited sample-and-hold
    for (w in 1:3) {   
      hasValue <- list()
      for (t in offset.t:day) {
        hasValue[[t]] <- !is.na(patient.data[t, ])
      }
      
      for (t in day:offset.t) {
        patient.data[t, !hasValue[[t]]] <- patient.data[t - 1, !hasValue[[t]]]
      }
    }
    
    complete.patient.data[[i]] <- patient.data
    
    
    j <- (1:nrow(df))[df$ID == patient] 
    
    # ignore data measured after severe aGVHD onset
    if (df$Severe.aGVHD.within.100.days[j] == 1) { 
      complete.patient.data[[i]][(df$Severe.aGVHD.onset.day[j]):end.t + offset.t, ] <- NA
      cat(".")
    }
    
    # ignore data measured after death date
    if (df$Death[j] == 1) { 
      complete.patient.data[[i]][(df$Last.followup.time.for.death[j]):end.t + offset.t, ] <- NA
      cat(".")
    }
    
    # ignore data measured before neutrophil engraftment if run in the 'engraftment-aware' mode
    if (engraftment.aware) { 
      j <- (1:nrow(df))[df$ID == patient]
      complete.patient.data[[i]][1:(df$Neutrophil.engraftment.days[j] + 14), ] <- NA
      cat(".")
    }
    
    complete.patient.data[[i]] <- complete.patient.data[[i]][1:day, ]
    data.density[i] <- sum(!is.na(complete.patient.data[[i]][, keys]))

  }
  
  
  big.severe.I <- c()
  for (i in 1:length(complete.severe)) {
    big.severe.I <- c(big.severe.I, (1:length(complete.patients))[complete.patients == complete.severe[i]])
  }
  big.severe.I <- big.severe.I[!is.na(big.severe.I)]
  big.mild.I <- setdiff(1:length(complete.patients), big.severe.I)
  
  
  out.pt <- list()
  for (t in start.t:end.t) {
    out1 <- (1:length(complete.patients))[ (df$Severe.aGVHD.within.100.days == 1) & (df$Severe.aGVHD.onset.day <= t) ]
    out2 <- (1:length(complete.patients))[ (df$Death == 1) & (df$Last.followup.time.for.death <= t) ]
    out.pt[[t]] <- union(out1, out2)
  }

  

  ## model -- daGOAT
  auc.daGOAT <- c()
  prc.daGOAT <- c()
  peak.daGOAT <- c()
  for (f in 1 : nfold) {
    
    test.set <- (1:length(complete.patients))[fold==f]
    severe.I <- setdiff(big.severe.I, test.set)
    mild.I <- setdiff(big.mild.I, test.set)
    
    out <- composite.xvalid(keys=keys, 
                            train=c(severe.I, mild.I), 
                            test=test.set, 
                            start.t=start.t, 
                            end.t=end.t,  
                            engraftment.aware=engraftment.aware, 
                            use.smoothing=TRUE, 
                            use.laplace=TRUE,
                            use.stationary=TRUE)
    
    test.set.start <- array(NA, length(test.set))
    for (i in 1:length(test.set)) {
      patient <- complete.patients[test.set[i]]
      test.set.start[i] <- df$Severe.aGVHD.onset.day[df$ID == patient]
    }
    test.set.start[is.na(test.set.start)] <- Inf
    
    auc <- array(NA, end.t)
    prc <- array(NA, end.t)
    
    for (t in start.t:end.t) {
      
      res <- data.frame(
        daGOAT.id = complete.patients[test.set],
        daGOAT.death = death.or.not[test.set], 
        daGOAT.death.date = death.date[test.set],
        ground.truth = as.factor(out[[1]]),
        Severe.aGVHD.time = test.set.start,
        pred = out[[2]][t + offset.t, ] 
      )
      
      res <- res[! (res$ground.truth==1 & res$Severe.aGVHD.time <= t), ]
      res <- res[! (res$daGOAT.death==1 & res$daGOAT.death.date <= t), ]
      res <- res[!is.na(res$pred), ]
      
      res$tw <- res$Severe.aGVHD.time - t
      res$ground.truth.fix <- res$ground.truth
      res$ground.truth.fix[ (res$ground.truth==1) & (res$tw > post.window.size)] <- 0 
      
      if(sd(as.numeric(res$ground.truth.fix)) != 0){
        
        auc.list <- roc(response=res$ground.truth.fix, predictor=res$pred, direction='<')
        auc[t] <- auc.list$auc
        
        if(sd(as.numeric(res$pred)) != 0){
          prc.list <- pr.curve(scores.class0=res[res$ground.truth.fix==1, ]$pred,
                               scores.class1=res[res$ground.truth.fix==0, ]$pred,
                               curve=TRUE, max.compute=TRUE, min.compute=TRUE)
          prc[t] <- prc.list$auc.integral
        } else{
          prc[t] <- length(res[res$ground.truth.fix==1, ]$pred) / dim(res)[1]
        }
      
      }
      
      if(t==peak.t){
        res$fold <- f
        peak.daGOAT <- rbind(peak.daGOAT, res)
      }
      
    }
    
    auc.daGOAT <- rbind(auc.daGOAT, auc)
    prc.daGOAT <- rbind(prc.daGOAT, prc)
    
  }
  
  
  
  ## Model -- Stationary Features Naive Bayes 
  auc.NB <- c()
  prc.NB <- c()
  peak.NB <- c()
  for (f in 1 : nfold) {
    
    test.set.raw <- (1:length(complete.patients))[fold==f]
    severe.I.raw <- setdiff(big.severe.I, test.set.raw)
    mild.I.raw <- setdiff(big.mild.I, test.set.raw)
    
    auc <- array(NA, end.t)
    prc <- array(NA, end.t)
    
    for (t in start.t:end.t) {
      
      test.set <- test.set.raw
      test.set <- setdiff(test.set, out.pt[[t]])
      
      mild.I <- mild.I.raw
      mild.I <- setdiff(mild.I, out.pt[[t]])
      
      severe.I <- severe.I.raw
      severe.I <- setdiff(severe.I, out.pt[[t]])
      
      train.data <- subset(df, ID %in% complete.patients[c(mild.I, severe.I)])
      test.data <- subset(df, ID %in% complete.patients[test.set])
      
      test.data$tw <- test.data$Severe.aGVHD.onset.day - t
      test.data$Severe.aGVHD.within.100.days.fix <- test.data$Severe.aGVHD.within.100.days
      test.data$Severe.aGVHD.within.100.days.fix[ (test.data$Severe.aGVHD.within.100.days==1) & 
                                                  (test.data$tw > Inf) ] <- 0
      
      if ( ( sd (as.numeric(train.data$Severe.aGVHD.within.100.days)) != 0 ) &
           ( sum(as.numeric(train.data$Severe.aGVHD.within.100.days)-1 ) > 1 ) & 
           ( sd (as.numeric(test.data$Severe.aGVHD.within.100.days.fix)) != 0 ) 
      ){
        model <- naiveBayes(x=train.data[, c(4:18)], y=train.data$Severe.aGVHD.within.100.days, laplace=10)
        prob <- predict(model, test.data[, c(4:18)], type = c("raw"))
        
        auc.list <- roc(response=test.data$Severe.aGVHD.within.100.days.fix,
                        predictor=log( (prob[, 2] + 0.01 ) / (prob[, 1] + 0.01) ),
                        direction = "<")
        auc[t] <- auc.list$auc
        
        temp <- data.frame(response=auc.list$response, pred=auc.list$predictor, nb.id=test.data$ID)
        prc.list <- pr.curve(scores.class0=temp[temp$response==1, ]$pred,
                             scores.class1=temp[temp$response==0, ]$pred,
                             curve=TRUE, max.compute=TRUE, min.compute=TRUE)
        prc[t] <- prc.list$auc.integral
      }
      
      if(t==peak.t){
        temp$fold <- f
        peak.NB <- rbind(peak.NB, temp)
      }
      
    }
    
    auc.NB <- rbind(auc.NB, auc)
    prc.NB <- rbind(prc.NB, prc)
    
  }
  
  
  
  ## Model -- Stationary Features Random Forest
  auc.RF <- c()
  prc.RF <- c()
  peak.RF <- c()
  for (f in 1 : nfold) {
    
    test.set.raw <- (1:length(complete.patients))[fold==f]
    severe.I.raw <- setdiff(big.severe.I, test.set.raw)
    mild.I.raw <- setdiff(big.mild.I, test.set.raw)
    
    auc <- array(NA, end.t)
    prc <- array(NA, end.t)
    
    for (t in start.t:end.t) {
      
      test.set <- test.set.raw
      test.set <- setdiff(test.set, out.pt[[t]])
      
      mild.I <- mild.I.raw
      mild.I <- setdiff(mild.I, out.pt[[t]])
      
      severe.I <- severe.I.raw
      severe.I <- setdiff(severe.I, out.pt[[t]])
      
      
      train.data <- subset(df, ID %in% complete.patients[c(mild.I, severe.I)])
      test.data <- subset(df, ID %in% complete.patients[test.set])
      
      test.data$tw <- test.data$Severe.aGVHD.onset.day - t
      test.data$Severe.aGVHD.within.100.days.fix <- test.data$Severe.aGVHD.within.100.days
      test.data$Severe.aGVHD.within.100.days.fix[ (test.data$Severe.aGVHD.within.100.days==1) &
                                                  (test.data$tw > Inf) ] <- 0
      
      
      if ( ( sd (as.numeric(train.data$Severe.aGVHD.within.100.days)) != 0 ) &
           ( sum(as.numeric(train.data$Severe.aGVHD.within.100.days) - 1) > 1 ) & 
           ( sd (as.numeric(test.data$Severe.aGVHD.within.100.days.fix)) != 0 ) 
      ){
        set.seed(100)
        model <- randomForest(Severe.aGVHD.within.100.days ~ Sex + Age + BMI + 
                                Primary.disease + Conditioning.regimen + Sex.match + ABO.match + 
                                TNC.count + CD34.count + Stem.cell.source + Donor.type + 
                                Year + GVHD.prophylaxis + Antithymocyte.globulin.in.conditioning + 
                                HLA.mismatch,
                              data=train.data,  ntree=ntree,  mtry=mtry, importance=TRUE, proximity = TRUE)
        prob <- predict(model, test.data[, c(4:18)], type = "prob")[, 2]
        
        auc.list <- roc(response=test.data$Severe.aGVHD.within.100.days.fix,
                        predictor=prob,
                        direction = "<")
        auc[t] <- auc.list$auc
        
        temp <- data.frame(response=auc.list$response, pred=auc.list$predictor, rf.id=test.data$ID)
        
        if(sd(as.numeric(temp$pred)) != 0){
          prc.list <- pr.curve(scores.class0=temp[temp$response==1, ]$pred,
                               scores.class1=temp[temp$response==0, ]$pred,
                               curve=TRUE, max.compute=TRUE, min.compute=TRUE)
          prc[t] <- prc.list$auc.integral
        } else{
          prc[t] <- length(temp[temp$response==1, ]$pred) / dim(temp)[1]
        }
        
        if(t==peak.t){
          temp$fold <- f
          peak.RF <- rbind(peak.RF, temp)
        }
        
      }
      
    }
    
    auc.RF <- rbind(auc.RF, auc)
    prc.RF <- rbind(prc.RF, prc)
    
  }
  
  
  
  ## Model -- XGBoost
  dynamic.data <- list()
  for (t in start.t:end.t) {
    temp <- c()
    for (j in 1:length(complete.patients)) {
      temp <- rbind(temp, complete.patient.data[[j]][t + offset.t, keys])
    }
    dynamic.data[[t]] <- temp
  }
  
  
  dynamic.data.tw <- list()
  for (t in start.t : end.t) {
    
    temp <- dynamic.data[[t]]
    
    if(t>1){
      for (u in (t-1) : max(1, t-pre.window.size+1)){
        temp <- cbind(temp, dynamic.data[[u]])
      }
    }
    
    temp <- cbind(temp, df[, c(4:18)])
    temp$num.feat <- dim(temp)[2]
    temp$Severe.aGVHD.onset.day <- df$Severe.aGVHD.onset.day
    temp$Severe.aGVHD.within.100.days <- df$Severe.aGVHD.within.100.days
    temp$Severe.aGVHD.within.100.days.fix <- df$Severe.aGVHD.within.100.days
    temp$tw <- temp$Severe.aGVHD.onset.day - t
    temp$Severe.aGVHD.within.100.days.fix[ (temp$Severe.aGVHD.within.100.days == 1) & 
                                           (temp$tw > post.window.size) ] <- 0
    
    dynamic.data.tw[[t]] <- temp
  }
  
  
  auc.XGB <- c()
  prc.XGB <- c()
  peak.XGB <- c()
  for (f in 1 : nfold) {
    
    test.set.raw <- (1:length(complete.patients))[fold==f]
    severe.I.raw <- setdiff(big.severe.I, test.set.raw)
    mild.I.raw <- setdiff(big.mild.I, test.set.raw)
    
    auc <- array(NA, end.t)
    prc <- array(NA, end.t)
    
    for (t in start.t:end.t) {
      
      test.set <- test.set.raw
      test.set <- setdiff(test.set, out.pt[[t]])
      
      mild.I <- mild.I.raw
      mild.I <- setdiff(mild.I, out.pt[[t]])
      
      severe.I <- severe.I.raw
      severe.I <- setdiff(severe.I, out.pt[[t]])
      
      train.data <- dynamic.data.tw[[t]][c(mild.I, severe.I), ]
      test.data <- dynamic.data.tw[[t]][c(test.set), ]
      
      feat <- dynamic.data.tw[[t]]$num.feat[1]
      
      train.x <- data.matrix(train.data[, c(1:feat)])
      train.y <- data.matrix(train.data$Severe.aGVHD.within.100.days.fix)
      
      test.x <- data.matrix(test.data[, c(1:feat)])
      test.y <- data.matrix(test.data$Severe.aGVHD.within.100.days.fix)
      
      xgbtrain <- xgb.DMatrix(data=train.x, label=train.y)
      xgbtest <- xgb.DMatrix(data=test.x, label=test.y)
      
      watchlist = list(train=xgbtrain, test=xgbtest)
      
      set.seed(200)
      model <- xgb.train(data=xgbtrain, 
                         watchlist=watchlist, 
                         objective = "binary:logistic",
                         nrounds=10,
                         verbose=FALSE)
      prob <- predict(model, test.x)
      
      res <- data.frame(pred=prob, 
                        ground.truth=test.y, 
                        xgb.id=complete.patients[test.set])
      if(sd(as.numeric(res$ground.truth)) != 0){
        
        auc.list <- roc(response=res$ground.truth, predictor=res$pred, direction='<')
        auc[t] <- auc.list$auc
        
        if(sd(as.numeric(res$pred)) != 0){
          prc.list <- pr.curve(scores.class0=res[res$ground.truth==1, ]$pred,
                               scores.class1=res[res$ground.truth==0, ]$pred,
                               curve=TRUE, max.compute=TRUE, min.compute=TRUE)
          prc[t] <- prc.list$auc.integral
        } else{
          prc[t] <- length(res[res$ground.truth==1, ]$pred) / dim(res)[1]
        }
        
      }
      
      if(t==peak.t){
        res$fold <- f
        peak.XGB <- rbind(peak.XGB, res)
      }
      
    }
    
    auc.XGB <- rbind(auc.XGB, auc)
    prc.XGB <- rbind(prc.XGB, prc)
    
  }
  
  
  
  ## Output
  if(Save.Files){
    auc.cv <- rbind(auc.daGOAT, auc.NB, auc.RF, auc.XGB)
    prc.cv <- rbind(prc.daGOAT, prc.NB, prc.RF, prc.XGB)
    auc.cv.mean <- rbind( colMeans(auc.daGOAT, na.rm=1),
                          colMeans(auc.NB, na.rm=1),
                          colMeans(auc.RF, na.rm=1),
                          colMeans(auc.XGB, na.rm=1) )
    prc.cv.mean <- rbind( colMeans(prc.daGOAT, na.rm=1),
                          colMeans(prc.NB, na.rm=1),
                          colMeans(prc.RF, na.rm=1),
                          colMeans(prc.XGB, na.rm=1) )
    
    peak.label.pred <- cbind(peak.daGOAT$daGOAT.id, peak.daGOAT$ground.truth.fix, peak.daGOAT$pred, peak.daGOAT$fold, 
                             peak.NB$nb.id, peak.NB$response, peak.NB$pred, peak.NB$fold,
                             peak.RF$rf.id, peak.RF$response, peak.RF$pred, peak.RF$fold,
                             peak.XGB$xgb.id, peak.XGB$ground.truth, peak.XGB$pred, peak.XGB$fold)
    colnames(peak.label.pred) <- c('daGOAT.ID', 'daGOAT.label', 'daGOAT.pred', 'daGOAT.fold',
                                   'NB.ID', 'NB.label', 'NB.pred', 'NB.fold',
                                   'RF.ID', 'RF.label', 'RF.pred', 'RF.fold',
                                   'XGB.ID', 'XGB.label', 'XGB.pred', 'XGB.fold')
    
    write.csv(rbind(auc.cv, prc.cv), paste(root,'auc_prc_cv.csv', sep = ""))
    write.csv(rbind(auc.cv.mean, prc.cv.mean), paste(root,'auc_prc_cv_mean.csv', sep = ""))
    write.csv(peak.label.pred, paste(root,'peak_cv.csv', sep = ""))
  }
  
  

}




# Data diversity was measured in training size (randomly selected)
if(CV.Training.Size){
  
  
  auc.ts <- c()
  prc.ts <- c()
  
  for(num.train in (1:5)){
    
    for (f in 1 : nfold) {
  
      test.set <- (1:length(complete.patients))[fold==f]
      severe.I <- setdiff(big.severe.I, test.set)
      mild.I <- setdiff(big.mild.I, test.set)
      
      train <- c(severe.I, mild.I)
      
      set.seed(100)
      q <- floor(length(train)/5)
      r <- length(train)%%5
      rand.indicator.pt <- sample(1:5, r)
      for (w in 1:q) {
        rand.indicator.pt <- c(rand.indicator.pt, sample(1:5, 5))
      }
      
      out <- composite.xvalid(keys=keys, 
                              train=train[rand.indicator.pt<=num.train], 
                              test=test.set, 
                              start.t=start.t, 
                              end.t=end.t,  
                              engraftment.aware=engraftment.aware, 
                              use.smoothing=TRUE, 
                              use.laplace=TRUE,
                              use.stationary=TRUE)
      
      test.set.start <- array(NA, length(test.set))
      for (i in 1:length(test.set)) {
        patient <- complete.patients[test.set[i]]
        test.set.start[i] <- df$Severe.aGVHD.onset.day[df$ID == patient]
      }
      test.set.start[is.na(test.set.start)] <- Inf
      
      auc <- array(NA, end.t)
      prc <- array(NA, end.t)
      
      for (t in start.t:end.t) {
        
        res <- data.frame(
          daGOAT.id = complete.patients[test.set],
          daGOAT.death = death.or.not[test.set], 
          daGOAT.death.date = death.date[test.set],
          ground.truth = as.factor(out[[1]]),
          Severe.aGVHD.time = test.set.start,
          pred = out[[2]][t + offset.t, ]
        )
        
        res <- res[! (res$ground.truth==1 & res$Severe.aGVHD.time <= t), ]
        res <- res[! (res$daGOAT.death==1 & res$daGOAT.death.date <= t), ]
        res <- res[!is.na(res$pred), ]
        
        res$tw <- res$Severe.aGVHD.time - t
        res$ground.truth.fix <- res$ground.truth
        res$ground.truth.fix[ (res$ground.truth==1) & (res$tw > post.window.size)] <- 0 
        
        if(sd(as.numeric(res$ground.truth.fix)) != 0){
          auc.list <- roc(response=res$ground.truth.fix, predictor=res$pred, direction='<')
          auc[t] <- auc.list$auc
          
          if(sd(as.numeric(res$pred)) != 0){
            prc.list <- pr.curve(scores.class0=res[res$ground.truth.fix==1, ]$pred,
                                 scores.class1=res[res$ground.truth.fix==0, ]$pred,
                                 curve=TRUE, max.compute=TRUE, min.compute=TRUE)
            prc[t] <- prc.list$auc.integral
          } else{
            prc[t] <- length(res[res$ground.truth.fix==1, ]$pred) / dim(res)[1]
          }
        }
        
      }
      
      auc.ts <- rbind(auc.ts, auc)
      prc.ts <- rbind(prc.ts, prc)
      
    }
    
  }
  
  
  if(Save.Files){
    auc.ts.mean <- rbind( colMeans(auc.ts[1:nfold, ], na.rm=T), 
                          colMeans(auc.ts[1:nfold + 1 * nfold, ], na.rm=T), 
                          colMeans(auc.ts[1:nfold + 2 * nfold, ], na.rm=T), 
                          colMeans(auc.ts[1:nfold + 3 * nfold, ], na.rm=T), 
                          colMeans(auc.ts[1:nfold + 4 * nfold, ], na.rm=T)
    )
    
    prc.ts.mean <- rbind( colMeans(prc.ts[1:nfold, ], na.rm=T), 
                          colMeans(prc.ts[1:nfold + 1 * nfold, ], na.rm=T), 
                          colMeans(prc.ts[1:nfold + 2 * nfold, ], na.rm=T), 
                          colMeans(prc.ts[1:nfold + 3 * nfold, ], na.rm=T), 
                          colMeans(prc.ts[1:nfold + 4 * nfold, ], na.rm=T)
    )
    
    write.csv(auc.ts, paste(root,'auc_ts.csv', sep = ""))
    write.csv(prc.ts, paste(root,'prc_ts.csv', sep = ""))
    write.csv(auc.ts.mean, paste(root,'auc_ts_mean.csv', sep = ""))
    write.csv(prc.ts.mean, paste(root,'prc_ts_mean.csv', sep = ""))
  }
  
  
}




# Key Importance was measured by AUPRC reduction when one key was not available for the daGOAT
if(CV.Key.Importance){
  
  
  # dynamic keys
  for (f in 1 : nfold) {
    
    auc.keys <- c()
    prc.keys <- c()
    
    test.set <- (1:length(complete.patients))[fold==f]
    severe.I <- setdiff(big.severe.I, test.set)
    mild.I <- setdiff(big.mild.I, test.set)
    
    for(k in 1:length(keys)){
      
      out <- composite.xvalid(keys=keys[-k], 
                              train=c(severe.I, mild.I), 
                              test=test.set, 
                              start.t=start.t, 
                              end.t=end.t,  
                              engraftment.aware=engraftment.aware, 
                              use.smoothing=TRUE, 
                              use.laplace=TRUE,
                              use.stationary=TRUE)
      
      test.set.start <- array(NA, length(test.set))
      for (i in 1:length(test.set)) {
        patient <- complete.patients[test.set[i]]
        test.set.start[i] <- df$Severe.aGVHD.onset.day[df$ID == patient]
      }
      test.set.start[is.na(test.set.start)] <- Inf
      
      auc <- array(NA, end.t)
      prc <- array(NA, end.t)
      
      for (t in start.t:end.t) {
        
        res <- data.frame(
          daGOAT.id = complete.patients[test.set],
          daGOAT.death = death.or.not[test.set], 
          daGOAT.death.date = death.date[test.set],
          ground.truth = as.factor(out[[1]]),
          Severe.aGVHD.time = test.set.start,
          pred = out[[2]][t + offset.t, ]
        )
        
        res <- res[! (res$ground.truth==1 & res$Severe.aGVHD.time <= t), ]
        res <- res[! (res$daGOAT.death==1 & res$daGOAT.death.date <= t), ]
        res <- res[!is.na(res$pred), ]
        
        res$tw <- res$Severe.aGVHD.time - t
        res$ground.truth.fix <- res$ground.truth
        res$ground.truth.fix[ (res$ground.truth==1) & (res$tw > post.window.size)] <- 0 
        
        if(sd(as.numeric(res$ground.truth.fix)) != 0){
          auc.list <- roc(response=res$ground.truth.fix, predictor=res$pred, direction='<')
          auc[t] <- auc.list$auc
          
          prc.list <- pr.curve(scores.class0=res[res$ground.truth.fix==1, ]$pred,
                               scores.class1=res[res$ground.truth.fix==0, ]$pred, 
                               curve=TRUE, max.compute=TRUE, min.compute=TRUE)
          prc[t] <- prc.list$auc.integral
        }
        
      }
      
      auc.keys <- rbind(auc.keys, auc)
      prc.keys <- rbind(prc.keys, prc)
      
    }
    
    if(Save.Files){
      write.csv(auc.keys, paste(root,'auc_TOP.csv', sep = ""))
      write.csv(prc.keys, paste(root,'prc_TOP.csv', sep = ""))
    }
    
  }
  
  
  # stationary keys
  auc.stationary.free <- c()
  prc.stationary.free <- c()
  
  for (f in 1 : nfold) {
    
    test.set <- (1:length(complete.patients))[fold==f]
    severe.I <- setdiff(big.severe.I, test.set)
    mild.I <- setdiff(big.mild.I, test.set)
    
    out <- composite.xvalid(keys=keys, 
                            train=c(severe.I, mild.I), 
                            test=test.set, 
                            start.t=start.t, 
                            end.t=end.t,  
                            engraftment.aware=engraftment.aware, 
                            use.smoothing=TRUE, 
                            use.laplace=TRUE,
                            use.stationary=FALSE)
    
    test.set.start <- array(NA, length(test.set))
    for (i in 1:length(test.set)) {
      patient <- complete.patients[test.set[i]]
      test.set.start[i] <- df$Severe.aGVHD.onset.day[df$ID == patient]
    }
    test.set.start[is.na(test.set.start)] <- Inf
    
    auc <- array(NA, end.t)
    prc <- array(NA, end.t)
    for (t in start.t:end.t) {
      
      res <- data.frame(
        daGOAT.id = complete.patients[test.set],
        daGOAT.death = death.or.not[test.set], 
        daGOAT.death.date = death.date[test.set],
        ground.truth = as.factor(out[[1]]),
        Severe.aGVHD.time = test.set.start,
        pred = out[[2]][t + offset.t, ] 
      )
      
      res <- res[! (res$ground.truth==1 & res$Severe.aGVHD.time <= t), ]
      res <- res[! (res$daGOAT.death==1 & res$daGOAT.death.date <= t), ]
      res <- res[!is.na(res$pred), ]
      
      res$tw <- res$Severe.aGVHD.time - t
      res$ground.truth.fix <- res$ground.truth
      res$ground.truth.fix[ (res$ground.truth==1) & (res$tw > post.window.size)] <- 0 
      
      if(sd(as.numeric(res$ground.truth.fix)) != 0){
        
        auc.list <- roc(response=res$ground.truth.fix, predictor=res$pred, direction='<')
        auc[t] <- auc.list$auc
        
        if(sd(as.numeric(res$pred)) != 0){
          prc.list <- pr.curve(scores.class0=res[res$ground.truth.fix==1, ]$pred,
                               scores.class1=res[res$ground.truth.fix==0, ]$pred,
                               curve=TRUE, max.compute=TRUE, min.compute=TRUE)
          prc[t] <- prc.list$auc.integral
        } else{
          prc[t] <- length(res[res$ground.truth.fix==1, ]$pred) / dim(res)[1]
        }
        
      }
      
      
    }
    
    auc.stationary.free <- rbind(auc.stationary.free, auc)
    prc.stationary.free <- rbind(prc.stationary.free, prc)
    
  }
  
  if(Save.Files){
    write.csv(auc.stationary.free, paste(root,'auc_no_smooth.csv', sep = ""))
    write.csv(prc.stationary.free, paste(root,'prc_no_smooth.csv', sep = ""))
  }
  
  
}




# Make top 10%, 20%, ... 100% keys available for the model in CV cohort
if(CV.Top.keys){
  
  if(population){
    sn <- 'adultsAUPRC'
  }else{
    sn <- 'pedsAUPRC'
  }
  
  top.keys <- read_excel('TOP.xlsx', sheet = sn)
  len <- length(top.keys$keys)
  
  auc.topkey <- c()
  prc.topkey <- c()
  
  for(num.feat in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)){
    
    for (f in 1 : nfold) {
      
      input.top.keys <- top.keys$keys[1: round( len * num.feat)]

      test.set <- (1:length(complete.patients))[fold==f]
      severe.I <- setdiff(big.severe.I, test.set)
      mild.I <- setdiff(big.mild.I, test.set)
      
      out <- composite.xvalid(keys=input.top.keys, 
                              train=c(severe.I, mild.I), 
                              test=test.set, 
                              start.t=start.t, 
                              end.t=end.t,  
                              engraftment.aware=engraftment.aware, 
                              use.smoothing=TRUE, 
                              use.laplace=TRUE,
                              use.stationary=TRUE)
      
      test.set.start <- array(NA, length(test.set))
      for (i in 1:length(test.set)) {
        patient <- complete.patients[test.set[i]]
        test.set.start[i] <- df$Severe.aGVHD.onset.day[df$ID == patient]
      }
      test.set.start[is.na(test.set.start)] <- Inf
      
      auc <- array(NA, end.t)
      prc <- array(NA, end.t)
      for (t in start.t:end.t) {
        
        res <- data.frame(
          daGOAT.id = complete.patients[test.set],
          daGOAT.death = death.or.not[test.set], 
          daGOAT.death.date = death.date[test.set],
          ground.truth = as.factor(out[[1]]),
          Severe.aGVHD.time = test.set.start,
          pred = out[[2]][t + offset.t, ]
        )
        
        res <- res[! (res$ground.truth==1 & res$Severe.aGVHD.time <= t), ]
        res <- res[! (res$daGOAT.death==1 & res$daGOAT.death.date <= t), ]
        res <- res[!is.na(res$pred), ]
        
        res$tw <- res$Severe.aGVHD.time - t
        res$ground.truth.fix <- res$ground.truth
        res$ground.truth.fix[ (res$ground.truth==1) & (res$tw > post.window.size)] <- 0 
        
        if(sd(as.numeric(res$ground.truth.fix)) != 0){
          auc.list <- roc(response=res$ground.truth.fix, predictor=res$pred, direction='<')
          auc[t] <- auc.list$auc
          
          prc.list <- pr.curve(scores.class0=res[res$ground.truth.fix==1, ]$pred,
                               scores.class1=res[res$ground.truth.fix==0, ]$pred, 
                               curve=TRUE, max.compute=TRUE, min.compute=TRUE)
          prc[t] <- prc.list$auc.integral
        }
        
      }
      
      auc.topkey <- rbind(auc.topkey, auc)
      prc.topkey <- rbind(prc.topkey, prc)
      
    }
    
  }
  
  if(Sava.Files){
    
    auc.topkey.mean <- rbind( colMeans(auc.topkey[1:nfold, ], na.rm=T), 
                              colMeans(auc.topkey[1:nfold + 1 * nfold, ], na.rm=T),
                              colMeans(auc.topkey[1:nfold + 2 * nfold, ], na.rm=T), 
                              colMeans(auc.topkey[1:nfold + 3 * nfold, ], na.rm=T),
                              colMeans(auc.topkey[1:nfold + 4 * nfold, ], na.rm=T), 
                              colMeans(auc.topkey[1:nfold + 5 * nfold, ], na.rm=T),
                              colMeans(auc.topkey[1:nfold + 6 * nfold, ], na.rm=T), 
                              colMeans(auc.topkey[1:nfold + 7 * nfold, ], na.rm=T),
                              colMeans(auc.topkey[1:nfold + 8 * nfold, ], na.rm=T), 
                              colMeans(auc.topkey[1:nfold + 9 * nfold, ], na.rm=T)
    )
    
    write.csv(auc.topkey, '')
    write.csv(auc.topkey.mean, '')
    
    prc.topkey.mean <- rbind( colMeans(prc.topkey[1:nfold, ], na.rm=T), 
                              colMeans(prc.topkey[1:nfold + 1 * nfold, ], na.rm=T),
                              colMeans(prc.topkey[1:nfold + 2 * nfold, ], na.rm=T), 
                              colMeans(prc.topkey[1:nfold + 3 * nfold, ], na.rm=T),
                              colMeans(prc.topkey[1:nfold + 4 * nfold, ], na.rm=T), 
                              colMeans(prc.topkey[1:nfold + 5 * nfold, ], na.rm=T),
                              colMeans(prc.topkey[1:nfold + 6 * nfold, ], na.rm=T),
                              colMeans(prc.topkey[1:nfold + 7 * nfold, ], na.rm=T),
                              colMeans(prc.topkey[1:nfold + 8 * nfold, ], na.rm=T), 
                              colMeans(prc.topkey[1:nfold + 9 * nfold, ], na.rm=T)
    )
    
    write.csv(prc.topkey, '')
    write.csv(prc.topkey.mean, '')
    
    
  }
  
}




# Holdout Validation
if(HO){
  
  
  # Load and wrangle data.
  encoded.file <- paste(root,'encoded.csv', sep = "")
  df <- read.csv(encoded.file)
  df <- df[, c(-2, -3)]
  df <- colApply(df, factor.col)
  

  complete.patients <- df$ID
  death.or.not <- df$Death
  death.date <- df$Last.followup.time.for.death
  fold <- c(df$fold)
  
  complete.severe <- df$ID[df$Severe.aGVHD.within.100.days == 1]
  complete.mild <- df$ID[df$Severe.aGVHD.within.100.days == 0]
  
  
  data.density <- c()
  complete.patient.data <- list()
  for (i in 1:length(complete.patients)) {
    
    patient <- complete.patients[i]
    file.name <- paste(root, "processed/pt_", patient, ".csv", sep = "")
    patient.data <- read.csv(file.name, fileEncoding = "gbk")
    patient.data <- patient.data[1:day, 2:ncol(patient.data)]
    
    # time-limited sample-and-hold
    for (w in 1:3) {   
      hasValue <- list()
      for (t in offset.t:day) {
        hasValue[[t]] <- !is.na(patient.data[t, ])
      }
      
      for (t in day:offset.t) {
        patient.data[t, !hasValue[[t]]] <- patient.data[t - 1, !hasValue[[t]]]
      }
    }
    
    complete.patient.data[[i]] <- patient.data
    
    
    j <- (1:nrow(df))[df$ID == patient] 
    
    # ignore data measured after severe aGVHD onset
    if (df$Severe.aGVHD.within.100.days[j] == 1) { 
      complete.patient.data[[i]][(df$Severe.aGVHD.onset.day[j]):end.t + offset.t, ] <- NA
      cat(".")
    }
    
    # ignore data measured after death date
    if (df$Death[j] == 1) { 
      complete.patient.data[[i]][(df$Last.followup.time.for.death[j]):end.t + offset.t, ] <- NA
      cat(".")
    }
    
    # ignore data measured before neutrophil engraftment if run in the 'engraftment-aware' mode
    if (engraftment.aware) { 
      j <- (1:nrow(df))[df$ID == patient]
      complete.patient.data[[i]][1:(df$Neutrophil.engraftment.days[j] + 14), ] <- NA
      cat(".")
    }
    
    complete.patient.data[[i]] <- complete.patient.data[[i]][1:day, ]
    data.density[i] <- sum(!is.na(complete.patient.data[[i]][, keys]))

    
  }
  
  
  big.severe.I <- c()
  for (i in 1:length(complete.severe)) {
    big.severe.I <- c(big.severe.I, (1:length(complete.patients))[complete.patients == complete.severe[i]])
  }
  big.severe.I <- big.severe.I[!is.na(big.severe.I)]
  big.mild.I <- setdiff(1:length(complete.patients), big.severe.I)
  
  
  out.pt <- list()
  for (t in start.t:end.t) {
    out1 <- (1:length(complete.patients))[ (df$Severe.aGVHD.within.100.days == 1) & (df$Severe.aGVHD.onset.day <= t) ]
    out2 <- (1:length(complete.patients))[ (df$Death == 1) & (df$Last.followup.time.for.death <= t) ]
    out.pt[[t]] <- union(out1, out2)
  }
  
  
  
  
  ## Model -- daGOAT
  auc.daGOAT <- c()
  prc.daGOAT <- c()
  
  test.set <- (1:length(complete.patients))[as.Date(df$HSCT.time) >= "2020/12/1"]
  severe.I <- setdiff(big.severe.I, test.set)
  mild.I <- setdiff(big.mild.I, test.set)
  
  out <- composite.xvalid(keys=keys, 
                          train=c(severe.I, mild.I), 
                          test=test.set, 
                          start.t=start.t, 
                          end.t=end.t,  
                          engraftment.aware=engraftment.aware, 
                          use.smoothing=TRUE, 
                          use.laplace=TRUE,
                          use.stationary=TRUE)
  
  test.set.start <- array(NA, length(test.set))
  for (i in 1:length(test.set)) {
    patient <- complete.patients[test.set[i]]
    test.set.start[i] <- df$Severe.aGVHD.onset.day[df$ID == patient]
  }
  test.set.start[is.na(test.set.start)] <- Inf
  
  auc <- array(NA, end.t)
  prc <- array(NA, end.t)
  
  for (t in start.t:end.t) {
    
    res <- data.frame(
      daGOAT.id = complete.patients[test.set],
      daGOAT.death = death.or.not[test.set], 
      daGOAT.death.date = death.date[test.set],
      ground.truth = as.factor(out[[1]]),
      Severe.aGVHD.time = test.set.start,
      pred = out[[2]][t + offset.t, ]
    )
    
    res <- res[! (res$ground.truth==1 & res$Severe.aGVHD.time <= t), ]
    res <- res[! (res$daGOAT.death==1 & res$daGOAT.death.date <= t), ]
    res <- res[!is.na(res$pred), ]
    
    res$tw <- res$Severe.aGVHD.time - t
    res$ground.truth.fix <- res$ground.truth
    res$ground.truth.fix[ (res$ground.truth==1) & (res$tw > post.window.size)] <- 0 
    
    if(sd(as.numeric(res$ground.truth.fix)) != 0){
      auc.list <- roc(response=res$ground.truth.fix, predictor=res$pred, direction='<')
      auc[t] <- auc.list$auc
      
      if(sd(as.numeric(res$pred)) != 0){
        prc.list <- pr.curve(scores.class0=res[res$ground.truth.fix==1, ]$pred,
                             scores.class1=res[res$ground.truth.fix==0, ]$pred,
                             curve=TRUE, max.compute=TRUE, min.compute=TRUE)
        prc[t] <- prc.list$auc.integral
      } else{
        prc[t] <- length(res[res$ground.truth.fix==1, ]$pred) / dim(res)[1]
      }
    }
    
    if(t==peak.t){
      peak.daGOAT <- res
    }
    
  }
  
  auc.daGOAT <- rbind(auc.daGOAT, auc)
  prc.daGOAT <- rbind(prc.daGOAT, prc)
  
  
  
  ## Model -- Stationary Features Naive Bayes 
  auc.NB <- c()
  prc.NB <- c()
  
  test.set.raw <- (1:length(complete.patients))[as.Date(df$HSCT.time) >= "2020/12/1"]
  severe.I.raw <- setdiff(big.severe.I, test.set.raw)
  mild.I.raw <- setdiff(big.mild.I, test.set.raw)
  
  auc <- array(NA, end.t)
  prc <- array(NA, end.t)
  
  for (t in start.t:end.t) {
    
    test.set <- test.set.raw
    test.set <- setdiff(test.set, out.pt[[t]])
    
    mild.I <- mild.I.raw
    mild.I <- setdiff(mild.I, out.pt[[t]])
    
    severe.I <- severe.I.raw
    severe.I <- setdiff(severe.I, out.pt[[t]])
    
    
    train.data <- subset(df, ID %in% complete.patients[c(mild.I, severe.I)])
    test.data <- subset(df, ID %in% complete.patients[test.set])
    
    test.data$tw <- test.data$Severe.aGVHD.onset.day - t
    test.data$Severe.aGVHD.within.100.days.fix <- test.data$Severe.aGVHD.within.100.days
    test.data$Severe.aGVHD.within.100.days.fix[ (test.data$Severe.aGVHD.within.100.days==1) & 
                                                (test.data$tw > Inf) ] <- 0
    
    
    if ( ( sd (as.numeric(train.data$Severe.aGVHD.within.100.days)) != 0 ) &
         ( sum(as.numeric(train.data$Severe.aGVHD.within.100.days) -1 ) > 1 ) & 
         ( sd (as.numeric(test.data$Severe.aGVHD.within.100.days.fix)) != 0 ) 
    ){
      model <- naiveBayes(x=train.data[, c(4:18)], y=train.data$Severe.aGVHD.within.100.days, laplace=10)
      prob <- predict(model, test.data[, c(4:18)], type = c("raw"))
      
      auc.list <- roc(response=test.data$Severe.aGVHD.within.100.days.fix,
                      predictor=log( (prob[, 2] + 0.01 ) / (prob[, 1] + 0.01) ),
                      direction = "<")
      auc[t] <- auc.list$auc
      
      temp <- data.frame(response=auc.list$response, pred=auc.list$predictor, nb.id=test.data$ID)
      prc.list <- pr.curve(scores.class0=temp[temp$response==1, ]$pred,
                           scores.class1=temp[temp$response==0, ]$pred,
                           curve=TRUE, max.compute=TRUE, min.compute=TRUE)
      prc[t] <- prc.list$auc.integral
      
      if(t==peak.t){
        peak.NB <- temp
      }
      
    }
    
  }
  
  auc.NB <- rbind(auc.NB, auc)
  prc.NB <- rbind(prc.NB, prc)
  
  
  
  ## Model -- Stationary Features Random Forest
  auc.RF <- c()
  prc.RF <- c()
  
  test.set.raw <- (1:length(complete.patients))[as.Date(df$HSCT.time) >= "2020/12/1"]
  severe.I.raw <- setdiff(big.severe.I, test.set.raw)
  mild.I.raw <- setdiff(big.mild.I, test.set.raw)
  
  auc <- array(NA, end.t)
  prc <- array(NA, end.t)
  
  for (t in start.t:end.t) {
    
    test.set <- test.set.raw
    test.set <- setdiff(test.set, out.pt[[t]])
    
    mild.I <- mild.I.raw
    mild.I <- setdiff(mild.I, out.pt[[t]])
    
    severe.I <- severe.I.raw
    severe.I <- setdiff(severe.I, out.pt[[t]])
    
    
    train.data <- subset(df, ID %in% complete.patients[c(mild.I, severe.I)])
    test.data <- subset(df, ID %in% complete.patients[test.set])
    
    test.data$tw <- test.data$Severe.aGVHD.onset.day - t
    test.data$Severe.aGVHD.within.100.days.fix <- test.data$Severe.aGVHD.within.100.days
    test.data$Severe.aGVHD.within.100.days.fix[ (test.data$Severe.aGVHD.within.100.days==1) & 
                                                (test.data$tw > Inf)] <- 0
    
    
    if ( ( sd (as.numeric(train.data$Severe.aGVHD.within.100.days)) != 0 ) &
         ( sum(as.numeric(train.data$Severe.aGVHD.within.100.days) - 1) > 1 ) & 
         ( sd (as.numeric(test.data$Severe.aGVHD.within.100.days.fix)) != 0 ) 
    ){
      set.seed(100)
      model <- randomForest(Severe.aGVHD.within.100.days ~ Sex + Age + BMI + 
                              Primary.disease + Conditioning.regimen + Sex.match + ABO.match + 
                              TNC.count + CD34.count + Stem.cell.source + Donor.type + 
                              Year + GVHD.prophylaxis + Antithymocyte.globulin.in.conditioning + 
                              HLA.mismatch,
                            data=train.data,  ntree=ntree,  mtry=mtry, importance=TRUE, proximity = TRUE)
      prob <- predict(model, test.data[, c(4:18)], type = "prob")[, 2]
      
      auc.list <- roc(response=test.data$Severe.aGVHD.within.100.days.fix,
                      predictor=prob,
                      direction = "<")
      auc[t] <- auc.list$auc
      
      temp <- data.frame(response=auc.list$response, pred=auc.list$predictor, rf.id=test.data$ID)
      if(sd(as.numeric(temp$pred)) != 0){
        prc.list <- pr.curve(scores.class0=temp[temp$response==1, ]$pred,
                             scores.class1=temp[temp$response==0, ]$pred,
                             curve=TRUE, max.compute=TRUE, min.compute=TRUE)
        prc[t] <- prc.list$auc.integral
      } else{
        prc[t] <- length(temp[temp$response==1, ]$pred) / dim(temp)[1]
      }
      
      
      if(t==peak.t){
        peak.RF <- temp
      }
      
    }
    
  }
  
  auc.RF <- rbind(auc.RF, auc)
  prc.RF <- rbind(prc.RF, prc)
  
  
  
  
  ## Model -- XGBoost
  dynamic.data <- list()
  for (t in start.t:end.t) {
    temp <- c()
    for (j in 1:length(complete.patients)) {
      temp <- rbind(temp, complete.patient.data[[j]][t + offset.t, keys])
    }
    dynamic.data[[t]] <- temp
  }
  
  
  dynamic.data.tw <- list()
  for (t in start.t : end.t) {
    
    temp <- dynamic.data[[t]]
    
    if(t>1){
      for (u in (t-1) : max(1, t-pre.window.size+1) ){
        temp <- cbind(temp, dynamic.data[[u]])
      }
    }
    temp <- cbind(temp, df[, c(4:18)])
    temp$num.feat <- dim(temp)[2]
    temp$Severe.aGVHD.onset.day <- df$Severe.aGVHD.onset.day
    temp$Severe.aGVHD.within.100.days <- df$Severe.aGVHD.within.100.days
    temp$Severe.aGVHD.within.100.days.fix <- df$Severe.aGVHD.within.100.days
    temp$tw <- temp$Severe.aGVHD.onset.day - t
    temp$Severe.aGVHD.within.100.days.fix[ (temp$Severe.aGVHD.within.100.days == 1) & 
                                             (temp$tw > post.window.size) ] <- 0
    
    dynamic.data.tw[[t]] <- temp
  }
  

  auc.XGB <- c()
  prc.XGB <- c()
  
  test.set.raw <- (1:length(complete.patients))[as.Date(df$HSCT.time) >= "2020/12/1"]
  severe.I.raw <- setdiff(big.severe.I, test.set.raw)
  mild.I.raw <- setdiff(big.mild.I, test.set.raw)
  
  auc <- array(NA, end.t)
  prc <- array(NA, end.t)
  
  for (t in start.t:end.t) {
    
    test.set <- test.set.raw
    test.set <- setdiff(test.set, out.pt[[t]])
    
    mild.I <- mild.I.raw
    mild.I <- setdiff(mild.I, out.pt[[t]])
    
    severe.I <- severe.I.raw
    severe.I <- setdiff(severe.I, out.pt[[t]])
    
    train.data <- dynamic.data.tw[[t]][c(mild.I, severe.I), ]
    test.data <- dynamic.data.tw[[t]][c(test.set), ]
    
    feat <- dynamic.data.tw[[t]]$num.feat[1]
    
    train.x <- data.matrix(train.data[, c(1:feat)])
    train.y <- data.matrix(train.data$Severe.aGVHD.within.100.days.fix)
    
    test.x <- data.matrix(test.data[, c(1:feat)])
    test.y <- data.matrix(test.data$Severe.aGVHD.within.100.days.fix)
    
    xgbtrain <- xgb.DMatrix(data=train.x, label=train.y)
    xgbtest <- xgb.DMatrix(data=test.x, label=test.y)
    
    watchlist = list(train=xgbtrain, test=xgbtest)
    set.seed(200)
    model <- xgb.train(data=xgbtrain, 
                       watchlist=watchlist, 
                       objective = "binary:logistic",
                       nrounds=10,
                       verbose=FALSE)
    prob <- predict(model, test.x)
    
    res <- data.frame(pred=prob, 
                      ground.truth=test.y, 
                      xgb.id=complete.patients[test.set])
    
    if(sd(as.numeric(res$ground.truth)) != 0){
      
      auc.list <- roc(response=res$ground.truth, predictor=res$pred, direction='<')
      auc[t] <- auc.list$auc
      
      if(sd(as.numeric(res$pred)) != 0){
        prc.list <- pr.curve(scores.class0=res[res$ground.truth==1, ]$pred,
                             scores.class1=res[res$ground.truth==0, ]$pred,
                             curve=TRUE, max.compute=TRUE, min.compute=TRUE)
        prc[t] <- prc.list$auc.integral
      } else{
        prc[t] <- length(res[res$ground.truth==1, ]$pred) / dim(res)[1]
      }
      
    }
    
    if(t==peak.t){
      peak.XGB <- res
    }
    
  }
  
  auc.XGB <- rbind(auc.XGB, auc)
  prc.XGB <- rbind(prc.XGB, prc)
  
  
  
  ## Output
  if(Save.Files){
    
    AUC <- rbind(auc.daGOAT, auc.NB, auc.RF, auc.XGB)
    PRC <- rbind(prc.daGOAT, prc.NB, prc.RF, prc.XGB) 
    
    peak.label.pred <- cbind(peak.daGOAT$daGOAT.id, peak.daGOAT$ground.truth.fix, peak.daGOAT$pred, 
                             peak.NB$nb.id, peak.NB$response, peak.NB$pred, 
                             peak.RF$rf.id, peak.RF$response, peak.RF$pred, 
                             peak.XGB$xgb.id, peak.XGB$ground.truth, peak.XGB$pred)
    colnames(peak.label.pred) <- c('daGOAT.ID', 'daGOAT.label', 'daGOAT.pred', 
                                   'NB.ID', 'NB.label', 'NB.pred', 
                                   'RF.ID', 'RF.label', 'RF.pred', 
                                   'XGB.ID', 'XGB.label', 'XGB.pred')
    
    write.csv(rbind(AUC, PRC), paste(root,'auc_prc_ho.csv', sep = "")) 
    
    write.csv(peak.label.pred, paste(root,'peak_ho.csv', sep = "")) 
    
    
  }
 
  
   
}




# Make top 10%, 20%, ... 100% keys available for the model in holdout cohort
if(HO.Top.Keys){
  
  if(population){
    sn <- 'adultsAUPRC'
  }else{
    sn <- 'pedsAUPRC'
  }
  
  top.keys <- read_excel('TOP.xlsx', sheet = sn)
  len <- length(top.keys$keys)
  
  auc.topkey <- c()
  prc.topkey <- c()
  
  for(num.feat in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)){
    
    input.top.keys <- top.keys$keys[1: round( len * num.feat)]
    
    test.set <- (1:length(complete.patients))[as.Date(df$HSCT.time) >= "2020/12/1"]
    severe.I <- setdiff(big.severe.I, test.set)
    mild.I <- setdiff(big.mild.I, test.set)
    
    out <- composite.xvalid(keys=input.top.keys, 
                            train=c(severe.I, mild.I), 
                            test=test.set, 
                            start.t=start.t, 
                            end.t=end.t,  
                            engraftment.aware=engraftment.aware, 
                            use.smoothing=TRUE, 
                            use.laplace=TRUE,
                            use.stationary=TRUE)
    
    test.set.start <- array(NA, length(test.set))
    for (i in 1:length(test.set)) {
      patient <- complete.patients[test.set[i]]
      test.set.start[i] <- df$Severe.aGVHD.onset.day[df$ID == patient]
    }
    test.set.start[is.na(test.set.start)] <- Inf
    
    auc <- array(NA, end.t)
    prc <- array(NA, end.t)
    for (t in start.t:end.t) {
      
      res <- data.frame(
        daGOAT.id = complete.patients[test.set],
        daGOAT.death = death.or.not[test.set], 
        daGOAT.death.date = death.date[test.set],
        ground.truth = as.factor(out[[1]]),
        Severe.aGVHD.time = test.set.start,
        pred = out[[2]][t + offset.t, ]
      )
      
      res <- res[! (res$ground.truth==1 & res$Severe.aGVHD.time <= t), ]
      res <- res[! (res$daGOAT.death==1 & res$daGOAT.death.date <= t), ]
      res <- res[!is.na(res$pred), ]
      
      res$tw <- res$Severe.aGVHD.time - t
      res$ground.truth.fix <- res$ground.truth
      res$ground.truth.fix[ (res$ground.truth==1) & (res$tw > post.window.size)] <- 0 
      
      if(sd(as.numeric(res$ground.truth.fix)) != 0){
        auc.list <- roc(response=res$ground.truth.fix, predictor=res$pred, direction='<')
        auc[t] <- auc.list$auc
        
        prc.list <- pr.curve(scores.class0=res[res$ground.truth.fix==1, ]$pred,
                             scores.class1=res[res$ground.truth.fix==0, ]$pred, 
                             curve=TRUE, max.compute=TRUE, min.compute=TRUE)
        prc[t] <- prc.list$auc.integral
      }
      
    }
    
    auc.topkey <- rbind(auc.topkey, auc)
    prc.topkey <- rbind(prc.topkey, prc)
      
    
    
  }
  
  
  
  if(Save.Files){
  
    write.csv(auc.topkey, paste(root,'auc_top_keys.csv', sep = ""))
    write.csv(prc.topkey, paste(root,'prc_top_keys.csv', sep = ""))
    
  }
  
  
}



