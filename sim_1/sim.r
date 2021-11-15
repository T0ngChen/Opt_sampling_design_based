library(survey)

## Exact Neyman allocation (Wright, 2017) https://doi.org/10.1016/j.spl.2017.04.026
integer.neyman.w1 = function(n.strata, NS, sample.size, upper){
  nc = max(upper+1)
  s = 1:nc
  arr = matrix(rep(NS, each = nc)/sqrt(s*(s+1)), nrow = n.strata, byrow = T)
  arr.list = as.list(as.data.frame(t(arr)))
  for(i in 1:length(arr.list)){
    arr.list[[i]][upper[i]:nc] = 0
  }
  arr = do.call(cbind, arr.list)[-1,]
  rk <- order(arr, na.last=TRUE, decreasing=TRUE)[1:(sample.size - 2 * n.strata)]
  re.om.zero = table(rk%/%(nc-1) + 1)
  re.zero = rep(2, n.strata)
  re.zero[1:n.strata %in% names(re.om.zero)] = 2 + re.om.zero
  re.zero
}

## raking by single imputation
impute.raking = function(data){
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, ~stra),
                              subset = ~insample, data = data, method="simple")
  svyest_ipw = svyglm(Y ~ X + Z, design = twophase.w2.imp)
  calformula <- make.formula("X_tilde")
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking")
  svyest_rak <- svyglm(Y ~ X + Z, design=cal_twophase_imp)
  list(svyest_ipw, svyest_rak)
}

## raking by single imputation for SRS
impute.raking.srs = function(data){
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, NULL),
                              subset = ~insample, data = data, method="simple")
  svyest_ipw = svyglm(Y ~ X + Z, design = twophase.w2.imp)
  calformula <- make.formula("X_tilde")
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking")
  svyest_rak <- svyglm(Y ~ X + Z, design=cal_twophase_imp)
  list(svyest_ipw, svyest_rak)
}


## take from https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
generateR <- function(y, rho) {
  x <- rnorm(length(y))
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}


one.sim = function(beta, N = 4000, n = 600, r){
  ## generate data
  df <- data.frame(X = rnorm(N, 0, 1))
  df$Z = rbinom(N, 1, prob=0.5)
  epsilon <- rnorm(N, 0, 1)
  df$Y <- beta[1] + beta[2]*df$X + beta[3]*df$Z + epsilon
  infl = dfbeta(lm(Y ~ X + Z, data=df))[, 2]
  df$X_tilde <- generateR(y = infl, rho = r)
  df$id <- 1:nrow(df)
  
  # stratification on X_tilde
  df$zstrat <- cut(df$X_tilde, c(-Inf, quantile(df$X_tilde, probs = c(0.35, 0.65)), Inf))
  
  df$stra = ifelse(as.numeric(df$zstrat) == 2, 1, 2)
  n1 = xtabs(~stra, df)
  
  ## define strata
  strata = list()
  for (index in 1:length(n1)){
    strata[[index]] = which(df$stra == index)  
  }
  
  ############################
  ##                        ##
  ## Optimal design  IF-IPW ##
  ##                        ##
  ############################
  sd.stra1 = numeric()
  for(i in 1:length(strata)){
    sd.stra1 = c(sd.stra1, sd(infl[strata[[i]]]))
  }
  ney = integer.neyman.w1(n.strata = length(strata), NS = n1 * sd.stra1, sample.size = n, upper = n1)
  
  ## sampling and calibration
  df0 = df
  s.ney = list()
  for(i in 1:length(strata)){
    s.ney[[i]] = sample(strata[[i]], ney[i])
  }
  df0$insample = (1:nrow(df0)) %in% unlist(s.ney)
  ipw.IFIPW = impute.raking(df0)[[1]]
  rak.IFIPW = impute.raking(df0)[[2]]
  
  ###############################
  ##                           ##
  ## Optimal design   IF-GR    ##
  ##                           ##
  ###############################
  df1 = df
  infl2 = residuals(lm(infl~df$X_tilde))
  sd.stra2 = numeric()
  for(kk in 1:length(strata)){
    sd.stra2 = c(sd.stra2, sd(infl2[strata[[kk]]]))
  }
  optimal.rak = integer.neyman.w1(n.strata = length(strata), NS = n1 * sd.stra2, sample.size = n, upper = n1)
  
  ## sampling and calibration
  s.rak = list()
  for(i in 1:length(strata)){
    s.rak[[i]] = sample(strata[[i]], optimal.rak[i])
  }
  df1$insample = (1:nrow(df1)) %in% unlist(s.rak)
  ipw.IFGR = impute.raking(df1)[[1]]
  rak.IFGR = impute.raking(df1)[[2]]
  
  ###############
  ##           ##  
  ##    SRS    ##     
  ##           ##
  ###############
  df2 = df
  s.srs = sample(1:nrow(df2), n)
  df2$insample = (1:nrow(df2)) %in% s.srs
  ipw.SRS = impute.raking.srs(df2)[[1]]
  rak.SRS = impute.raking.srs(df2)[[2]]
  
  ##############
  ##          ##
  ##   BSS    ##
  ##          ##
  ##############
  df3 = df
  s.bss = list()
  for(i in 1:length(strata)){
    s.bss[[i]] = sample(strata[[i]], n/length(strata))
  }
  
  df3$insample = (1:nrow(df3)) %in% unlist(s.bss)
  ipw.BSS = impute.raking(df3)[[1]]
  rak.BSS = impute.raking(df3)[[2]]
  
  ##############
  ##          ##
  ##   PSS    ##
  ##          ##
  ##############
  df4 = df
  s.pss = list()
  for(i in 1:length(strata)){
    s.pss[[i]] = sample(strata[[i]], (n * n1/sum(n1))[i])
  }
  
  df4$insample = (1:nrow(df4)) %in% unlist(s.pss)
  ipw.PSS = impute.raking(df4)[[1]]
  rak.PSS = impute.raking(df4)[[2]]
  
  ####################################
  ##                                ##
  ##  phase-two data only with SRS  ##
  ##                                ##
  ####################################
  df5 = df2[df2$insample,]
  ph.two.srs = lm(Y~X+Z, data = df5)
  
  
  coef.ipw = c((coef(ipw.SRS) - beta)[2], (coef(ipw.BSS) - beta)[2], 
               (coef(ipw.PSS) - beta)[2], (coef(ipw.IFIPW) - beta)[2],
               (coef(ipw.IFGR) - beta)[2], (coef(ph.two.srs) - beta)[2])
  
  coef.rak = c((coef(rak.SRS) - beta)[2], (coef(rak.BSS) - beta)[2], 
               (coef(rak.PSS) - beta)[2], (coef(rak.IFIPW) - beta)[2],
               (coef(rak.IFGR) - beta)[2], (coef(ph.two.srs) - beta)[2])
  
  
  sd.ipw = c(summary(ipw.SRS)$coefficients[2,2],
             summary(ipw.BSS)$coefficients[2,2],
             summary(ipw.PSS)$coefficients[2,2],
             summary(ipw.IFIPW)$coefficients[2,2],
             summary(ipw.IFGR)$coefficients[2,2],
             summary(ph.two.srs)$coefficients[2,2])
  
  sd.rak = c(summary(rak.SRS)$coefficients[2,2],
             summary(rak.BSS)$coefficients[2,2],
             summary(rak.PSS)$coefficients[2,2],
             summary(rak.IFIPW)$coefficients[2,2],
             summary(rak.IFGR)$coefficients[2,2],
             summary(ph.two.srs)$coefficients[2,2])
  
  sample.size = cbind(ney, optimal.rak)
  names(coef.ipw) = names(coef.rak) = names(sd.ipw) = names(sd.rak) = c("SRS", "BSS", "PSS", "IFIPW", "IFGR", "phase.2")
  
  
  list(coef.ipw = coef.ipw, coef.rak = coef.rak,
       sd.rak = sd.rak, sd.ipw = sd.ipw,
       sample.size = sample.size)
}


a = replicate(2000, one.sim(beta = c(1,1,1), r = 0.99))


