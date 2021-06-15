impute.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
name = paste0("sim", impute.id, ".rds")

if (impute.id == 1){
  beta = c(-1.5,0,1)
  err = 0.5
} else if (impute.id == 2){
  beta = c(-1.5, 0.5,1)
  err = 0.5
} else if (impute.id == 3){
  beta = c(-1.5,1,1)
  err = 0.5
} else if (impute.id == 4){
  beta = c(-1.5,0,1)
  err = 0.75
} else if (impute.id == 5){
  beta = c(-1.5, 0.5,1)
  err = 0.75
} else if (impute.id == 6){
  beta = c(-1.5,1,1)
  err = 0.75
} else if (impute.id == 7){
  beta = c(-1.5,0,1)
  err = 1
} else if (impute.id == 8){
  beta = c(-1.5, 0.5,1)
  err = 1
} else if (impute.id == 9){
  beta = c(-1.5,1,1)
  err = 1
}



library(survey)
library(mice)

## influence function for logistic regression
inf.fun <- function(fit) {
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) / nrow(dm)
  ## influence function
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
  infl
}


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

integer.neyman.w2 = function(n.strata, NS, sample.size, upper){
  nc = max(upper + 1)
  s = 1:nc
  arr = matrix(rep(NS, each = nc)/sqrt(s*(s+1)), nrow = n.strata, byrow = T)
  arr.list = as.list(as.data.frame(t(arr)))
  for(i in 1:length(arr.list)){
    arr.list[[i]][upper[i]:nc] = 0
  }
  arr = do.call(cbind, arr.list)
  rk <- order(arr, na.last=TRUE, decreasing=TRUE)[1:(sample.size - n.strata)]
  re.om.zero = table(rk%/%nc + 1) + 1
  re.zero = rep(0, n.strata)
  re.zero[1:n.strata %in% names(re.om.zero)] = re.om.zero 
  re.zero[which.max(upper - re.zero)] =  re.zero[which.max(upper - re.zero)] + length(which(re.zero == 0))
  re.zero
}

## perform raking and IPW

impute.raking = function(data){
  twophase.w2 <- twophase(id = list(~1, ~1), strata = list(NULL, ~stra), 
                          subset = ~insample, data = data)
  svyest_ipw = svyglm(Y ~ X + Z, family = quasibinomial, design = twophase.w2)
  impmodel <- svyglm(X ~ X_tilde + Z, design = twophase.w2)
  data$imputex <- as.vector(predict(impmodel, newdata=data, type="response", se.fit=FALSE))
  phase1model_imp <- glm(Y ~ imputex + Z, family = binomial, data=data)
  inffun_imp <- inf.fun(phase1model_imp)
  colnames(inffun_imp)<-paste0("if",1:ncol(inffun_imp))
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, ~stra), 
                              subset = ~insample, data = cbind(data, inffun_imp), method="simple")
  calformula <- make.formula(colnames(inffun_imp)) 
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking")
  svyest_rak<-svyglm(Y ~ X + Z, family = quasibinomial, design=cal_twophase_imp)
  list(svyest_ipw, svyest_rak)
}

pps.prob = function(n1 = n1, w1.n = n){
  if(any(n1/sum(n1) * w1.n < 3)){
    prob = n1/sum(n1)
    prob[prob<0.05] = 0.05
    prob[which.max(prob)] = prob[which.max(prob)] - (sum(prob) - 1)
  } else {
    prob = n1/sum(n1)
  }
  prob
}


one.sim<-function(beta, N = 4000, n = 600, err){
  ## generate data
  df <- data.frame(Z = rbinom(N, 1, .5))
  df$X <- rnorm(N, 0, 1)
  modelterm = beta[1] + beta[2]*df$X + beta[3]*df$Z
  p_y <- exp(modelterm) / (1 + exp(modelterm))
  df$Y <- rbinom(N, size = 1, p_y)
  
  U = rnorm(N, 0, err)
  df$X_tilde <- df$X + U
  df$id<-1:nrow(df)
  
  # stratification on X_tilde
  df$zstrat<-cut(df$X_tilde, c(-Inf, quantile(df$X_tilde, probs = c(0.25, 0.75)), Inf))
  df$stra = as.numeric(interaction(df$zstrat, df$Y))
  n1 = xtabs(~stra, df)
  strata = list()
  for (index in 1:length(n1)){
    strata[[index]] = which(df$stra == index)  
  }
  
  ###########################
  ##                       ##
  ## Optimal design IF-IPW ##
  ##                       ##
  ###########################
  infl = inf.fun(glm(Y~X+Z,data=df, family = binomial))[, 2]
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
  IFIPW = impute.raking(df0)
  ipw.IFIPW = IFIPW[[1]]
  rak.IFIPW = IFIPW[[2]]
  
  ###########################
  ##                       ##
  ## Optimal design IF-GR  ##
  ##                       ##
  ###########################
  df1 = df
  suppressWarnings(ini <- mice(df1, maxit = 0))
  ini$predictorMatrix = matrix(0, nrow = nrow(ini$predictorMatrix), ncol = ncol(ini$predictorMatrix),
                               dimnames = list(row.names(ini$predictorMatrix), row.names(ini$predictorMatrix)))
  ini$predictorMatrix["X", c(1,3:4)] = 1
  meth = ini$method
  meth[meth!=""] = ""
  meth["X"] = "norm"
  mip = matrix(FALSE, nrow = nrow(df1), ncol = ncol(df1))
  mip[,2] = TRUE
  
  M = 50
  imp = mice(df1, meth = meth, pred = ini$predictorMatrix, print = F, m = M,
             where = mip)
  inffun_mi = matrix(0, nrow = N, ncol = 3)
  
  for(i in 1:M){
    phase1model.imp <- glm(Y~X+Z, family = binomial, data=complete(imp, i))
    inffun_mi <- inffun_mi+inf.fun(phase1model.imp)
  }
  inffun_mi<-inffun_mi/M
  colnames(inffun_mi)<-paste0("if",1:ncol(inffun_mi))  
  
  infl.imp = inffun_mi[,2]
  
  
  infl2 = residuals(lm(infl~infl.imp))
  
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
  IFGR = impute.raking(df1)
  ipw.IFGR = IFGR[[1]]
  rak.IFGR = IFGR[[2]]
  
  ###############
  ##           ##  
  ##    SRS    ##     
  ##           ##
  ###############
  df2 = df
  s.srs = sample(1:nrow(df2), n)
  df2$insample = (1:nrow(df2)) %in% s.srs
  SRS = impute.raking(df2)
  ipw.SRS = SRS[[1]]
  rak.SRS = SRS[[2]]
  
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
  BSS = impute.raking(df3)
  ipw.BSS = BSS[[1]]
  rak.BSS = BSS[[2]]
  
  ##############
  ##          ##
  ##   PSS    ##
  ##          ##
  ##############
  df4 = df
  s.pss = list()
  pps.num = round(600 * pps.prob(n1 = n1, w1.n = 600))
  for(i in 1:length(strata)){
    s.pss[[i]] = sample(strata[[i]], pps.num[i])
  }
  
  df4$insample = (1:nrow(df4)) %in% unlist(s.pss)
  PSS = impute.raking(df4)
  ipw.PSS = PSS[[1]]
  rak.PSS = PSS[[2]]
  
  
  
  
  coef.ipw = c((coef(ipw.SRS) - beta)[2], (coef(ipw.BSS) - beta)[2], 
               (coef(ipw.PSS) - beta)[2], (coef(ipw.IFIPW) - beta)[2],
               (coef(ipw.IFGR) - beta)[2])
  
  coef.rak = c((coef(rak.SRS) - beta)[2], (coef(rak.BSS) - beta)[2], 
               (coef(rak.PSS) - beta)[2], (coef(rak.IFIPW) - beta)[2],
               (coef(rak.IFGR) - beta)[2])
  
  
  sample.size = cbind(ney, optimal.rak)
  names(coef.ipw) = names(coef.rak) = c("SRS", "BSS", "PSS", "IF-IPW", "IF-GR")
  
  
  list(coef.ipw = coef.ipw, coef.rak = coef.rak,
       sample.size = sample.size)
}



out = replicate(2000, one.sim(beta = beta, err = err))

