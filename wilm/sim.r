library(survey)


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

########################
##                    ##
## influence function ##
##                    ## 
########################

inf.fun <- function(fit) {
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) / nrow(dm)
  ## influence function
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
  infl
}


########################
##                    ##
## Neyman allocation  ##
##                    ## 
########################

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


#################################
##                             ##
## raking by single imputation ##
##                             ##
#################################

impute.raking = function(data){
  twophase.w2 <- twophase(id = list(~1, ~1), strata = list(NULL, ~stra), 
                          subset = ~insample, data = data)
  svyest_ipw = svyglm(relaps~histol+age1+age2+stage1*tumdiam, family=quasibinomial, design=twophase.w2)
  impmodel <- svyglm(histol~instit+age3+study*stage2,family=quasibinomial,design=twophase.w2)
  data$imputex <- as.vector(predict(impmodel,newdata=data,type="response",se.fit=FALSE))
  phase1model_imp <- glm(relaps~imputex+age1+age2+stage1*tumdiam, family=binomial, data=data)
  inffun_imp <- inf.fun(phase1model_imp)
  colnames(inffun_imp)<-paste0("if",1:ncol(inffun_imp))
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, ~stra), 
                              subset = ~insample, data = cbind(data, inffun_imp), method="simple")
  calformula <- make.formula(colnames(inffun_imp)) 
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking", force = T)
  svyest_rak<-svyglm(relaps~histol+age1+age2+stage1*tumdiam, family=quasibinomial, design=cal_twophase_imp)
  list(svyest_ipw, svyest_rak)
}



impute.raking.srs = function(data){
  twophase.w2 <- twophase(id = list(~1, ~1), strata = list(NULL, NULL), 
                          subset = ~insample, data = data)
  svyest_ipw = svyglm(relaps~histol+age1+age2+stage1*tumdiam, family=quasibinomial, design=twophase.w2)
  impmodel <- svyglm(histol~instit+age3+study*stage2,family=quasibinomial,design=twophase.w2)
  data$imputex <- as.vector(predict(impmodel,newdata=data,type="response",se.fit=FALSE))
  phase1model_imp <- glm(relaps~imputex+age1+age2+stage1*tumdiam, family=binomial, data=data)
  inffun_imp <- inf.fun(phase1model_imp)
  colnames(inffun_imp)<-paste0("if",1:ncol(inffun_imp))
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, NULL), 
                              subset = ~insample, data = cbind(data, inffun_imp), method="simple")
  calformula <- make.formula(colnames(inffun_imp)) 
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking", force = T)
  svyest_rak<-svyglm(relaps~histol+age1+age2+stage1*tumdiam, family=quasibinomial, design=cal_twophase_imp)
  list(svyest_ipw, svyest_rak)
}

## simulation
nwts <- read.table("nwts-share.txt", header=TRUE)
nwts$age1 <- with(nwts, pmin(age, 1))
nwts$age2 <- with(nwts, pmax(age, 1))
nwts$age3 <- ifelse(nwts$age > 10, 1, 0)
nwts$stage1 <- ifelse(nwts$stage > 2, 1, 0)
nwts$stage2 <- ifelse(nwts$stage > 3, 1, 0)
nwts$study <- nwts$study - 3


## simulation
## define strata
## define strata 
#nwts$stra = 1 + 4 * nwts$relaps + 2 * nwts$instit + nwts$study
nwts$stra = as.numeric(interaction(nwts$instit, nwts$relaps))

## population strata size
n1 = xtabs(~stra, nwts)

strata = list()
for (index in 1:length(n1)){
  strata[[index]] = which(nwts$stra == index)  
}

n = sum(n1[3:4]*2)


############################
##                        ##
## Optimal design for IPW ##
##                        ##
############################
full <- glm(relaps~histol+age1+age2+stage1*tumdiam, family=binomial, data=nwts)
beta = coef(full)

infl <- inf.fun(full)[, 2]

sd.stra1 = numeric()
for(i in 1:length(strata)){
  sd.stra1 = c(sd.stra1, sd(infl[strata[[i]]]))
}

ney = integer.neyman.w1(n.strata = length(strata), NS = n1 * sd.stra1, sample.size = n, upper = n1)




one.sim = function(){
  
  ## sampling and calibration
  df0 = nwts
  s.ney = list()
  for(i in 1:length(strata)){
    s.ney[[i]] = sample(strata[[i]], ney[i]) 
  }
  df0$insample = (1:nrow(df0)) %in% unlist(s.ney)
  
  IFIPW = impute.raking(df0)
  ipw.IFIPW = IFIPW[[1]]
  rak.IFIPW = IFIPW[[2]]
  ###############################
  ##                           ##
  ## Optimal design for raking ##
  ##                           ##
  ###############################
  df1 = nwts
  
  imp.full <- glm(histol~instit+age3+study*stage2,family=binomial, data=df1)
  inffun_mi<-matrix(0,nrow=nrow(df1),ncol=7)
  M<-50
  for(i in 1:M){
    df1$impx <- rbinom(nrow(df1), 1, as.vector(predict(imp.full, newdata=df1, type="response", se.fit=FALSE)))
    phase1model_out_rak <- glm(relaps~impx+age1+age2+stage1*tumdiam, family = binomial, data=df1)
    inffun_mi <- inffun_mi+inf.fun(phase1model_out_rak)
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
  df2 = nwts
  s.srs = sample(1:nrow(df2), n)
  df2$insample = (1:nrow(df2)) %in% s.srs
  nsrs = xtabs(~df2$stra[df2$insample == 1]) 
  SRS = impute.raking.srs(df2)
  ipw.SRS = SRS[[1]]
  rak.SRS = SRS[[2]]
  
  
  
  
  
  
  
  ##############
  ##          ##
  ##   PSS    ##
  ##          ##
  ##############
  df4 = nwts
  s.pss = list()
  pps.num = round(n * pps.prob(n1 = n1, w1.n = n))
  for(i in 1:length(strata)){
    s.pss[[i]] = sample(strata[[i]], pps.num[i])
  }
  
  df4$insample = (1:nrow(df4)) %in% unlist(s.pss)
  PSS = impute.raking(df4)
  ipw.PSS = PSS[[1]]
  rak.PSS = PSS[[2]]
  
  
  ##############
  ##          ##
  ##   SCC    ##
  ##          ##
  ##############
  df5 = nwts
  s.scc = list()
  scc.num = rep(n1[3:4], 2)
  for(i in 1:length(strata)){
    s.scc[[i]] = sample(strata[[i]], scc.num[i])
  }
  
  df5$insample = (1:nrow(df5)) %in% unlist(s.scc)
  SCC = impute.raking(df5)
  ipw.SCC = SCC[[1]]
  rak.SCC = SCC[[2]]
  
  
  ####################################
  ##                                ##
  ##  phase-two data only with SRS  ##
  ##                                ##
  ####################################
  df6 = df2[df2$insample,]
  ph.two.srs = glm(relaps~histol+age1+age2+stage1*tumdiam, data = df6, family = binomial)
  
  
  
  
  coef.ipw = c((coef(ipw.SRS) - beta)[2], 
               (coef(ipw.SCC) - beta)[2],
               (coef(ipw.PSS) - beta)[2], 
               (coef(ipw.IFIPW) - beta)[2],
               (coef(ipw.IFGR) - beta)[2],
               (coef(ph.two.srs) - beta)[2])
  
  coef.rak = c((coef(rak.SRS) - beta)[2], 
               (coef(rak.SCC) - beta)[2],
               (coef(rak.PSS) - beta)[2], 
               (coef(rak.IFIPW) - beta)[2],
               (coef(rak.IFGR) - beta)[2],
               (coef(ph.two.srs) - beta)[2])
  
  sd.ipw = c(summary(ipw.SRS)$coefficients[2,2],
             summary(ipw.SCC)$coefficients[2,2],
             summary(ipw.PSS)$coefficients[2,2],
             summary(ipw.IFIPW)$coefficients[2,2],
             summary(ipw.IFGR)$coefficients[2,2],
             summary(ph.two.srs)$coefficients[2,2])
  
  sd.rak = c(summary(rak.SRS)$coefficients[2,2],
             summary(rak.SCC)$coefficients[2,2],
             summary(rak.PSS)$coefficients[2,2],
             summary(rak.IFIPW)$coefficients[2,2],
             summary(rak.IFGR)$coefficients[2,2],
             summary(ph.two.srs)$coefficients[2,2])
  
  
  sample.size = cbind(ney, optimal.rak, nsrs)
  names(coef.ipw) = names(coef.rak) = names(sd.ipw) = names(sd.rak) = c("SRS", "SCC", "PSS", "IF-IPW", "IF-GR", "phase.2")
  
  
  
  
  list(coef.ipw = coef.ipw, coef.rak = coef.rak,
       sd.rak = sd.rak, sd.ipw = sd.ipw,
       sample.size = sample.size)
  
}


out = replicate(2000, one.sim())


