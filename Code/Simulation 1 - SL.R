
rm(list = ls())
WD.PATH = paste(getwd(),"/Functions", sep = "")

source(paste(WD.PATH, '/Additional.r', sep = ""))
source(paste(WD.PATH, '/normal-mix.r', sep = ""))
source(paste(WD.PATH, '/Function/SL-mix.r', sep = ""))


nrun = 500
nsize = c(100, 500, 1000, 2000)
g = 2 

Beta10 = Beta11 = Beta12 = Beta13 = Beta20 = Beta21 = 
  Beta22 = Beta23= sigma1 = sigma22 = lambda1 = lambda2 = pI1 = matrix(NA, nrow = nrun)

Beta10m = Beta11m = Beta12m = Beta13m = Beta20m = Beta21m = 
  Beta22m = Beta23m = sigma1m = sigma2m = lambda1m = lambda2m = 
  pI1m = pI2m = matrix(NA, nrow = nrun, ncol = length(nsize))

for (m in 1:length(nsize)) {
  
  n = nsize[m]
  P = 4 # number of covariate x
  
  #######------------data generating -------------#######
  x1 = runif(n, 1, 5)
  x2 = runif(n, -2, 2)
  x3 = runif(n, 1, 4)
  x = cbind(1, x1, x2, x3)
  
  beta01 = c(0.5, -1, -2, -3)
  beta02 = c(-1, 1, 2, 3)
  beta0 = cbind(beta01, beta02)
  sigma0 = c(1, 2)
  pI0 = c(0.4, 0.6)
  lambda0 = c(2, 3)
  nu0 = cbind(2, 5)
  
  #######---------runing the main function---------#######
  
  Betai = matrix(NA, nrow = nrun, ncol = g * P)
  sigmai = pIi = nui = lambdai = matrix(NA, nrow = nrun, ncol = g)
  
  aici = bici = edci = aic_ci = abici = loglikei = mcri = RIi = ARIi = JIi = c()
  
  mis = 0
  l = 0
  while(l <= (nrun - 1)){
    l = l + 1
    out.N = r.mix.reg.NMV(x, beta0, sigma0, lambda0, nu0, pI0, family = "SL")
    
    Class = out.N$class
    y = out.N$y
    
    Out.N = mix.Reg.norm.EM(y, x, betas = NULL, sigma2 = NULL, g = g, pI = NULL, 
                            Class = Class, error = 0.00001, iter.max = 100, 
                            Stp.rule = "Log.like", per = 100, print = F)
    betas = Out.N$betas
    sigma2 = Out.N$sigma2
    pI = Out.N$pI
    
    Out.SL = mix.Reg.SLap.EM(y, x, betas, sigma2, lambda = rep(0.1, g), 
                             pI, g = g, Class = Class, error = 0.00001, 
                             iter.max = 1000, Stp.rule = "Log.like", 
                             error.est = T, per = 100, print = F)
    tab.out = t(as.matrix(c(as.vector(Out.SL$betas), as.vector(Out.SL$sigma2), 
                            as.vector(Out.SL$lambda), as.vector(Out.SL$pI),
                            Out.SL$aic, Out.SL$bic, Out.SL$edc,
                            Out.SL$aic_c, Out.SL$abic, Out.SL$loglike)))
    
    write.table(tab.out, '/SL/tabout.txt', append = TRUE,row.names = F,col.names = F)
    
    Betai[l, ] = as.vector(Out.SL$betas)
    sigmai[l, ] = as.vector(Out.SL$sigma2)
    lambdai[l, ] = as.vector(Out.SL$lambda)
    pIi[l, ] = as.vector(Out.SL$pI)
    
    aici[l] = Out.SL$aic
    bici[l] = Out.SL$bic
    edci[l] = Out.SL$edc
    aic_ci[l] = Out.SL$aic_c
    abici[l] = Out.SL$abic
    loglikei[l] = Out.SL$loglike
    print(l)
  }
  
  Beta10m[, m] = Betai[, 1]
  Beta11m[, m] = Betai[, 2]
  Beta12m[, m] = Betai[, 3]
  Beta13m[, m] = Betai[, 4]
  Beta20m[, m] = Betai[, 5]
  Beta21m[, m] = Betai[, 6]
  Beta22m[, m] = Betai[, 7]
  Beta23m[, m] = Betai[, 8]
  sigma1m[, m] = sigmai[, 1]
  sigma2m[, m] = sigmai[, 2]
  lambda1m[, m] = lambdai[, 1]
  lambda2m[, m] = lambdai[, 2]
  pI1m[, m] = pIi[, 1]
  
  
  BetaN = matrix(colMeans(Betai, na.rm = T), nrow = P, ncol = g)
  sigmaN = colMeans(sigmai, na.rm = T)
  lambdaN = colMeans(lambdai, na.rm = T)
  pIN = colMeans(pIi, na.rm = T)
  
  Bias.BetaN = matrix(colMeans(abs(t( t(Betai)/c(beta01, beta02) - 1)), na.rm = T), nrow = P, ncol = g) 
  Bias.sigmaN = colMeans(t(abs(t(sigmai)/sigma0 -  1)), na.rm = T)
  Bias.lambdaN = colMeans(t(abs(t(lambdai)/c(lambda0) -  1)), na.rm = T)
  Bias.pIN = colMeans(t(abs(t(pIi)/pI0 -  1)), na.rm = T) 
  
  
  MSE.BetaN = matrix(colMeans((t(t(Betai) - c(beta01,beta02)))^2, na.rm = T), nrow = P, ncol = g)
  MSE.sigmaN = colMeans(t((t(sigmai) - sigma0) ^ 2), na.rm = T)
  MSE.lambdaN = colMeans((t(t(lambdai) - lambda0))^ 2, na.rm = T)
  MSE.pIN = colMeans((t(t(pIi) - pI0))^ 2, na.rm = T)
  
  aicN = mean(aici); bicN = mean(bici)
  edcN = mean(edci); aic_cN = mean(aic_ci)
  abicN = mean(abici); loglikeN = mean(loglikei)
  
  sim.out =
    list(betas = BetaN, Bias.betas = Bias.BetaN,
      MSE.betas = MSE.BetaN, sigma2 = sigmaN,
      Bias.sigma2 = Bias.sigmaN, MSE.sigma2 = MSE.sigmaN,
      lambda = lambdaN, Bias.lambda = Bias.lambdaN, MSE.lambda = MSE.lambdaN,
      pI = pIN, Bias.pI = Bias.pIN, MSE.pI = MSE.pIN,
      aic = aicN, bic = bicN, edc = edcN,
      aic_c = aic_cN, abic = abicN, loglike = loglikeN)
  
  tab1 = round(t(as.matrix(
    c(n, sim.out$aic, sim.out$bic, sim.out$edc,sim.out$aic_c, sim.out$abic,sim.out$loglike))), 4)
  
  write.table(tab1, '/SL/table1.txt',  append = TRUE,row.names = F,col.names = F)
  
  tab2 = round(t(as.matrix(
    c(n, as.vector(sim.out$betas),
      as.vector(sim.out$Bias.betas),
      as.vector(sim.out$MSE.betas)))), 4)
  
  write.table(tab2,'/SL/table2.txt', append = TRUE,row.names = F,col.names = F)
  
  tab2c1 = round(t(as.matrix(
    c(n, as.vector(sim.out$sigma2),
      as.vector(sim.out$Bias.sigma2),
      as.vector(sim.out$MSE.sigma2)))), 4)
  
  write.table(tab2c1, 'SL/table2c1.txt', append = TRUE,row.names = F,col.names = F)
  
  tab2c2 = round(t(as.matrix(
    c(n, as.vector(sim.out$lambda),
      as.vector(sim.out$Bias.lambda),
      as.vector(sim.out$MSE.lambda)))), 4)
  
  write.table(tab2c2, 'SL/table2c2.txt', append = TRUE,row.names = F,col.names = F)
  
  tab2c3 = round(t(as.matrix(
    c(n, as.vector(sim.out$pI),
      as.vector(sim.out$Bias.pI),
      as.vector(sim.out$MSE.pI)))), 4)
  
  write.table(tab2c3, 'SL/table2c3.txt', append = TRUE,row.names = F,col.names = F)
  
}


Beta10 = cbind(Beta10, Beta10m)
Beta11 = cbind(Beta11, Beta11m)
Beta12 = cbind(Beta12, Beta12m)
Beta13 = cbind(Beta13, Beta13m)
Beta20 = cbind(Beta20, Beta20m)
Beta21 = cbind(Beta21, Beta21m)
Beta22 = cbind(Beta22, Beta22m)
Beta23 = cbind(Beta23, Beta23m)
sigma1 = cbind(sigma1, sigma1m)
sigma22 = cbind(sigma22, sigma2m)
lambda1 = cbind(lambda1, lambda1m)
lambda2 = cbind(lambda2, lambda2m)
pI1 = cbind(pI1, pI1m)


write.table(Beta10, '/SL/Beta10.txt', append = TRUE,row.names = F,col.names = F)
write.table(Beta11,'/SL/Beta11.txt', append = TRUE,row.names = F,col.names = F)
write.table(Beta12,'/SL/Beta12.txt', append = TRUE,row.names = F,col.names = F)
write.table(Beta13,'/SL/Beta13.txt', append = TRUE,row.names = F,col.names = F)
write.table(Beta20, '/SL/Beta20.txt', append = TRUE,row.names = F,col.names = F)
write.table(Beta21, '/SL/Beta21.txt', append = TRUE,row.names = F,col.names = F)
write.table(Beta22,'/SL/Beta22.txt', append = TRUE,row.names = F,col.names = F)
write.table(Beta23, '/SL/Beta23.txt', append = TRUE,row.names = F,col.names = F)
write.table(sigma1, '/SL/sigma1.txt', append = TRUE,row.names = F,col.names = F)
write.table(sigma22, '/SL/sigma22.txt', append = TRUE,row.names = F,col.names = F)
write.table(lambda1, '/SL/lambda1.txt', append = TRUE,row.names = F,col.names = F)
write.table(lambda2, '/SL/lambda2.txt', append = TRUE,row.names = F,col.names = F)
write.table(pI1, '/SL/pI1.txt', append = TRUE,row.names = F,col.names = F)
