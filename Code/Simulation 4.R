
rm(list = ls())

WD.PATH = paste(getwd(),"/Functions", sep = "")

source(paste(WD.PATH, '/Additional.r', sep = ""))
source(paste(WD.PATH, '/normal-mix.r', sep = ""))
source(paste(WD.PATH, '/GHST-mix.r', sep = ""))
source(paste(WD.PATH, '/NIG-mix.r', sep = ""))
source(paste(WD.PATH, '/NMVBS-mix.r', sep = ""))
source(paste(WD.PATH, '/SMSN-fiting.R', sep=""))

set.seed(777)

nrun = 100
n = 500
g = 3

#######------------data generating -------------#######
x1 = runif(n, 1, 5)
x2 = runif(n, -1, 1)
x3 = runif(n, 0, 1)
x4 = runif(n, -2, 2)
x = cbind(1, x1, x2, x3, x4)

beta01 = c(-1,-1,-4,-3,-2)
beta02 = c(1, 1, 4, 3, 2)
beta03 = c(-1, 1,-4, 3, 2)
beta0 = cbind(beta01, beta02, beta03)
sigma0 = c(1, 1, 1)
pI0 = c(0.3, 0.5, 0.2)
lambda0 = c(2, 3, 4)
nu0 = cbind(2, 3, 5)


BIC.GHST = ARI.GHST = c()
BIC.NMVBS = ARI.NMVBS = c()
BIC.NIG = ARI.NIG = c()

BIC.ST = ARI.ST = c()
BIC.SCN = ARI.SCN = c()
BIC.SSL = ARI.SSL = c()

l = 0
while (l <= (nrun - 1)) {
  l = l + 1
  out.N = r.mix.reg.NMV(x, beta0, sigma0, lambda0, nu0, pI0, family = "GHST", SN = TRUE)
  
  Class = out.N$class
  y = out.N$y
  
  Out.N = mix.Reg.norm.EM(y, x, betas = NULL, sigma2 = NULL, g = g, pI = NULL, 
                          Class = Class, error = 0.00001, iter.max = 100, 
                          Stp.rule = "Log.like", per = 100, print = F)
  betas = Out.N$betas
  sigma2 = Out.N$sigma2
  pI = Out.N$pI
  
  ### GHST
  Out.GHST = mix.Reg.ST.EM(y, x, betas, sigma2, lambda = rep(0.1, g), nu = rep(1, g), 
                           pI, g = g, Class = Class, error = 0.00001, iter.max = 1000, 
                           Stp.rule = "Log.like", per = 100, print = F)
  BIC.GHST[l] = Out.GHST$bic
  ARI.GHST[l] = Out.GHST$MCR
  
  ### NMVBS 
  Out.NMVBS = mix.Reg.NMVBS.EM(y, x, betas, sigma2, lambda = rep(0.1, g), alpha = rep(1, g),
                               pI, g = g, Class = Class, error = 0.00001, iter.max = 1000, 
                               Stp.rule = "Log.like", per = 100, print = F)
  BIC.NMVBS[l] = Out.NMVBS$bic
  ARI.NMVBS[l] = Out.NMVBS$MCR
  
  ### NIG 
  Out.NIG = mix.Reg.nig.EM(y, x, betas, sigma2, lambda = rep(0.1, g), chi = rep(1, g),
                           psi = rep(1, g), pI, g = g, Class = Class, error = 0.00001, 
                           iter.max = 1000, Stp.rule = "Log.like", per = 100, print = F)
  BIC.NIG[l] = Out.NIG$bic
  ARI.NIG[l] = Out.NIG$MCR
  
  
  #### ST 
  Out.ST = smsn.mixReg.numesma(y, 
                               x, 
                               nu = rep(3, g), 
                               Abetas = betas,  
                               sigma2 = sigma2, 
                               shape = rep(0.1, g), 
                               pii = pI, 
                               g = g, 
                               family = "Skew.t", 
                               error = 0.00001, 
                               iter.max = 1000)
  BIC.ST[l] = Out.ST$bic
  tab = table(Class, Out.ST$group)
  MCR = 1 - sum(diag(tab))/sum(tab)
  ARI.ST[l] = MCR
  
  ### SCN 
  Out.SCN = smsn.mixReg.numesma(y, 
                                x, 
                                nu = matrix(0.1, 2, g), 
                                Abetas = betas,  
                                sigma2 = sigma2, 
                                shape = rep(0.1, g), 
                                pii = pI, 
                                g = g,
                                family = "Skew.cn", 
                                error = 0.00001, 
                                iter.max = 100)
  BIC.SCN[l] = Out.SCN$bic
  tab = table(Class, Out.SCN$group)
  MCR = 1 - sum(diag(tab))/sum(tab)
  ARI.SCN[l] = MCR
  
  ### SSL 
  Out.SSL = smsn.mixReg.numesma(y, 
                                x, 
                                nu = rep(1.3, g), 
                                Abetas = betas,  
                                sigma2 = sigma2, 
                                shape = rep(0.1, g), 
                                pii = pI, 
                                g = g, 
                                family = "Skew.slash", 
                                error = 0.00001, 
                                iter.max = 1000)
  BIC.SSL[l] = Out.SSL$bic
  tab = table(Class, Out.SSL$group)
  MCR = 1 - sum(diag(tab))/sum(tab)
  ARI.SSL[l] = MCR
  
  print(l)
}



