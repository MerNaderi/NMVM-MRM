
rm(list = ls(all=TRUE))

WD.PATH = paste(getwd(),"/Functions", sep = "")

source(paste(WD.PATH, '/Additional.r', sep = ""))
source(paste(WD.PATH, '/normal-mix.r', sep = ""))
source(paste(WD.PATH, '/SL-mix.r', sep = ""))
source(paste(WD.PATH, '/GHST-mix.r', sep = ""))
source(paste(WD.PATH, '/NIG-mix.r', sep = ""))
source(paste(WD.PATH, '/NMVBS-mix.r', sep = ""))
source(paste(WD.PATH, '/NMVL-mix.r', sep = ""))
source(paste(WD.PATH, '/VG-mix.r', sep = ""))

nrun = 200
n = 500
g = 3

P = 2 # number of covariate x

#######------------data generating -------------#######
x1 = runif(n, -1, 5)
xx = cbind(1, x1)

beta01 = c(-3, 1) #for Conc    c(-1, 3)
beta02 = c(3, 1) #for Conc    c(3, -1)
beta03 = c(0, 1) #for Conc    c(3, 1)
beta0 = cbind(beta01, beta02, beta03)
sigma0 = c(0.2, 0.4, 0.6) #for Conc    c(0.6, 0.9, 0.3)
pI0 = c(0.3, 0.3, 0.4)
lambda0 = c(2, 3, 4)
nu0 = cbind(2, 3, 5)

l = 0
while(l <= (nrun - 1)){
  l = l + 1
  out.N = r.mix.reg.SMSN(xx, beta0, sigma0, lambda0, nu0, pI0, family = "ESN")
  
  #Class = out.N$class
  #y = out.N$y
  Class1 = out.N$class
  y1 = out.N$y
  
  
  #### adding noise
  y2 = runif(50, -4, 8)
  x2 = runif(50, -1, 5)
  class2 = sample(1:3, 50, replace = T, prob = pI0)
  
  y = c(y1, y2)
  x = cbind(1, c(x1, x2))
  Class = c(Class1, class2)
  
  Out.N = mix.Reg.norm.EM(y, x, betas = NULL, sigma2 = NULL, g = g, pI = NULL, 
                          Class = Class, error = 0.00001, iter.max = 100, 
                          Stp.rule = "Log.like", per = 100, print = F)
  betas = Out.N$betas
  sigma2 = Out.N$sigma2
  pI = Out.N$pI
  
  ### SL
  Out.SL = mix.Reg.SLap.EM(y, x, betas, sigma2, lambda = rep(0.1, g), 
                           pI, g = g, Class = Class, error = 0.00001, 
                           iter.max = 1000, Stp.rule = "Log.like", 
                           per = 100, print = F)
  
  tab.outSL = t(as.matrix(c(Out.SL$aic, Out.SL$bic, Out.SL$edc, Out.SL$aic_c, Out.SL$abic, 
                            Out.SL$loglike, Out.SL$similarity)))
  ### GHST
  Out.GHST = mix.Reg.ST.EM(y, x, betas, sigma2, lambda = rep(0.1, g), nu = rep(1, 3), 
                           pI, g = g, Class = Class, error = 0.00001, iter.max = 1000, 
                           Stp.rule = "Log.like", per = 100, print = F)
  
  tab.outGHST = t(as.matrix(c(Out.GHST$aic, Out.GHST$bic, Out.GHST$edc, Out.GHST$aic_c, 
                              Out.GHST$abic, Out.GHST$loglike, Out.GHST$similarity)))
  
  ### NMVBS 
  Out.NMVBS = mix.Reg.NMVBS.EM(y, x, betas, sigma2, lambda = rep(0.1, g), alpha = rep(1, 3),
                               pI, g = g, Class = Class, error = 0.00001, iter.max = 1000, 
                               Stp.rule = "Log.like", per = 100, print = F)
  
  tab.outNMVBS = t(as.matrix(c(Out.NMVBS$aic, Out.NMVBS$bic, Out.NMVBS$edc, Out.NMVBS$aic_c, 
                               Out.NMVBS$abic, Out.NMVBS$loglike, Out.NMVBS$similarity)))
  
  ### NMVL
  Out.NMVL = mix.Reg.NMVL.EM(y, x, betas, sigma2, lambda = rep(0.1, g), alpha = rep(1, 3),
                             pI, g = g, Class = Class, error = 0.00001, iter.max = 1000, 
                             Stp.rule = "Log.like", per = 1, print = F)
  
  tab.outNMVL = t(as.matrix(c(Out.NMVL$aic, Out.NMVL$bic, Out.NMVL$edc, Out.NMVL$aic_c, 
                              Out.NMVL$abic, Out.NMVL$loglike, Out.NMVL$similarity)))
  
  ### NIG 
  Out.NIG = mix.Reg.nig.EM(y, x, betas, sigma2, lambda = rep(0.1, g), chi = rep(1, 3),
                           psi = rep(1, 3), pI, g = g, Class = Class, error = 0.00001, 
                           iter.max = 1000, Stp.rule = "Log.like", per = 100, print = F)
  
  tab.outNIG = t(as.matrix(c(Out.NIG$aic, Out.NIG$bic, Out.NIG$edc, Out.NIG$aic_c, Out.NIG$abic, 
                             Out.NIG$loglike, Out.NIG$similarity)))
  
  ### VG
  Out.VG = mix.Reg.VG.EM(y, x, betas, sigma2, lambda = rep(0.1, g), 
                         pI, g = g, Class = Class, error = 0.00001, 
                         iter.max = 1000, Stp.rule = "Log.like", 
                         per = 100, print = F)
  
  tab.outVG = t(as.matrix(c(Out.VG$aic, Out.VG$bic, Out.VG$edc, Out.VG$aic_c, Out.VG$abic, 
                            Out.VG$loglike, Out.VG$similarity)))
  
  write.table(tab.outSL,'/outSL1.txt', append = TRUE,row.names = F,col.names = F)
  write.table(tab.outGHST, '/outGHST1.txt', append = TRUE,row.names = F,col.names = F)
  write.table(tab.outNMVBS, '/outNMVBS1.txt',  append = TRUE,row.names = F,col.names = F)
  write.table(tab.outNMVL, '/outNMVL1.txt', append = TRUE,row.names = F,col.names = F)
  write.table(tab.outNIG, '/outNIG1.txt', append = TRUE,row.names = F,col.names = F)
  write.table(tab.outVG, '/outVG1.txt', append = TRUE,row.names = F,col.names = F)
  print(l)
}



out.N = r.mix.reg.SMSN(x, beta0, sigma0, lambda0, nu0, pI0, family = "ST")
Class = out.N$class
y = out.N$y

plot(x1, y, col = Class)
abline(beta01)
abline(beta02)
abline(beta03)



