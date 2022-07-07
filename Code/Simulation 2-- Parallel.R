
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

nrun = 100
n = 1000
g = 3

P = 2 # number of covariate x

#######------------data generating -------------#######
x1 = runif(n, -1, 5)
x = cbind(1, x1)

beta01 = c(-3, 1)
beta02 = c(3, 1)
beta03 = c(0.5, 1)
beta0 = cbind(beta01, beta02, beta03)
sigma0 = c(0.2, 0.4, 0.6)
pI0 = c(0.3, 0.3, 0.4)
lambda0 = c(2, 3, 4)
nu0 = cbind(0, 0, 0)

l = 0
while(l <= (nrun - 1)){
  l = l + 1
  out.N = r.mix.reg.SMSN(x, beta0, sigma0, lambda0, nu0, pI0, family = "ESN")
  
  Class = out.N$class
  y = out.N$y
  
  Out.N = mix.Reg.norm.EM(y, x, betas = NULL, sigma2 = NULL, g = g, pI = NULL, 
                          Class = Class, error = 0.00001, iter.max = 1000, 
                          Stp.rule = "Log.like", per = 100, print = F)
  betas = Out.N$betas
  sigma2 = Out.N$sigma2
  pI = Out.N$pI
  
  ### VG
  Out.VG = mix.Reg.VG.EM(y, x, betas, sigma2, lambda = rep(0.1, g), kappa = rep(1, g),
                         psi = rep(1, g), pI, g = g, Class = NULL, error = 0.00001, 
                         iter.max = 1000, Stp.rule = "Log.like", 
                         per = 100, print = F)
  betasVG = Out.VG$betas
  sigma2VG = Out.VG$sigma2
  lambdaVG = Out.VG$lambda
  
  ### SL
  Out.SL = mix.Reg.SLap.EM(y, x, betas, sigma2, lambda = rep(0.1, g), 
                           pI, g = g, Class = NULL, error = 0.00001, 
                           iter.max = 1000, Stp.rule = "Log.like", 
                           per = 100, print = F)
  betasSL = Out.SL$betas
  sigma2SL = Out.SL$sigma2
  lambdaSL = Out.SL$lambda
  
  ### GHST
  Out.GHST = mix.Reg.ST.EM(y, x, betas, sigma2, lambda = rep(0.1, g), nu = rep(1, g), 
                           pI, g = g, Class = NULL, error = 0.00001, iter.max = 1000, 
                           Stp.rule = "Log.like", per = 100, print = F)
  betasGHST = Out.GHST$betas
  sigma2GHST = Out.GHST$sigma2
  lambdaGHST = Out.GHST$lambda
  
  ### NMVBS 
  Out.NMVBS = mix.Reg.NMVBS.EM(y, x, betas, sigma2, lambda = rep(0.1, g), alpha = rep(1, g),
                               pI, g = g, Class = NULL, error = 0.00001, iter.max = 1000, 
                               Stp.rule = "Log.like", per = 100, print = F)
  betasNMVBS = Out.NMVBS$betas
  sigma2NMVBS = Out.NMVBS$sigma2
  lambdaNMVBS = Out.NMVBS$lambda
  
  ### NMVL
  Out.NMVL = mix.Reg.NMVL.EM(y, x, betas, sigma2, lambda = rep(0.1, g), alpha = rep(1, g),
                             pI, g = g, Class = NULL, error = 0.00001, iter.max = 1000, 
                             Stp.rule = "Log.like", per = 1, print = F)
  betasNMVL = Out.NMVL$betas
  sigma2NMVL = Out.NMVL$sigma2
  lambdaNMVL = Out.NMVL$lambda
  
  ### NIG 
  Out.NIG = mix.Reg.nig.EM(y, x, betas, sigma2, lambda = rep(0.1, g), chi = rep(1, g),
                           psi = rep(1, g), pI, g = g, Class = NULL, error = 0.00001, 
                           iter.max = 1000, Stp.rule = "Log.like", per = 100, print = F)
  betasNIG = Out.NIG$betas
  sigma2NIG = Out.NIG$sigma2
  lambdaNIG = Out.NIG$lambda
  
  
  MMRE.VG.In = MMRE.NIG.In = MMRE.GHST.In = MMRE.NMVBS.In =
    MMRE.NMVL.In = MMRE.SL.In = numeric(length = 10)
  
  outli = seq(1, 20, length.out = 10)
  
  for(j in 1:10){
    
    y[which.max(y)] = y[which.max(y)] + outli[j]
    ### VG
    Out.VG = mix.Reg.VG.EM(y, x, betas, sigma2, lambda = rep(0.1, g), kappa = rep(1, g),
                           psi = rep(1, g), pI, g = g, Class = NULL, error = 0.00001, 
                           iter.max = 1000, Stp.rule = "Log.like", 
                           per = 100, print = F)
    MMRE.VG.In[j] = mean(abs((Out.VG$betas - betasVG) / betasVG)) +
      mean(abs((Out.VG$sigma2 - sigma2VG)/sigma2VG)) + mean(abs((Out.VG$lambda - lambdaVG)/lambdaVG))
    
    ### SL
    Out.SL = mix.Reg.SLap.EM(y, x, betas, sigma2, lambda = rep(0.1, g), 
                             pI, g = g, Class = NULL, error = 0.00001, 
                             iter.max = 1000, Stp.rule = "Log.like", 
                             per = 100, print = F)
    MMRE.SL.In[j] = mean(abs((Out.SL$betas - betasSL) / betasSL)) +
      mean(abs((Out.SL$sigma2 - sigma2SL)/sigma2SL)) + mean(abs((Out.SL$lambda - lambdaSL)/lambdaSL))
    
    ### GHST
    Out.GHST = mix.Reg.ST.EM(y, x, betas, sigma2, lambda = rep(0.1, g), nu = rep(1, g), 
                             pI, g = g, Class = NULL, error = 0.00001, iter.max = 1000, 
                             Stp.rule = "Log.like", per = 100, print = F)
    MMRE.GHST.In[j] = mean(abs((Out.GHST$betas - betasGHST)/betasGHST)) +
      mean(abs((Out.GHST$sigma2 - sigma2GHST)/sigma2GHST)) + mean(abs((Out.GHST$lambda - lambdaGHST)/lambdaGHST))
    
    ### NMVBS 
    Out.NMVBS = mix.Reg.NMVBS.EM(y, x, betas, sigma2, lambda = rep(0.1, g), alpha = rep(1, g),
                                 pI, g = g, Class = NULL, error = 0.00001, iter.max = 1000, 
                                 Stp.rule = "Log.like", per = 100, print = F)
    MMRE.NMVBS.In[j] = mean(abs((Out.NMVBS$betas - betasNMVBS) / betasNMVBS)) +
      mean(abs((Out.NMVBS$sigma2 - sigma2NMVBS)/sigma2NMVBS)) + mean(abs((Out.NMVBS$lambda - lambdaNMVBS)/lambdaNMVBS))
    
    ### NMVL
    Out.NMVL = mix.Reg.NMVL.EM(y, x, betas, sigma2, lambda = rep(0.1, g), alpha = rep(1, g),
                               pI, g = g, Class = NULL, error = 0.00001, iter.max = 1000, 
                               Stp.rule = "Log.like", per = 1, print = F)
    MMRE.NMVL.In[j] = mean(abs((Out.NMVL$betas - betasNMVL) / betasNMVL)) +
      mean(abs((Out.NMVL$sigma2 - sigma2NMVL)/sigma2NMVL)) +mean(abs((Out.NMVL$lambda - lambdaNMVL)/lambdaNMVL))
    
    ### NIG 
    Out.NIG = mix.Reg.nig.EM(y, x, betas, sigma2, lambda = rep(0.1, g), chi = rep(1, g),
                             psi = rep(1, g), pI, g = g, Class = NULL, error = 0.00001, 
                             iter.max = 1000, Stp.rule = "Log.like", per = 100, print = F)
    MMRE.NIG.In[j] = mean(abs((Out.NIG$betas - betasNIG) / betasNIG)) +
      mean(abs((Out.NIG$sigma2 - sigma2NIG)/sigma2NIG)) + mean(abs((Out.NIG$lambda - lambdaNIG)/lambdaNIG))
    
    y[which.max(y)] = y[which.max(y)] - outli[j]
  }
  
  write.table(t(as.matrix(MMRE.SL.In)), '/SL.txt', append = TRUE,row.names = F,col.names = F)
  write.table(t(as.matrix(MMRE.GHST.In)), '/GHST.txt', append = TRUE,row.names = F,col.names = F)
  write.table(t(as.matrix(MMRE.NMVBS.In)), '/NMVBS.txt', append = TRUE,row.names = F,col.names = F)
  write.table(t(as.matrix(MMRE.NMVL.In)), '/NMVL.txt', append = TRUE,row.names = F,col.names = F)
  write.table(t(as.matrix(MMRE.VG.In)), '/VG.txt', append = TRUE,row.names = F,col.names = F)
  write.table(t(as.matrix(MMRE.NIG.In)), '/NIG.txt', append = TRUE,row.names = F,col.names = F)
  print(l)
}




MMRE.GHST = colMeans(GHST, na.rm = T)
MMRE.NMVBS = colMeans(NMVBS, na.rm = T)
MMRE.NMVL = colMeans(NMVL, na.rm = T)
MMRE.SL =  colMeans(SL, na.rm = T)
MMRE.VG = colMeans(VG, na.rm = T)
MMRE.NIG = colMeans(NIG, na.rm = T)

outli = seq(1, 20, length.out = 10)
plot(c(0,20), c(0, 1.4), type = "n", xlab = expression(upsilon),
     ylab = "MMRE", yaxt = "n", yaxs = "i", bty = "n",
     main = "")
axis(2, at = seq(0, 1.4, by = 0.2), labels = c(NA, seq(0.2, 1.4, by = 0.2)), pos = 0, las = 2)

lines(outli, MMRE.GHST, type='l')
lines(outli, MMRE.NMVBS, lty=2, col=2)
lines(outli, MMRE.NMVL, lty=3, col=3)
lines(outli, MMRE.SL, lty=4, col=4)
lines(outli, MMRE.NIG, lty=5, col=5)
lines(outli, MMRE.VG, lty=6, col=6)


legend(1, 1.4, legend=c("GHST", "NMVBS", "NMVL", "SL", "NIG", "VG"),
       lty = 1:6, col = 1:6, bty ='n', cex = 1)

