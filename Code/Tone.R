rm (list = ls ())

WD.PATH = paste(getwd(),"/Functions", sep = "")

source(paste(WD.PATH, '/Additional.r', sep=""))
source(paste(WD.PATH, '/normal-mix.r', sep=""))
source(paste(WD.PATH, '/SL-mix.r', sep=""))
source(paste(WD.PATH, '/GHST-mix.r', sep=""))
source(paste(WD.PATH, '/NIG-mix.r', sep=""))
source(paste(WD.PATH, '/NMVBS-mix.r', sep=""))
source(paste(WD.PATH, '/NMVL-mix.r', sep=""))
source(paste(WD.PATH, '/VG-mix.r', sep=""))


library(fpc)
data("tonedata")
y = tonedata[, 2]
x = cbind(1, tonedata[, 1])

outliers = rep(1, length(tonedata[, 2]))
outliers[60] = 2
plot(tonedata[, 1], tonedata[, 2], xlab="X",ylab="Y",pch=16,col=outliers)

g = 2

Out.N = mix.Reg.norm.EM(y, x, betas = NULL, sigma2 = NULL, g = g, pI = NULL, Class = NULL,
                        error = 0.00001, iter.max = 5000, Stp.rule = "Atiken", 
                        per = 100, print = T, fix.sigma = F)


plot(tonedata[, 1], tonedata[, 2], xlab="X",ylab="Y",pch=16,col=outliers)
abline(Out.N$betas[,1])
abline(Out.N$betas[,2])


betas = cbind(Out.N$betas[,1],Out.N$betas[,2])
sigma2 = Out.N$sigma2
pI = Out.N$pI


lambda = c(0, 0)
Out.SLap = mix.Reg.SLap.EM(y, x, betas+0.01, sigma2, lambda = c(1, 2), 
                           pI, g = g, Class = NULL, error = 0.00001, 
                           iter.max = 5000, Stp.rule = "Atiken", 
                           per = 1, print = T, fix.sigma = F)


plot(tonedata[, 1], tonedata[, 2], xlab="X",ylab="Y",pch=16,col=outliers)
abline(Out.SLap$betas[,1])
abline(Out.SLap$betas[,2])

Out.ST = mix.Reg.ST.EM(y, x, betas, sigma2, lambda = rep(0.1, g), nu = rep(10, g),
                       pI = pI, g = g, Class = NULL, error = 0.00001, 
                       iter.max = 5000, Stp.rule = "Atiken", 
                       per = 100, print = T, fix.sigma = F)


plot(tonedata[, 1], tonedata[, 2], xlab="X",ylab="Y",pch=16,col=outliers)
abline(Out.ST$betas[,1])
abline(Out.ST$betas[,2])

Out.ST$lambda
Out.ST$nu
Out.ST$aic
Out.ST$bic

Out.nig = mix.Reg.nig.EM(y, x, betas, sigma2, lambda = rep(0, g), chi = rep(0.2, g),
                         psi = rep(0.2, g), pI, g = 2, Class = NULL, 
                         error = 0.00001, iter.max = 5000, Stp.rule = 
                           "Atiken", per = 100, print = T, fix.sigma = F)


plot(tonedata[, 1], tonedata[, 2], xlab="X",ylab="Y",pch=16,col=outliers)
abline(Out.nig$betas[,1])
abline(Out.nig$betas[,2])

Out.nig$lambda
Out.nig$aic
Out.nig$bic


Out.NMVBS = mix.Reg.NMVBS.EM(y, x, betas, sigma2, lambda, alpha = rep(2, g), 
                             pI, g, Class = NULL, error = 0.00001, iter.max = 5000, 
                             Stp.rule = "Atiken", per = 100, print = T, fix.sigma = F)


plot(tonedata[, 1], tonedata[, 2], xlab="X",ylab="Y",pch=16,col=outliers)
abline(Out.NMVBS$betas[,1])
abline(Out.NMVBS$betas[,2])
Out.NMVBS$lambda
Out.NMVBS$loglike
Out.NMVBS$aic
Out.NMVBS$bic

Out.NMVL = mix.Reg.NMVL.EM(y, x, betas, sigma2, lambda, alpha = rep(2, g), 
                           pI, g, Class = NULL, error = 0.00001, 
                           iter.max = 5000, Stp.rule = "Atiken", per = 100, 
                           print = T, fix.sigma = F)


plot(tonedata[, 1], tonedata[, 2], xlab="X",ylab="Y",pch=16,col=outliers)
abline(Out.NMVL$betas[,1])
abline(Out.NMVL$betas[,2])

Out.NMVL$aic
Out.NMVL$bic

Out.VG = mix.Reg.VG.EM(y, x, betas, sigma2, lambda, kappa = rep(3, g),
                       psi = rep(6, g), pI, g = 2, Class = NULL, 
                       error = 0.00001, iter.max = 39, Stp.rule = 
                         "Atiken", per = 1, print = T, fix.sigma = F)


plot(tonedata[, 1], tonedata[, 2], xlab="X",ylab="Y",pch=16,col=outliers)
abline(Out.VG$betas[,1])
abline(Out.VG$betas[,2])

Out.VG$aic
Out.VG$bic