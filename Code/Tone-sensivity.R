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
                        error = 0.00001, iter.max = 100, Stp.rule = "Log.like", 
                        per = 100, print = T)
N.without = Out.N$loglike

betas = cbind(Out.N$betas[,1],Out.N$betas[,2])
sigma2 = Out.N$sigma2
pI = Out.N$pI


lambda = c(0, 0)

Out.ST = mix.Reg.ST.EM(y, x, betas, sigma2, lambda = rep(0.1, g), nu = rep(10, g),
                       pI = pI, g = g, Class = NULL, error = 0.00001, 
                       iter.max = 500, Stp.rule = "Log.like", 
                       per = 100, print = T)
ST.without = Out.ST$loglike

Out.nig = mix.Reg.nig.EM(y, x, betas, sigma2, lambda, chi = rep(2, g),
                         psi = rep(2, g), pI, g = 2, Class = NULL, 
                         error = 0.00001, iter.max = 500, Stp.rule = 
                           "Log.like", per = 100, print = T)
nig.without = Out.nig$loglike

Out.NMVBS = mix.Reg.NMVBS.EM(y, x, betas, sigma2, lambda, alpha = rep(2, g), 
                             pI, g, Class = NULL, error = 0.00001, iter.max = 500, 
                             Stp.rule = "Log.like", per = 100, print = T)
NMVBS.without = Out.NMVBS$loglike


#### likelihood distance ##
LD.N = LD.ST = LD.NIG = LD.NMVBS = numeric(length(y))

for(j in 1:length(y)){
  yy = y[-j]; xx = x[-j, ]
  Out.N = mix.Reg.norm.EM(yy, xx, betas = NULL, sigma2 = NULL, g = g, pI = NULL, Class = NULL,
                          error = 0.00001, iter.max = 100, Stp.rule = "Log.like", 
                          per = 100, print = T)
  LD.N[j] = 2 * (N.without - Out.N$loglike)
  
  betas = cbind(Out.N$betas[,1],Out.N$betas[,2])
  sigma2 = Out.N$sigma2
  pI = Out.N$pI
  
  lambda = c(0, 0)
  
  Out.ST = mix.Reg.ST.EM(yy, xx, betas, sigma2, lambda = rep(0.1, g), nu = rep(10, g),
                         pI = pI, g = g, Class = NULL, error = 0.00001, 
                         iter.max = 500, Stp.rule = "Log.like", 
                         per = 100, print = T)
  LD.ST[j] = 2 *(ST.without - Out.ST$loglike)
  
  Out.nig = mix.Reg.nig.EM(yy, xx, betas, sigma2, lambda, chi = rep(2, g),
                           psi = rep(2, g), pI, g = 2, Class = NULL, 
                           error = 0.00001, iter.max = 500, Stp.rule = 
                             "Log.like", per = 100, print = T)
  LD.NIG[j] = 2 *(nig.without - Out.nig$loglike)
  
  Out.NMVBS = mix.Reg.NMVBS.EM(yy, xx, betas, sigma2, lambda, alpha = rep(2, g), 
                               pI, g, Class = NULL, error = 0.00001, iter.max = 500, 
                               Stp.rule = "Log.like", per = 100, print = T)
  LD.NMVBS[j] = 2 *(NMVBS.without - Out.NMVBS$loglike)
  
}

Index = 1:150
V = data.frame(Index, LD = abs(LD.N))
f <- ggplot(V, aes(Index, LD))
f + geom_line( color="grey") + 
  geom_point(shape=21, color="black", fill="#69b3a2", size=2)+
  ggtitle("Normal")


postscript(paste(WD.PATH, '/fig2.eps', sep=''), width=15, height=27)
library(ggpubr)
#theme_set(theme_pubr())
library("gridExtra")

V1 = data.frame(Index, LD = abs(LD.N))
f1 <- ggplot(V1, aes(Index, LD))+ geom_line( color="grey") + 
  geom_point(shape=21, color="black", fill="#69b3a2", size=2)+
  ggtitle("Normal")

V2 = data.frame(Index, LD = abs(LD.NIG))
f2 <- ggplot(V2, aes(Index, LD)) + geom_line( color="grey") + 
  geom_point(shape=21, color="black", fill="#69b3a2", size=2)+
  ggtitle("NIG")

V3 = data.frame(Index, LD = abs(LD.ST))
f3 <- ggplot(V3, aes(Index, LD)) + geom_line( color="grey") + 
  geom_point(shape=21, color="black", fill="#69b3a2", size=2)+
  ggtitle("GHST")


V4 = data.frame(Index, LD = abs(LD.NMVBS)) 
f4 <- ggplot(V4, aes(Index, LD)) + geom_line( color="grey") + 
  geom_point(shape=21, color="black", fill="#69b3a2", size=2)+
  ggtitle("NMVBS")


ggarrange(ggarrange(f1, f2, ncol = 2), # First row with scatter plot
          ggarrange(f3, f4, ncol = 2), # Second row with box and dot plots
          nrow = 2)   #labels = "A"       # Labels of the scatter plot 


dev.off()


