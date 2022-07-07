rm(list=ls(all=TRUE))

WD.PATH = paste(getwd(),"/Functions", sep = "")

source(paste(WD.PATH, '/Additional.r', sep=""))
n=500
x1 = runif(n, -1, 5)
xx = cbind(1, x1)

beta01 = c(-3, 1)
beta02 = c(3, 1)
beta03 = c(0, 1)
beta0 = cbind(beta01, beta02, beta03)
sigma0 = c(0.2, 0.4, 0.6)
pI0 = c(0.3, 0.3, 0.4)
lambda0 = c(2, 3, 4)
nu0 = cbind(2, 3, 5)

par(mfrow=c(2,2))
out.N = r.mix.reg.SMSN(xx, beta0, sigma0, lambda0, nu0, pI0, family = "ESN")

#Class = out.N$class
#y = out.N$y
Class1 = out.N$class
y1 = out.N$y

y2 = runif(50, -4, 8)
x2 = runif(50, -1, 5)
class2 = sample(1:3, 50, replace = T, prob = pI0)

y = c(y1, y2)
x = cbind(1, c(x1, x2))
Class = c(Class1, class2)

plot(c(x1, x2), y, col = Class, xlab="X",main="ESN", pch = 21, bg = c("red", "green3", "blue")[unclass(Class)], cex.axis = 1.1,
      cex = 0.8, cex.main = 1)
abline(beta01, lty = 1, lwd = 2)
abline(beta02, lty = 1, lwd = 2)
abline(beta03, lty = 1, lwd = 2)





out.N = r.mix.reg.SMSN(xx, beta0, sigma0, lambda0, nu0, pI0, family = "ST")

#Class = out.N$class
#y = out.N$y
Class1 = out.N$class
y1 = out.N$y

y2 = runif(50, -4, 8)
x2 = runif(50, -1, 5)
class2 = sample(1:3, 50, replace = T, prob = pI0)

y = c(y1, y2)
x = cbind(1, c(x1, x2))
Class = c(Class1, class2)

plot(c(x1, x2), y, col = Class, xlab="X",main="ST", pch = 21, bg = c("red", "green3", "blue")[unclass(Class)], cex.axis = 1.1,
      cex = 0.8, cex.main = 1)
abline(beta01, lty = 1, lwd = 2)
abline(beta02, lty = 1, lwd = 2)
abline(beta03, lty = 1, lwd = 2)



beta01 = c(-1, 3)
beta02 = c(3, -1)
beta03 = c(3, 1)
beta0 = cbind(beta01, beta02, beta03)
sigma0 = c(0.6, 0.9, 0.3)
pI0 = c(0.3, 0.3, 0.4)
lambda0 = c(2, 3, 4)
nu0 = cbind(2, 3, 5)

out.N = r.mix.reg.SMSN(xx, beta0, sigma0, lambda0, nu0, pI0, family = "ESN")

#Class = out.N$class
#y = out.N$y
Class1 = out.N$class
y1 = out.N$y

y2 = runif(50, -4, 8)
x2 = runif(50, -1, 5)
class2 = sample(1:3, 50, replace = T, prob = pI0)

y = c(y1, y2)
x = cbind(1, c(x1, x2))
Class = c(Class1, class2)

plot(c(x1, x2), y, col = Class, xlab="X",main="ESN", pch = 21, bg = c("red", "green3", "blue")[unclass(Class)], cex.axis = 1.1,
      cex = 0.8, cex.main = 1)
abline(beta01, lty = 1, lwd = 2)
abline(beta02, lty = 1, lwd = 2)
abline(beta03, lty = 1, lwd = 2)



out.N = r.mix.reg.SMSN(xx, beta0, sigma0, lambda0, nu0, pI0, family = "ST")

#Class = out.N$class
#y = out.N$y
Class1 = out.N$class
y1 = out.N$y

y2 = runif(50, -4, 8)
x2 = runif(50, -1, 5)
class2 = sample(1:3, 50, replace = T, prob = pI0)

y = c(y1, y2)
x = cbind(1, c(x1, x2))
Class = c(Class1, class2)

plot(c(x1, x2), y, col = Class, xlab="X",main="ST", pch = 21, bg = c("red", "green3", "blue")[unclass(Class)], cex.axis = 1.1,
      cex = 0.8, cex.main = 1)
abline(beta01, lty = 1, lwd = 2)
abline(beta02, lty = 1, lwd = 2)
abline(beta03, lty = 1, lwd = 2)



