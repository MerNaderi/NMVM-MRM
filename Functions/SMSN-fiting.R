
dt.lsC <- function(x, loc, sigma2 = 1, shape = 1, nu = 4){
  k1 <- sqrt(nu / 2) * gamma((nu-1)/2) / gamma(nu/2)
  b <-  -sqrt(2/pi)*k1
  delta <- shape / (sqrt(1 + shape^2))
  Delta <- sqrt(sigma2)*delta
  med <- loc+b*Delta
  d <- (x - med)/sqrt(sigma2)
  dens <- 2 * dt(d, df = nu) * pt(sqrt((1+nu)/(d^2 + nu)) * d * shape, 1 + nu)/sqrt(sigma2)
  return(dens)
}

dSNCC <- function(y, mu, sigma2, shape, nu){
  k1 <- nu[1] / nu[2]^(1/2) + 1 - nu[1]
  b <-  -sqrt(2/pi) * k1
  delta <- shape / (sqrt(1 + shape^2))
  Delta <- sqrt(sigma2)*delta
  med<- mu + b * Delta
  dens <- 2 * (nu[1] * dnorm(y, med, sqrt(sigma2/nu[2])) *
                 pnorm(sqrt(nu[2]) * shape * sigma2^(-1/2) * (y - med)) + 
                 (1 - nu[1]) * dnorm(y, med, sqrt(sigma2)) * pnorm(shape * sigma2^(-1/2) * (y - med)))
  return(dens)
}

dSSC <- function(y, mu, sigma2, shape,nu){
  k1<- 2 * nu/(2 * nu-1)
  b<-  -sqrt(2/pi) * k1
  delta <- shape / (sqrt(1 + shape^2))
  Delta <- sqrt(sigma2) * delta
  med<- mu + b * Delta
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    f <- function(u) 2 * nu * u^(nu - 1) * dnorm(y[i], med[i], sqrt(sigma2/u)) *
      pnorm(u^(1/2) * shape * (sigma2^(-1/2)) * (y[i] - med[i]))
    resp[i] <- integrate(f, 0, 0.9999)$value
  }
  return(resp)
}

d.mixedSTC <- function(x, pi1, mu, sigma2, shape, nu){
  g <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j] * dt.lsC(x, mu[,j], sigma2[j], shape[j], nu[j])
  return(dens)
}

d.mixedSNCC <- function(x, pi1, mu, sigma2, shape, nu){
  g <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j] * dSNCC(x, mu[,j], sigma2[j], shape[j], nu[(2*j-1) : (2*j)])
  return(dens)
}

d.mixedSSC <- function(x, pi1, mu, sigma2, shape, nu){
  g <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j] * dSSC(x, mu[, j], sigma2[j], shape[j], nu[j])
  return(dens)
}


smsn.mixReg.numesma <- function(y, 
                                x, 
                                nu, 
                                Abetas = NULL,  
                                sigma2 = NULL, 
                                shape = NULL, 
                                pii = NULL, 
                                g = NULL,
                                family = "Skew.normal", 
                                error = 0.00001, 
                                iter.max = 100){
  
  ################################################################################
  ##                                   Skew-t                                   ##
  ################################################################################
  
  if (family == "Skew.t"){
    k1 <- sqrt(nu/2) * gamma((nu - 1)/2) / gamma(nu/2)
    b <- -sqrt(2/pi) * k1
    
    n <- length(y)
    p <- ncol(x)
    
    delta <- Delta  <- Gama <- Gamaaux1 <- Gamaaux2 <- rep(0,g)
    
    mu <- matrix(0,n,g)
    media <- matrix(0,n,g)
    
    for (k in 1:g){
      delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
      Delta[k] <- sqrt(sigma2[k]) * delta[k]
      Gama[k] <- sigma2[k] - Delta[k]^2
      media[,k] <- x %*% Abetas[,k]
      mu[,k] <- media[, k] + b[k] * Delta[k]
    }
    
    teta <- c(as.vector(Abetas), Delta, Gama, pii, nu)
    
    Abetas.old <- Abetas
    Delta.old <- Delta
    Gama.old <- Gama
    
    criterio <- 1
    count <- 0
    
    lk <- sum(log(d.mixedSTC(y, pii, media, sigma2, shape, nu)))
    
    while((criterio > error) && (count <= iter.max)){
      count <- count + 1
      tal <- matrix(0, n, g)
      S1 <- matrix(0, n, g)
      S2 <- matrix(0, n, g)
      S3 <- matrix(0, n, g)
      for (j in 1:g){
        dj <- ((y - mu[,j]) / sqrt(sigma2[j]))^2
        Mtij2 <- 1/(1 + (Delta[j]^2) * (Gama[j]^(-1)))
        Mtij <- sqrt(Mtij2)
        mutij <- Mtij2 * Delta[j] * (Gama[j]^(-1))*(y - mu[,j]) + b[j]
        A <- (mutij - b[j]) / Mtij
        
        E = (2 * (nu[j])^(nu[j]/2) * gamma((2 + nu[j])/2) *
               ((dj + nu[j] + A^2))^(-(2 + nu[j])/2)) / (gamma(nu[j]/2) * 
                                                         pi * sqrt(sigma2[j]) * 
                                                         dt.lsC(y, media[,j], sigma2[j], shape[j], nu[j]))
        u = ((4*(nu[j])^(nu[j]/2) * gamma((3 + nu[j])/2) * 
                (dj + nu[j])^(-(nu[j] + 3)/2)) / (gamma(nu[j]/2) * sqrt(pi) * 
                                                    sqrt(sigma2[j]) * 
                                                    dt.lsC(y, media[,j], sigma2[j], shape[j], nu[j])) ) *
          pt(sqrt((3 + nu[j])/(dj + nu[j])) * A, 3 + nu[j])
        
        d1 <- dt.lsC(y, media[,j], sigma2[j], shape[j], nu[j])
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
        d2 <- d.mixedSTC(y, pii, media, sigma2, shape, nu)
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
        
        tal[,j] <- d1 * pii[j] / d2
        S1[,j] <- tal[,j] * u
        S2[,j] <- tal[,j] * (mutij * u + Mtij*E)
        S3[,j] <- tal[,j] * (mutij^2  *u + Mtij2 + Mtij * (mutij + b[j]) * E)
        
        S1aux<-as.vector(S1[,j])
        S2aux<-as.vector(S2[,j]/S1[,j])
        
        ### M-step: 
        pii[j] <- (1/n) * sum(tal[,j])
        Delta[j] <- sum(S2[, j] * (y - mu[,j] + b[j] * Delta[j])) / sum(S3[,j])
        Gamaaux1[j] <- sum(S1[,j]*(y - mu[,j] + b[j] * Delta[j])^2 - 2*(y - mu[,j] + b[j] * Delta[j])*
                             Delta[j] * S2[,j] + Delta[j]^2 * S3[,j]) / sum(tal[,j])
        Gamaaux2[j]<-Gamaaux1[j] * sum(tal[,j])
        Abetas[,j] <- solve(t(x) %*% diag(S1aux) %*% x) %*% t(x) %*% diag(S1aux) %*% (y - Delta[j] * S2aux)
        media[,j] <- x %*% Abetas[,j]
        mu[,j] <- media[,j] + b[j]*Delta[j]
      }
      
      Gamafim <- (1/n) * sum(Gamaaux2)
      
      for (jj in 1:g){
        Gama[jj]<- Gamafim
        sigma2[jj] <- Gama[jj] + Delta[jj]^2
        shape[jj] <- ((sigma2[jj]^(-1/2)) * Delta[jj] )/(sqrt(1 - (Delta[jj]^2) * (sigma2[jj]^(-1))))
        }
      
      logvero.ST <- function(nu) -sum(log( d.mixedSTC(y, pii, media, sigma2, shape, nu)))
      nu <- optim(nu, logvero.ST, method = "L-BFGS-B", lower = 2, upper = 50)$par
      
      lk1 <- sum(log( d.mixedSTC(y, pii, media, sigma2, shape, nu) ))
      pii[g] <- 1 - (sum(pii) - pii[g])
      
      zero.pos <- NULL
      zero.pos <- which(pii == 0)
      if(length(zero.pos) != 0){
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }
      
      param <- teta
      teta <- c(as.vector(Abetas), Delta, Gama, pii, nu)
      criterio <- abs(lk1/lk-1)
      
      Abetas.old <- Abetas
      Delta.old <- Delta
      Gama.old <- Gama
      lk <- lk1
    } 
    
    cl <- apply(tal, 1, which.max)
    icl <- 0
    for (j in 1:g) icl <- icl + sum(log(pii[j ]* dt.lsC(y, media[,j], sigma2[j], shape[j], nu[j])))
  }
  
  if (family == "Skew.cn"){
    k1 <- nu[1, ]/ nu[2, ]^(1/2) + 1 - nu[1, ]
    b <- -sqrt(2/pi) * k1
    n <- length(y)
    p <- ncol(x)
    
    delta <- Delta <- Gama <- rep(0,g)
    
    media <- mu <- matrix(0, n, g)
    
    for (k in 1:g){
      delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
      Delta[k] <- sqrt(sigma2[k]) * delta[k]
      Gama[k] <- sigma2[k] - Delta[k]^2
      media[, k] <- x %*% Abetas[,k]
      mu[, k] <- media[,k] + b[k] * Delta[k]
    }
    
    teta <- c(as.vector(Abetas), Delta, Gama, pii, nu)
    Abetas.old<- Abetas
    Delta.old <- Delta
    Gama.old <- Gama
    nu.old<-nu
    
    criterio <- 1
    count <- 0
    
    lk <- sum(log( d.mixedSNCC(y, pii, media, sigma2, shape, as.vector(nu))))
    
    while((criterio > error) && (count <= iter.max)){
      count <- count + 1
      tal <- matrix(0, n, g)
      S1 <- matrix(0, n, g)
      S2 <- matrix(0, n, g)
      S3 <- matrix(0, n, g)
      for (j in 1:g){
        ### E-step: 
        dj <- ((y - mu[,j])/sqrt(sigma2[j]))^2
        Mtij2 <- 1/(1 + (Delta[j]^2) * (Gama[j]^(-1)))
        Mtij <- sqrt(Mtij2)
        mutij <- Mtij2 * Delta[j] * (Gama[j]^(-1)) * (y - mu[,j]) + b[j]
        A <- (mutij - b[j]) / Mtij
        
        u=(2/dSNCC(y, media[,j], sigma2[j], shape[j], nu[, j]))*
          (nu[1, j] * nu[2, j] * dnorm(y, mu[,j], sqrt(sigma2[j] / nu[2, j])) *
             pnorm(sqrt(nu[2,j]) * A) + (1 - nu[1,j]) * dnorm(y, mu[,j], sqrt(sigma2[j])) * pnorm(A))
        E=(2/dSNCC(y, media[,j], sigma2[j], shape[j], nu[, j])) *
          (nu[1] * sqrt(nu[2, j]) * dnorm(y, mu[, j], sqrt(sigma2[j] / nu[2, j])) *
             dnorm(sqrt(nu[2, j]) * A)+(1 - nu[1, j]) * dnorm(y, mu[,j], sqrt(sigma2[j])) * dnorm(A))
        
        d1 <-dSNCC(y, media[,j], sigma2[j], shape[j], nu[, j])
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
        d2 <- d.mixedSNCC(y, pii, media, sigma2, shape, as.vector(nu))
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
        
        tal[,j] <- d1 * pii[j] / d2
        S1[,j] <- tal[,j] * u
        S2[,j] <- tal[,j] * (mutij*u + Mtij * E)
        S3[,j] <- tal[,j] * (mutij^2*u + Mtij2 + Mtij*(mutij + b[j]) * E)
        
        S1aux<-as.vector(S1[,j])
        S2aux<-as.vector(S2[,j]/S1[,j])
        
        ### M-step:
        pii[j] <- (1/n) * sum(tal[,j])
        Delta[j] <- sum(S2[,j] * (y - mu[,j] + b[j] * Delta[j])) / sum(S3[,j])
        Gama[j] <- sum(S1[,j] * (y - mu[,j] + b[j] * Delta[j])^2 - 2*(y - mu[,j] + b[j]*Delta[j]) * 
                         Delta[j] * S2[,j] + Delta[j]^2 * S3[,j]) / sum(tal[,j])
        sigma2[j] <- Gama[j] + Delta[j]^2
        shape[j] <- ((sigma2[j]^(-1/2)) * Delta[j] )/(sqrt(1 - (Delta[j]^2) * (sigma2[j]^(-1))))
        Abetas[,j] <- solve(t(x) %*% diag(S1aux) %*% x) %*% t(x) %*% diag(S1aux) %*% (y - Delta[j] * S2aux)
        media[,j]<- x %*% Abetas[,j]
        mu[,j]<- media[,j] + b[j] * Delta[j]
      }
      
      logvero.SNC <- function(nu) -sum(log( d.mixedSNCC(y, pii, media, sigma2, shape, nu) ))
      NUOPT = try(nlminb(as.vector(nu), logvero.SNC, lower = rep(0.01, 2), upper = rep(0.49, 0.99)))
      if (!('try-error' %in% class(NUOPT))) {
        nu = matrix(NUOPT$par, 2, g) 
      } else {nu = nu}
      lk1 <- sum(log(d.mixedSNCC(y, pii, media, sigma2, shape, as.vector(nu))))
      
      pii[g] <- 1 - (sum(pii) - pii[g])
      
      zero.pos <- NULL
      zero.pos <- which(pii == 0)
      if(length(zero.pos) != 0){
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }
      
      param <- teta
      teta <- c(as.vector(Abetas), Delta, Gama, pii, nu)
      criterio <- abs(lk1/lk-1)
      
      Abetas.old <- Abetas
      Delta.old <- Delta
      Gama.old <- Gama
      nu.old <- nu
      lk<-lk1
      
    }
    
    cl <- apply(tal, 1, which.max)
    icl <- 0
    for (j in 1:g) icl <- icl+sum(log(pii[j]*dSNCC(y, media[,j], sigma2[j], shape[j], nu)))
  }
  
  if (family == "Skew.slash"){
    
    k1<- 2 * nu/(2 * nu - 1)
    b<- -sqrt(2/pi) * k1
    n <- length(y)
    p<- ncol(x)
    
    delta <- Delta  <- Gama<- Gamaaux1 <- Gamaaux2<-rep(0,g)
    
    mu<-matrix(0,n,g)
    media<-matrix(0,n,g)
    
    for (k in 1:g){
      delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
      Delta[k] <- sqrt(sigma2[k]) * delta[k]
      Gama[k] <- sigma2[k] - Delta[k]^2
      media[, k] <- x %*% Abetas[, k]
      mu[, k] <- media[,k] + b[k] * Delta[k]
    }
    
    teta <- c(as.vector(Abetas), Delta, Gama, pii, nu)
    
    Abetas.old<- Abetas
    Delta.old <- Delta
    Gama.old <- Gama
    nu.old<-nu
    criterio <- 1
    count <- 0
    
    lk <- sum(log( d.mixedSSC(y, pii, media, sigma2, shape, nu) ))
    
    while((criterio > error) && (count <= iter.max)){
      count <- count + 1
      tal <- matrix(0, n, g)
      S1 <- matrix(0, n, g)
      S2 <- matrix(0, n, g)
      S3 <- matrix(0, n, g)
      for (j in 1:g){
        dj <- ((y - mu[, j])/sqrt(sigma2[j]))^2
        Mtij2 <- 1/(1 + (Delta[j]^2) * (Gama[j]^(-1)))
        Mtij <- sqrt(Mtij2)
        mutij <- Mtij2 * Delta[j] * (Gama[j]^(-1)) * (y - mu[,j]) + b[j]
        A <- (mutij - b[j]) / Mtij
        u <- vector(mode="numeric", length = n)
        E <- vector(mode="numeric", length = n)
        for(i in 1:n){
          E[i] <- (((2^(nu[j] + 1))*nu[j]*gamma(nu[j] + 1))/(dSSC(y[i], media[i,j], sigma2[j], shape[j], nu[j]) *
                                                               pi * sqrt(sigma2[j]))) * 
            ((dj[i] + A[i]^2)^(-nu[j]-1)) * pgamma(1, nu[j] + 1,(dj[i] + A[i]^2)/2)
          faux <- function(u) u^(nu[j]+0.5) * exp(-u * dj[i]/2) * pnorm(u^(1/2) * A[i])
          aux22 <- integrate(faux,0 , 1)$value
          u[i] <- ((sqrt(2) * nu[j]) / (dSSC(y[i], media[i,j], sigma2[j], shape[j], nu[j]) * 
                                          sqrt(pi) * sqrt(sigma2[j]))) * aux22
        }
        d1 <- dSSC(y, media[,j], sigma2[j], shape[j], nu[j])
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
        d2 <- d.mixedSSC(y, pii, media, sigma2, shape, nu)
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
        
        tal[,j] <- d1 * pii[j] / d2
        S1[,j] <- tal[,j] * u
        S2[,j] <- tal[,j] * (mutij * u + Mtij * E)
        S3[,j] <- tal[,j] * (mutij^2 * u + Mtij2 + Mtij * (mutij + b[j]) * E)
        
        S1aux<-as.vector(S1[,j])
        S2aux<-as.vector(S2[,j]/S1[,j])
        
        ### M-step: 
        pii[j] <- (1/n) * sum(tal[,j])
        Delta[j] <- sum(S2[,j] * (y - mu[,j] + b[j] * Delta[j])) / sum(S3[,j])
        Gamaaux1[j] <- sum(S1[,j]*(y - mu[,j] + b[j] * Delta[j])^2 - 
                             2 * (y - mu[,j] + b[j] * Delta[j]) * Delta[j] * S2[,j] + Delta[j]^2 * S3[,j]) / sum(tal[,j])
        Gamaaux2[j]<-Gamaaux1[j] * sum(tal[,j])
        Abetas[,j] <- solve(t(x) %*% diag(S1aux) %*% x) %*% t(x) %*% diag(S1aux) %*% (y - Delta[j] * S2aux)
        media[,j]<-x %*% Abetas[,j]
        mu[,j]<- media[,j] + b[j] * Delta[j]
      }
      
      Gamafim<-(1/n) * sum(Gamaaux2)
      
      for (jj in 1:g){
        Gama[jj]<- Gamafim
        sigma2[jj] <- Gama[jj] + Delta[jj]^2
        shape[jj] <- ((sigma2[jj]^(-1/2)) * Delta[jj] )/(sqrt(1 - (Delta[jj]^2) * (sigma2[jj]^(-1))))
      }
      
      logvero.SS <- function(nu) -sum(log( d.mixedSSC(y, pii, media, sigma2, shape, nu) ))
      nu <- optim(nu, logvero.SS, method = "L-BFGS-B", lower = 1, upper = 50)$par
      
      lk1 <- sum(log( d.mixedSSC(y, pii, media, sigma2, shape, nu) ))
      
      pii[g] <- 1 - (sum(pii) - pii[g])
      
      zero.pos <- NULL
      zero.pos <- which(pii == 0)
      if(length(zero.pos) != 0){
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }
      
      param <- teta
      teta <- c(mu, Delta, Gama, pii, nu)
      criterio<-abs(lk1/lk-1)
      
      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama
      nu.old<-nu
      lk<-lk1
    }
    
    cl <- apply(tal, 1, which.max)
    icl <- 0
    for (j in 1:g) icl <- icl+sum(log(pii[j]*dSSC(y, media[,j], sigma2[j], shape[j], nu[j])))
  }
  
  
    d <- g * (p + 4) - 1
    if(family == "Skew.cn") d = d + g
    aic <- -2*lk + 2*d
    bic <- -2*lk + log(n)*d
    edc <- -2*lk + 0.2*sqrt(n)*d
    icl <- -2*icl + log(n)*d
    obj.out <- list(Abetas = Abetas, sigma2 = sigma2, shape = shape, pii = pii, 
                    nu = nu, aic = aic, bic = bic, edc = edc, icl = icl,
                    iter = count, n = length(y), group = cl)
  
  class(obj.out) <- family
  obj.out
}

