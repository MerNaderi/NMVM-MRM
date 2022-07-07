
##------------------------------------------------------##
##----- mixture of NMVBS linear regression model ------##
##------------------------------------------------------##
mix.Reg.NMVBS.EM <- function(y, x, betas = NULL, sigma2 = NULL, 
                             lambda = NULL, alpha = NULL, 
                             pI = NULL, g = NULL, Class = NULL, 
                             error = 0.00001, iter.max = 100, 
                             Stp.rule = c("Log.like", "Atiken"), 
                             error.est = F, per = 1, print = T, 
                             fix.sigma = F)
{
  if(isTRUE(print)){
    cat(paste(rep("-", 70), sep = "", collapse = ""), "\n")
    cat('Finite mixture of NMVBS-based linear regression model','\n')
  }
  
  p = ncol(x); n = length(y)
  
  if(is.null(sigma2) || is.null(betas) || is.null(pI) || is.null(lambda)){
    
    class = Class
    if(is.null(g) && is.null(class) ) 
      stop("The model is not specified correctly.\n")
    if(!is.null(class)){
      tt = table(Class)
      if(g == 0 || max(as.numeric(labels(tt)$Class)) != g) 
        g = max(as.numeric(labels(tt)$Class))
    }
    
    if(is.null(class)) class = kmeans(y, g)$cluster
    
    betas = matrix(0, p, g)
    sigma2 = numeric(g)
    for(j in 1:g){
      LM = lm(y[class == j] ~ 1 + x[class == j, 2:p] )
      betas[,j] = LM$coefficients
      sigma2[j] = mean(LM$residuals^2)
    }
    pI = table(class)/n
    lambda = rep(1, g)
    alpha = rep(1, g)
  }
  
  g = length(pI)
  
  
  ##-----------------------------------------------##
  ##---- PDF of mixture of NMVBS regression    ----##
  ##----------------------------------------------##
  BesselK = function(x, kappa, log.val = T){
    F1 = log( besselK(x, kappa, expon.scaled = TRUE) ) - x
    F2 = besselK(x, kappa)
    ifelse(log.val == T, return(F1), return(F2))
  }
  
 fGH = function(x, mu, tau, sigma2, kappa, chi, psi, log.val = F){
    chii = chi + (x-mu)^2 / sigma2
    psii = psi + tau^2 / sigma2
    out1 = BesselK(sqrt(chii * psii), kappa - 0.5, log.val = T) -
      BesselK(sqrt(chi * psi), kappa, log.val = T)
    out2 = (kappa/2 - 1/4) * log(chii /  psii)
    out3 = (x - mu) * tau / sigma2 - 0.5 * log(sigma2 * 2 * pi)
    F1 = out1 + out2 + out3
    ifelse(log.val == T, return(F1), return(exp(F1)))
  }
  
  f.NMVBS = function(x, mu, sigma2, lambda, alpha) {
    0.5 * (fGH(x, mu, lambda, sigma2, 0.5, alpha^(-2), alpha^(-2), log.val = F)+
             fGH(x, mu, lambda, sigma2, -0.5, alpha^(-2), alpha^(-2), log.val = F))
  }
  
  d.NMVBS.mix.reg = function(y, mu, sigma2, lambda, alpha, pI, log = F)
  {
    ## y: a vetor of respond observations
    ## mu[, ] a matrix = X^top * beta
    ## sigma2: a vector of sigma2
    g = length(sigma2)
    PDF = 0
    for(j in 1:g) PDF = PDF + pI[j] * f.NMVBS(y, mu[, j], sigma2[j], lambda[j], alpha[j])
    Out = PDF
    Out[which(Out == 0)] = .Machine$double.xmin
    ifelse(log == T, return(log(Out)), return(Out))
  }
  
  ##------------------------------------------------------##
  ##----------------- Expectations -----------------------##
  ##------------------------------------------------------##
  EUY.NMVBS = function(y, mu, sigma2, lambda, alpha)
  {
    R = function(cc, k, a) BesselK(cc, k+a, log.val = T) - BesselK(cc, k, log.val = T)
    
    psii = alpha^(-2) + (lambda)^2 / sigma2
    
    pp = 1/ (1 + exp(fGH(y, mu, lambda, sigma2, -0.5, alpha^(-2), alpha^(-2), log.val = T) -
                       fGH(y, mu, lambda, sigma2, 0.5, alpha^(-2), alpha^(-2), log.val = T)))
    chiij = alpha^(-2) + (y - mu)^2 / sigma2
    
    cc = sqrt(psii * chiij)
    w1 = exp(0.5 * log(chiij / psii) + R(cc, 0, 1))
    t1 = exp(-0.5 * log(chiij/psii) + R(cc, 0, -1))
    w2 = exp(0.5 * log(chiij / psii) + R(cc, -1, 1))
    t2 = exp(-0.5 * log(chiij/psii) + R(cc, -1, -1))
    
    w = pp * w1 + (1 - pp) * w2
    t = pp * t1 + (1 - pp) * t2
    return(list(w = w, t = t))
  }
  
  start.time = Sys.time()
  
  mu = x %*% betas
  lk = lk.old = sum(d.NMVBS.mix.reg(y, mu, sigma2, lambda, alpha, pI, log = T))
  
  criterio = 1; count = 0
  if(isTRUE(print)){
    cat(paste(rep("-", 70), sep = "", collapse = ""), "\n")
    cat("iter =", count, "\t logli.old=",  lk.old, "\n")
    cat(paste(rep("-", 70), sep = "", collapse = ""), "\n")
  }
  
  repeat
  {
    count = count + 1
    tal = matrix(0, n, g)
    
    for(j in 1:g){
      tal[, j] = pI[j] * f.NMVBS(y, mu[, j], sigma2[j], lambda[j], alpha[j])
    }
    
    for(k in 1:n) if(all(tal[k,] == 0)) tal[k,] = .Machine$double.xmin
    tal = tal/rowSums(tal)
    
    SS = 0
    for (j in 1:g)
    {
      ### E-step: 
      EUTT = EUY.NMVBS(y, mu[, j], sigma2[j], lambda[j], alpha[j])
      w = EUTT$w; t = EUTT$t
      
      sij = w + t - 2
      alpha[j] = sqrt(sum(sij * tal[, j])/sum(tal[, j]))
      
      ### M-step: 
      pI[j] = sum(tal[, j]) / n
      Cent = y - mu[, j]
      lambda[j] = sum(tal[, j] * Cent)/ sum(w * tal[, j])
      
      Sxy = t(x) %*% diag(tal[, j]) %*% (y * t - lambda[j])
      Sxx = t(x) %*% diag(t * tal[, j] ) %*% x
      betas[, j] = solve(Sxx) %*% Sxy
      
      mu[, j] = x %*% betas[, j]
      Cent = y - mu[, j]
      bb = t * Cent^2 + (lambda[j])^2 * w - 2 * lambda[j] * Cent
      sigma2[j] = sum(tal[, j] * bb) / sum(tal[, j])
      if(fix.sigma) SS = sum(tal[, j] * bb) + SS
    }
    if(fix.sigma) sigma2 = rep(SS/n, g)
    
    lk.new = sum(d.NMVBS.mix.reg(y, mu, sigma2, lambda, alpha, pI, log = T))
    
    if(is.nan(lk.new)) {
      lk.new = lk.old
      break
    }
    
    lk = c(lk, lk.new)
    if(Stp.rule == "Log.like") criterio = (lk.new - lk.old)/abs(lk.old)
    else { criterio = Stop.rule(lk) }
    
    diff = lk.new -lk.old
    if(count %% per == 0 || is.na(diff))
    {
      if(isTRUE(print)){
        cat('iter =', count, '\t logli =', lk.new, '\t diff =', 
            diff, Stp.rule, "'s diff =", criterio, '\n')
        cat(paste(rep("-", 60), sep = "", collapse = ""), "\n")
      }
    }
    
    if(criterio < error | count == iter.max) break
    lk.old = lk.new
  }
  
  # End of the estimation process
  lk = lk.new
  end.time = Sys.time()
  time.taken = end.time - start.time
  
  
  m <- g * (p + 3) + (g - 1) # Abeta + Sigma + alpha
  if(fix.sigma) m = m - g + 1
  
  aic = -2 * lk + 2 * m
  bic = -2 * lk + log(n) * m
  edc = -2 * lk + 0.2 * sqrt(n) * m
  aic_c = -2 * lk + 2 * n * m / (n - m - 1)
  abic = -2 * lk + m * log((n + 2) / 24)
  
  Group = apply(tal, 1, which.max)
  
  if(!is.null(Class)){
    true.clus = Class
    km.clus = Group
    tab = table(true.clus, km.clus)
    MCR = 1 - sum(diag(tab))/sum(tab)
    RII = aricode :: clustComp(km.clus, true.clus)
    
    if(error.est == T){
      SS = 0
      for(i in 1:n){
        sp = ssigma = slambda = salpha = numeric(g); sbeta = matrix(0, p, g)
        for (j in 1:g)
        {
          EUTT = EUY.NMVBS(y[i], mu[i, j], sigma2[j], lambda[j], alpha[j])
          w = EUTT$w; t = EUTT$t
          sij = w + t - 2
          Cent = y[i] - mu[i, j]
          
          sp[j] = (tal[i, j]/pI[j]) - (tal[i, g]/pI[g])
          slambda[j] = tal[i, j] * (Cent - w*lambda[j]) / sigma2[j]
          ssigma[j] = - tal[i, j] * (sigma2[j] - t * Cent^2 + 
                                       2 * lambda[j] * Cent - w*lambda[j]^2)/(sigma2[j])^2 / 2
          sbeta[, j] = c(tal[i, j] * (t * Cent - lambda[j]) * x[i, ]) / sigma2[j]
          salpha[j] = tal[i, j] * (sij/alpha[j]^3 - 1/alpha[j])
        }
        SSS = c(sp[1:(g-1)], slambda, salpha, ssigma, c(sbeta))
        SS = SS + SSS %*% t(SSS)
      }
      Error = sqrt(diag(inv.mat(SS)))
      
      obj.out = list( time = time.taken, group = Group, m = m, betas = betas,
                      sigma2 = sigma2, lambda = lambda, alpha = alpha, pI = pI, 
                      loglike = lk, aic = aic, bic = bic, edc = edc, 
                      aic_c = aic_c, abic = abic, iter = count, MCR = MCR, 
                      similarity = t(as.matrix(RII)),  
                      cross_class = tab, Zij = tal, Error = Error, 
                      convergence = criterio < error, crite = criterio)
    }
    
    if(error.est == F){
      obj.out = list( time = time.taken, group = Group, m = m, betas = betas,
                      sigma2 = sigma2, lambda = lambda, alpha = alpha, pI = pI, 
                      loglike = lk, aic = aic, bic = bic, edc = edc, 
                      aic_c = aic_c, abic = abic, iter = count, MCR = MCR, 
                      similarity = t(as.matrix(RII)), 
                      cross_class = tab, Zij = tal, convergence = criterio < error, 
                      crite = criterio)}
  }
  if(is.null(Class)){
    
    if(error.est == T){
      SS = 0
      for(i in 1:n){
        sp = ssigma = slambda = salpha = numeric(g); sbeta = matrix(0, p, g)
        for (j in 1:g)
        {
          EUTT = EUY.NMVBS(y[i], mu[i, j], sigma2[j], lambda[j], alpha[j])
          w = EUTT$w; t = EUTT$t
          sij = w + t - 2
          Cent = y[i] - mu[i, j]
          
          sp[j] = (tal[i, j]/pI[j]) - (tal[i, g]/pI[g])
          slambda[j] = tal[i, j] * (Cent - w*lambda[j]) / sigma2[j]
          ssigma[j] = - tal[i, j] * (sigma2[j] - t * Cent^2 + 
                                       2 * lambda[j] * Cent - w*lambda[j]^2)/(sigma2[j])^2 / 2
          sbeta[, j] = c(tal[i, j] * (t * Cent - lambda[j]) * x[i, ]) / sigma2[j]
          salpha[j] = tal[i, j] * (sij/alpha[j]^3 - 1/alpha[j])
        }
        SSS = c(sp[-g], slambda, salpha, ssigma, c(sbeta))
        SS = SS + SSS %*% t(SSS)
      }
      Error = sqrt(diag(inv.mat(SS)))
      
      obj.out = list( time = time.taken, group = Group, m = m, betas = betas,
                      sigma2 = sigma2, lambda = lambda, alpha = alpha, pI = pI, 
                      loglike = lk, aic = aic, bic = bic, edc = edc, 
                      aic_c = aic_c, abic = abic, iter = count, Error = Error, 
                      convergence = criterio < error, crite = criterio)
    }
    
    if(error.est == F){
      obj.out = list( time = time.taken, group = Group, m = m, betas = betas,
                      sigma2 = sigma2, lambda = lambda, alpha = alpha, pI = pI, 
                      loglike = lk, aic = aic, bic = bic, edc = edc, 
                      aic_c = aic_c, abic = abic, iter = count, 
                      convergence = criterio < error, crite = criterio)
    }
  }
  obj.out
}

