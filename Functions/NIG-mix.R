
##------------------------------------------------------##
##------- mixture of NIG linear regression model -------##
##------------------------------------------------------##
mix.Reg.nig.EM = function(y, x, betas = NULL, sigma2 = NULL,
                          lambda = NULL, chi = NULL, psi = NULL, 
                          pI = NULL, g = NULL, Class = NULL,
                          error = 0.00001, iter.max = 100, 
                          Stp.rule = c("Log.like", "Atiken"), 
                          per = 1, print = T, fix.sigma = F)
{
  if(isTRUE(print)){
    cat(paste(rep("-", 70), sep = "", collapse = ""), "\n")
    cat('Finite mixture of NIG-based linear regression model','\n')
  }
  
  # x            <- matrix n * p
  # betas        <- matrix p * G
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
    chi = psi = rep(2, g)
  }
  g = length(pI)
  
  ##------------------------------------------------------##
  ##--------PDF of mixture of NIG regression -------------##
  ##------------------------------------------------------##
  BesselK = function(x, kappa, log.val = T){
    F1 = log( besselK(x, kappa, expon.scaled = TRUE) ) - x
    F2 = besselK(x, kappa)
    ifelse(log.val == T, return(F1), return(F2))
  }
  
  f.nig = function(y, mu, sigma2, lambda, chi, psi, log.val = F){
    ob = ghyp :: NIG(chi = chi, psi = psi, mu = mu, 
                     sigma = sqrt(sigma2), gamma = lambda)
    Val= ghyp :: dghyp(y, object = ob, logvalue = log.val)
    Val
  }
  
  d.nig.mix.reg = function(y, mu, sigma2, lambda, chi, psi, pI, log = F)
  {
    ## y: a vetor of respond observations
    ## mu[, ] a matrix = X^top * beta
    ## sigma2: a vector of sigma2
    g = length(sigma2)
    PDF = 0
    for(j in 1:g){
      cen = y - mu[, j] 
      PDF = PDF + pI[j] * f.nig(cen, 0, sigma2[j], lambda[j], chi[j], psi[j], log.val = F)
    } 
    Out = PDF
    Out[which(Out == 0)] = .Machine$double.xmin
    ifelse(log == T, return(log(Out)), return(Out))
  }
  
  ##------------------------------------------------------##
  ##----------------- Expectations -----------------------##
  ##------------------------------------------------------##
  EUY.nig = function(y, mu, sigma2, lambda, chi, psi)
  {
    R = function(cc, k, a) BesselK(cc, k+a, log.val = T) - BesselK(cc, k, log.val = T)
    
    chii = chi + (y - mu)^2 / sigma2
    psii = psi + lambda^2 / sigma2
    kappa = -1
    cc = sqrt(psii * chii)
    w = exp(0.5 * log(chii / psii) + R(cc, kappa, 1))
    t = exp(-0.5 * log(chii / psii) + R(cc, kappa, -1))
    return(list(w = w, t = t))
  }
  
  start.time = Sys.time()
  
  mu = x %*% betas
  lk = lk.old = sum(d.nig.mix.reg(y, mu, sigma2, lambda, chi, psi, pI, log = T))
  
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
      cen = y - mu[, j]
      tal[, j] = pI[j] * f.nig(cen, 0, sigma2[j], lambda[j], chi[j], psi[j], log.val = F)
    }
    
    for(k in 1:n) if(all(tal[k,] == 0)) tal[k,] = .Machine$double.xmin
    tal = tal/rowSums(tal)
    
    SS = 0
    for (j in 1:g)
    {
      ### E-step: 
      EUTT = EUY.nig(y, mu[, j], sigma2[j], lambda[j], chi[j], psi[j])
      w = EUTT$w; t = EUTT$t
      
      ### M-step: 
      nj = sum(tal[, j])
      pI[j] = nj / n
      Cent = y - mu[, j]
      lambda[j] = sum(tal[, j] * Cent)/ sum(w * tal[, j])
      
      Sxy = t(x) %*% diag(tal[, j]) %*% (y * t - lambda[j])
      Sxx = t(x) %*% diag(t * tal[, j] ) %*% x
      betas[, j] = solve(Sxx) %*% Sxy
      
      mu[, j] = x %*% betas[, j]
      Cent = y - mu[, j]
      bb = t * Cent^2 + (lambda[j])^2 * w - 2 * lambda[j] * Cent
      sigma2[j] = sum( tal[, j] * bb ) / sum(tal[, j])
      if(fix.sigma) SS = sum(tal[, j] * bb) + SS
      
      thetaa = -1 / (1 - sum(w * tal[, j]) * sum(t * tal[, j]) / nj^2 ) 
      chi[j] = thetaa * sum(w * tal[, j])/nj 
      psi[j] = thetaa^2 / chi[j]
    }
    if(fix.sigma) sigma2 = rep(SS/n, g)
    
    lk.new = sum(d.nig.mix.reg(y, mu, sigma2, lambda, chi, psi, pI, log = T))
    
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
    
    obj.out = 
      list( time = time.taken, group = Group, m = m, betas = betas,
            sigma2 = sigma2, lambda = lambda, chi = chi, psi = psi, 
            pI = pI, loglike = lk, aic = aic, bic = bic, edc = edc, 
            aic_c = aic_c, abic = abic, iter = count, MCR = MCR, 
            similarity = t(as.matrix(RII)),  
            cross_class = tab, Zij = tal, convergence = 
              criterio < error, crite = criterio)
  }
  else{obj.out <-
    list( time = time.taken, group = Group, m = m, betas = betas,
          sigma2 = sigma2, lambda = lambda, chi = chi, psi = psi, 
          pI = pI, loglike = lk, aic = aic, bic = bic, edc = edc, 
          aic_c = aic_c, abic = abic, iter = count, convergence = 
            criterio < error, crite = criterio)} 
  obj.out
}

