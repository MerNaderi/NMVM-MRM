Stop.rule <- function(log.lik) {
  if (length(log.lik) >= 3) { 
    n = length(log.lik)
    l.new = log.lik[n]
    l.old  = log.lik[(n-1)]
    l.old2  = log.lik[(n-2)]
    ait = (l.new - l.old)/(l.old - l.old2)
    ln.Inf = l.old + (l.new - l.old)/(1 - ait)
    out = ln.Inf - l.new	
    if (!is.na(out) ) out = abs(out)	
    else out = 0
    
  } else out = (log.lik[2] - log.lik[1])/abs(log.lik[1])		
  return( out )
}



mix.sample = function(n, param, family =  c("BS", "Lindly", "gig", "Exp"))
  {
  
  if(family == "BS") {
    alpha = param[1]
    U = rnorm(n)
    out = (alpha * U + sqrt(( alpha * U)^2 + 4))^2 / 4
  }
  
  if(family == "Lindly") {
    alpha = param[1]
    out = VGAM::rlind(n, alpha)
  }
  
  if(family == "gig") {
    kappa = param[1]
    chi = param[2]
    psi = param[3]
    out = GIGrvg::rgig(n, lambda = kappa, chi = chi, psi = psi)
  }
  
  if(family == "Exp") {
    out = rexp(n, 0.5)
  }
  return(out)
}

r.NMV.family = function(n, mu, sigma2, lambda, theta,
                        family = c("normal", "SL", "GHST", 
                                   "VG", "NMVBS", "NMVL",
                                   "NIG", "GH"),
                        SN = FALSE){
  z = rnorm(n, mean = 0, sd = sqrt(sigma2))
  if(family == "SL") W = rexp(n, 0.5)
  if(family == "GHST") W = mix.sample(n, c(-theta/2, theta, 0), family =  "gig")
  if(family == "VG" || 
     family == "NIG"|| 
     family == "GH") W = mix.sample(n, theta, family =  "gig")
  if(family == "NMVBS") W = mix.sample(n, theta, family =  "BS")
  if(family == "NMVL") W = mix.sample(n, theta, family =  "Lindly")
  if(family == "normal") W = 1; lambda = 0
  out = mu + lambda * W + sqrt(W) * z
  if(SN){
    z = sn :: rsn(n, xi = 0, omega = 1, alpha = 1, tau = 0)
    out = mu + lambda * W + sqrt(W) * z
  }
  return(out)
}


r.mix.reg.NMV = function(x, beta, sigma2, lambda, theta, pI,
                         family = c("normal", "SL", "GHST", 
                                    "VG", "NMVBS", "NMVL",
                                    "NIG", "GH"),
                         SN = FALSE){
  n = nrow(x); Zi = y = c()
  
  for(i in 1 : n){
    ni = c(rmultinom(1, 1, pI))
    j = Zi[i] = which.max(ni)
    epsi = r.NMV.family(1, 0, sigma2[j], lambda[j], theta[, j], family = family, SN = SN)
    y[i] = x[i, ] %*% beta[, j] + epsi
  }
  out = list(y = y, class = Zi)
  return(out)
}



r.mix.reg.SMSN = function(x, beta, sigma2, lambda, nu, pI,
                         family = c("ESN", "ST")){
  n = nrow(x); Zi = y = c()
  
  if(family == "ESN"){
    for(i in 1 : n){
      ni = c(rmultinom(1, 1, pI))
      j = Zi[i] = which.max(ni)
      epsi = sn :: rsn(n = 1, xi = 0, omega = sigma2[j], alpha = lambda[j], tau = nu[j])
      y[i] = x[i, ] %*% beta[, j] + epsi
    }
  }
  
  if(family == "ST"){
    for(i in 1 : n){
      ni = c(rmultinom(1, 1, pI))
      j = Zi[i] = which.max(ni)
      epsi = sn :: rst(n = 1, xi = 0, omega = sigma2[j], alpha = lambda[j], nu = nu[j])
      y[i] = x[i, ] %*% beta[, j] + epsi
    }
  }
  
  out = list(y = y, class = Zi)
  return(out)
}


inv.mat = function(M){
  eg = eigen(M)
  val = diag(1/eg$val)
  vec = cbind(eg$vec)
  vec %*% val %*% solve(vec)
}
