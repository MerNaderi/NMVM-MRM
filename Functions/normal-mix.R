
##------------------------------------------------------##
##----- mixture of normal linear regression model ------##
##------------------------------------------------------##
mix.Reg.norm.EM <- function(y, x, betas = NULL, sigma2 = NULL, 
                            g = NULL, pI = NULL, Class = NULL,
                            error = 0.00001, iter.max = 100, 
                            Stp.rule = c("Log.like", "Atiken"), 
                            per = 1, print = T, fix.sigma = F)
{
  if(isTRUE(print)){
    cat(paste(rep("-", 70), sep = "", collapse = ""), "\n")
    cat('Finite mixture of normally-based linear regression model','\n')
  }
  
  p <- ncol(x); n = length(y)

  if(is.null(sigma2) || is.null(betas) || is.null(pI) ){
    
    class = Class
    if(is.null(g) && is.null(class) ) 
      stop("The model is not specified correctly.\n")
    if(!is.null(class)){
      tt = table(Class)
      if(g == 0 || max(as.numeric(labels(tt)$Class)) != g) 
        g = max(as.numeric(labels(tt)$Class))
    }
    
    if(is.null(class)) class = ClusterR::Cluster_Medoids(as.matrix(y), g)$clusters 
    #kmeans(y, g)$cluster
    
    betas = matrix(0, p, g)
    sigma2 = numeric(g)
    for(j in 1:g){
      LM = lm(y[class == j] ~ 1 + x[class == j, 2:p] )
      betas[,j] = LM$coefficients
      sigma2[j] = mean(LM$residuals^2)
    }
    pI = table(class)/n
  }
  
  g = length(pI)
  
  
  ##------------------------------------------------------##
  ##-------- PDF of mix-normal distribution --------------##
  ##------------------------------------------------------##
  d.Norm.mix.reg <- function(y, mu, sigma2, pI, log = F)
  {
    ## y: a vetor of respond observations
    ## mu[, ] a matrix = X^top * beta
    ## sigma2: a vector of sigma2
    g = length(sigma2)
    sigma = sqrt(sigma2)
    PDF <- 0
    for(j in 1:g) PDF <- PDF + pI[j] * dnorm(y, mu[, j], sigma[j], log = F)
    Out = PDF
    Out[which(Out == 0)] <- .Machine$double.xmin
    ifelse(log == T, return(log(Out)), return(Out))
  }

  start.time = Sys.time()
  
  mu = x %*% betas
  lk = lk.old = sum(d.Norm.mix.reg(y, mu, sigma2, pI, log = T))
  
  criterio  <- 1
  count  <- 0
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
      tal[, j] = pI[j] * dnorm(y, mu[, j], sqrt(sigma2[j])) 
    }
    
    for(k in 1:n) if(all(tal[k,] == 0)) tal[k,] = .Machine$double.xmin
    tal = tal/rowSums(tal)
    
    SS = 0
    for (j in 1:g)
    {
      ### M-step: 
      pI[j] = sum(tal[, j]) / n 
      Sxy = t(x) %*% diag(tal[, j]) %*% y
      Sxx = t(x) %*% diag(tal[, j]) %*% x
      betas[, j] = inv.mat(Sxx) %*% Sxy
      
      mu[, j] = x %*% betas[, j]
      bb = (y - mu[, j])^2
      sigma2[j] = sum(tal[, j] * bb) / sum(tal[, j])
      if(fix.sigma) SS = sum(tal[, j] * bb) + SS
    }
    
    if(fix.sigma) sigma2 = rep(SS/n, g)
    
    lk.new = sum(d.Norm.mix.reg(y, mu, sigma2, pI, log = T))
    
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
  lk         <- lk.new
  end.time        <- Sys.time()
  time.taken      <- end.time - start.time
  
  
  m <- g * (p + 1) + (g - 1) # Abeta + Sigma + alpha
  if(fix.sigma) m = m - g + 1
  
  aic = -2 * lk + 2 * m
  bic = -2 * lk + log(n) * m
  edc = -2 * lk + 0.2 * sqrt(n) * m
  aic_c = -2 * lk + 2 * n * m / (n - m - 1)
  abic  = -2 * lk + m * log((n + 2) / 24)
  
  Group = apply(tal, 1, which.max)
  
  if(!is.null(Class)){
    true.clus = Class
    km.clus = Group
    tab = table(true.clus, km.clus)
    MCR = 1 - sum(diag(tab))/sum(tab)
    RII = aricode :: clustComp(km.clus, true.clus)
    
    obj.out = list( time = time.taken, group = Group, m = m, betas = betas,
                    sigma2 = sigma2, pI = pI, loglike = lk, aic = aic,
                    bic = bic, edc = edc, aic_c = aic_c, abic = abic,
                    iter = count, MCR = MCR, similarity = t(as.matrix(RII)), 
                    cross_class = tab, Zij = tal, 
                    convergence = criterio < error, crite = criterio)
  }
  else{obj.out <-
    list( time = time.taken, group = Group, m = m, betas = betas,
          sigma2 = sigma2, pI = pI, loglike = lk, aic = aic, 
          bic = bic, edc = edc, aic_c = aic_c, abic = abic,
          iter = count, convergence = criterio < error, 
          crite = criterio)} 
  obj.out
}

