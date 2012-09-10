# Square root process
#  dv(t) = \kappa (\theta - v(t))dt + \sigma \sqrt{v(t)}dW(t)
#  v(0) = v_0.

library(fOptions)

pintsqr <- function(x, kappa, theta, sigma, u, vu, t, vt) {
	.Call("pintsqr_c", kappa, theta, sigma, u, vu, t, vt, x, PACKAGE="volmodels")
}

qintsqrbridge <- function(p, kappa, theta, sigma, u, vu, tau, vt) {
  # quantile function of \int_u^t v(s)ds given v(u)=vu and v(t)=vt
  uniroot(function(x) {pintsqr(x, kappa, theta, sigma, u, vu, tau, vt) - p}, 
          lower=0.0, upper=1.0)$root
}

qsqrexact <- function(p, kappa, theta, sigma, u, vu, tau) {
  # quantile function of square-root process v(t) given v(u) = vu
  df <- 4*kappa*theta / (sigma^2)
  ncp <- vu * 4*kappa*exp(-kappa*(tau-u))/(sigma^2 *(1- exp(-kappa*(tau-u))))
  k <- sigma^2 * (1-exp(-kappa*(tau-u))) / (4*kappa)
  k * qchisq(p, df, ncp)
}

qhestonexact <- function(p1, p2, p3, r, rho, kappa, theta, sigma, u, vu, Su, tau) {
  vt <- qsqrexact(p1, kappa, theta, sigma, u, vu, tau)
  intv <- mapply(function(p,vt) qintsqrbridge(p, kappa, theta, sigma, u, vu, tau, vt), p2, vt)
  intsqrvdW <- (vt - vu - kappa*theta*(tau-u) + kappa * intv) / sigma
  m <- (r*(tau-u) - 0.5*intv + rho*intsqrvdW)
  s <- sqrt((1-rho^2)*intv)
  Su * exp(m + s*qnorm(p3))
}

rhestonexact <- function(r, rho, kappa, theta, sigma, u, vu, Su, t, n) {
  p <- matrix(runif(n*3),ncol=3)
  qhestonexact(p[,1], p[,2], p[,3], r, rho, kappa, theta, sigma, u, vu, Su, t)
}

rhestonexact.sobol <- function(r, rho, kappa, theta, sigma, u, vu, Su, t, n) {
  p <- runif.sobol(n,3,scrambling=1)
  qhestonexact(p[,1], p[,2], p[,3], r, rho, kappa, theta, sigma, u, vu, Su, t)
}

#hestoncall.mc <- function(r, rho, kappa, theta, sigma, u, vu, Su, t, K, n, m, seed=1) {

hestoncall <- function(lambda, vbar, eta, rho, v0, r, tau, S0, K) {
  PIntegrand <- function(u, lambda, vbar, eta, rho, v0, r, tau, S0, K, j) {
	  F <- S0*exp(r*tau)
	  x <- log(F/K)
	  a <- lambda * vbar
	  
	  if (j == 1) {
	    b <- lambda - rho* eta
	    alpha <- - u^2/2 - u/2 * 1i + 1i * u
	    beta <- lambda - rho * eta - rho * eta * 1i * u
	  } else { # j ==0
	    b <- lambda
	    alpha <- - u^2/2 - u/2 * 1i
	    beta <- lambda - rho * eta * 1i * u
	  }
	
	  gamma <- eta^2/2
	  d <- sqrt(beta^2 - 4*alpha*gamma)
	  rplus <- (beta + d)/(2*gamma)
	  rminus <- (beta - d)/(2*gamma)
	  g <- rminus / rplus
	
	  D <- rminus * (1 - exp(-d*tau))/(1-g*exp(-d*tau))
	  C <- lambda * (rminus * tau - 2/(eta^2) * log( (1-g*exp(-d*tau))/(1-g) ) )
	  
	  top <- exp(C*vbar + D*v0 + 1i*u*x)
	  bottom <- (1i * u)
	  Re(top/bottom)
	}
	
	P <- function(lambda, vbar, eta, rho, v0, r, tau, S0, K, j) {
	  value <- integrate(PIntegrand, lower = 0, upper = Inf, lambda, vbar, eta, rho, v0, r, tau, S0, K, j, subdivisions=2000)$value
	  0.5 + 1/pi * value
	}

  A <- S0*P(lambda, vbar, eta, rho, v0, r, tau, S0, K, 1)
  B <- K*exp(-r*tau)*P(lambda, vbar, eta, rho, v0, r, tau, S0, K, 0)
  A-B
}

hestoncall.mc <- function(r, rho, kappa, theta, sigma, u, vu, Su, t, K, n, m, seed=1) {
  payoff <- function(s) exp(-r * tau) * max(s-K, 0)
  est <- rep(0.0, m)
  for (i in 1:m) {
    p <- matrix(runif(n*3),ncol=3)
    V <- mapply(payoff, qhestonexact(p[,1], p[,2], p[,3], r, rho, kappa, theta, sigma, u, vu, Su, t))
    est[i] <- mean(V)
  }
  m <- mean(est)
  stderr <- sqrt(sum((est-m)^2)/(m*(m-1)))
  return(list(value=m, stderr=stderr)) 
}

hestoncall.sobol <- function(r, rho, kappa, theta, sigma, u, vu, Su, t, K, n, m, seed=1) {
  payoff <- function(s) exp(-r * tau) * max(s-K, 0)
  est <- rep(0.0, m)
  for (i in 1:m) {
    p <- runif.sobol(n,3,scrambling=1, seed=seed+i)
    V <- mapply(payoff, qhestonexact(p[,1], p[,2], p[,3], r, rho, kappa, theta, sigma, u, vu, Su, t))
    est[i] <- mean(V)
  }
  m <- mean(est)
  stderr <- sqrt(sum((est-m)^2)/(m*(m-1)))
  return(list(value=m, stderr=stderr)) 
}

rhestonnaive <- function(r, rho, kappa, theta, sigma, u, vu, Su, t, nPaths, nSteps=30, vneg=2) {
  v0 <- vu
  S0 <- Su
  tau <- (t-u)
  
  n <- nSteps
  N <- nPaths
  
  dt <- tau / n
  
  negCount <- 0
  
  S <- rep(S0,N)
  v <- rep(v0,N)
  
  for (i in 1:n) 
  {
    W1 <- rnorm(N);
    W2 <- rnorm(N);
    W2 <- rho*W1 + sqrt(1 - rho^2)*W2;

    sqvdt <- sqrt(v*dt)
    S <- S*exp((r-v/2)*dt + sqrt(v * dt) * W1)
    
    if ((vneg == 3) & (2*kappa*theta/(sigma^2) <= 1)) {
        cat("Variance not guaranteed to be positive with choice of kappa, theta, and sigma\n")
        cat("Defaulting to Reflection + Milstein method\n")
        vneg = 2
    }

    if (vneg == 0){
      ## Absorbing condition
      v <- v + kappa*(theta - v)* dt + sigma * sqvdt * W2
      negCount <- negCount + length(v[v < 0])
      v[v < 0] <- 0
    }
    if (vneg == 1){
      # Reflecting condition
      sqvdt <- sqrt(v*dt)
      v <- v + kappa*(theta - v)* dt + sigma * sqvdt * W2
      negCount <- negCount + length(v[v < 0])
      v <- ifelse(v<0, -v, v)
    }
    if (vneg == 2) {
      # Reflecting condition + Milstein
      v <- (sqrt(v) + sigma/2*sqrt(dt)*W2)^2 - kappa*(v-theta)*dt - sigma^2/4*dt
      negCount <- negCount + length(v[v < 0])
      v <- ifelse(v<0, -v, v)     
    }
    if (vneg == 3) {
      # Alfonsi - See Gatheral p.23
      v <- v -kappa*(v-theta)*dt +sigma*sqrt(v*dt)*W2 - sigma^2/2*dt      
    }
  }

  S
}

test_histograms <- function(n=2048) {
	kappa <- 6.21
	theta <- 0.019
	sigma <- 0.61
	
	u <- 0
	vu <- 0.010201
	tau <- 1.0
	vt <- 0.02
	r <- 0.0319
	Su <- 100
	K <- 100
	rho <- -0.7
	
	# Heston using naive euler paths
	rhestonnaive(r, rho, kappa, theta, sigma, u, vu, Su, tau, n) -> S1
	hist(S1,breaks=100,freq=F,col="cyan3",main="Heston law - Euler")
	
	# Heston using exact simulation
	rhestonexact(r, rho, kappa, theta, sigma, u, vu, Su, tau, n) -> S2
	hist(S2,breaks=100,freq=F,col="cyan3",main="Heston law - Exact")
	
	# Heston using exact simulation with Sobol
	rhestonexact.sobol(r, rho, kappa, theta, sigma, u, vu, Su, tau, n) -> S3
	hist(S3,breaks=100,freq=F,col="cyan3",main="Heston law - Exact Sobol")
	
	# compare histograms
	hist(S2,breaks=100,freq=T,col="wheat3",main=NA)
	hist(S1,breaks=100,freq=T,col="cyan3",main=NA,add=T)
	hist(S3,breaks=100,freq=T,col="darkred",main=NA,add=T)
}

test_eurocalls <- function(n=2048) {
	kappa <- 6.21
	theta <- 0.019
	sigma <- 0.61
	
	u <- 0
	vu <- 0.010201
	tau <- 1.0
	vt <- 0.02
	r <- 0.0319
	Su <- 100
	K <- 100
	rho <- -0.7

	hestoncall(lambda=kappa, vbar=theta,  eta=sigma, rho=rho, v0=vu, r=r, tau=tau, S0=Su, K=K) -> sol1
	hestoncall.mc(r,rho,kappa,theta,sigma,u,vu,Su,tau,K,512,4, seed=2) -> sol2
	hestoncall.sobol(r,rho,kappa,theta,sigma,u,vu,Su,tau,K,512,4, seed=2) -> sol3
}

test_root <- function(p) {
  kappa <- 6.21
	theta <- 0.019
	sigma <- 0.61
	
	u <- 0
	vu <- 0.010201
	tau <- 1.0
	vt <- 0.02
	r <- 0.0319
	Su <- 100
	K <- 100
	rho <- -0.7
  
  qintsqrbridge(p, kappa, theta, sigma, u, vu, tau, vt)
}


test_qhestonexact <- function(p1,p2,p3) {
  kappa <- 6.21
  theta <- 0.019
	sigma <- 0.61
	
	u <- 0
	vu <- 0.010201
	tau <- 1.0
	vt <- 0.02
	r <- 0.0319
	Su <- 100
	K <- 100
	rho <- -0.7
  
  qhestonexact(p1, p2, p3, r, rho, kappa, theta, sigma, u, vu, Su, tau)
}

#test_qhestonexact(0.5,0.5,0.5)

