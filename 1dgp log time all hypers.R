library(deSolve)

sir = function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS = -beta * S * I
    dI =  beta * S * I - gamma * I
    dR =                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# create initial values and times
init = c(S=1-1e-6, I=1e-6, R=0.0)
times = seq(0, 100, by=0.1)

library(plgp)

mylhs <- function(n, m)
{
  ## generate the Latin hypercube 
  l <- (-(n - 1)/2):((n - 1)/2)
  L <- matrix(NA, nrow=n, ncol=m)
  for(j in 1:m) L[,j] <- sample(l, n)
  
  ## draw the random uniforms and turn the hypercube into a sample
  U <- matrix(runif(n*m), ncol=m)
  X <- (L + (n - 1)/2 + U)/n
  colnames(X) <- paste0("x", 1:m)
  
  ## return the design and the grid it lives on for visualization
  return(list(X=X, g=c((l + (n - 1)/2)/n,1)))
}

#n = 10
#X = 8.5*mylhs(n,1)$X + 1.5
n = 10
X = c(seq(1.5, 10, length=n))
Y = c()
for (b in X) {
  out = ode(y=init, times=times, func=sir, parms=c(beta=b, gamma=1))
  out = as.data.frame(out)
  #out$time = NULL
  Y = c(Y, log(out$time[out$I==max(out$I)]))
}

eps = sqrt(.Machine$double.eps)

X = matrix(X)
X = rbind(X, X) # rbinding X more than twice improves the estimate
n = nrow(X)
Y = Y + rnorm(n, sd=0.05) # how to choose this value of sd???
D = distance(X)

nl <- function(par, D, Y) 
{
  theta <- par[1]                                       ## change 1
  g <- par[2]
  n <- length(Y)
  K <- exp(-D/theta) + diag(g, n)                       ## change 2
  Ki <- solve(K)
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
  counter <<- counter + 1
  return(-ll)
}

gradnl <- function(par, D, Y)
{
  ## extract parameters
  theta <- par[1]
  g <- par[2]
  
  ## calculate covariance quantities from data and parameters
  n <- length(Y)
  K <- exp(-D/theta) + diag(g, n)
  Ki <- solve(K)
  dotK <- K*D/theta^2
  KiY <- Ki %*% Y
  
  ## theta component
  dlltheta <- (n/2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - 
    (1/2)*sum(diag(Ki %*% dotK))
  
  ## g component
  dllg <- (n/2) * t(KiY) %*% KiY / (t(Y) %*% KiY) - (1/2)*sum(diag(Ki))
  
  ## combine the components into a gradient vector
  return(-c(dlltheta, dllg))
}

outg <- optim(c(1.375, 0.01*var(Y)), nl, gradnl, method="L-BFGS-B", 
              lower=eps, upper=c(1.4375, var(Y)), D=D, Y=Y)

K <- exp(- D/outg$par[1]) + diag(outg$par[2], n)
Ki <- solve(K)
tau2hat <- drop(t(Y) %*% Ki %*% Y / n)

XX = matrix(seq(1.1, 11, length=100), ncol=1)
DXX <- distance(XX)
KXX <- exp(-DXX/outg$par[1]) + diag(outg$par[2], ncol(DXX))
DX <- distance(XX, X)
KX <- exp(-DX/outg$par[1])

mup <- KX %*% Ki %*% Y
Sigmap <- tau2hat*(KXX - KX %*% Ki %*% t(KX))
sdp <- sqrt(diag(Sigmap))

trueY = c()
for (b in XX) {
  out = ode(y=init, times=times, func=sir, parms=c(beta=b, gamma=1))
  out = as.data.frame(out)
  #out$time = NULL
  trueY = c(trueY, log(out$time[out$I==max(out$I)]))
}

YY = rmvnorm(100, mup, Sigmap)
q1 = mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 = mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))

matplot(XX, t(YY), type="l", lty=1, col="gray", xlab=expression(beta), ylab="log(time of peak infections)")
points(X, Y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, trueY, lwd=2, col="green")
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)

legend("bottomleft", c("Emulator mean", "Emulations", "True values", "95% Uncertainty"),
       lty =c(1, 1, 1, 2), lwd = c(2, 1, 2, 2), col=c("black", "gray", "green", "red"),
       bty = "n")

title(main="1D Emulator with Hyperparameters")

# sometimes gives a really nice plot with the emulations very close
# to the true value, when g approx 4 and theta approx 0.1
# but sometimes looks horrible