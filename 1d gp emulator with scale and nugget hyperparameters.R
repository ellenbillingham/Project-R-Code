# SIR model (code from webpage)

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
times = seq(0, 15, by=1)


#-------- 1D EMULATION - FIX VALUE OF GAMMA AND TEST AT DIFFERENT BETAS ----
# follow similar approach to section 5.1.1 in surrogates book
library(plgp)

n = 10
X = c(seq(1, 10, length=n))
Y = c()
for (b in X) {
  out = ode(y=init, times=times, func=sir, parms=c(beta=b, gamma=1))
  out = as.data.frame(out)
  out$time = NULL
  Y = c(Y, out$I[6])
}
D = distance(X)
eps = sqrt(.Machine$double.eps)
Sigma = exp(-D) #+ diag(eps, ncol(D))

XX = matrix(seq(0, 11, length=100), ncol=1)

nlg = function(g, D1, Y1)
{
  m = length(Y1)
  K = exp(-D1) + diag(g, m)
  Ki = solve(K)
  ldetK = determinant(K, logarithm=TRUE)$modulus
  ll = - (n/2)*log(t(Y1) %*% Ki %*% Y1) - (1/2)*ldetK
  counter <- counter + 1
  return(-ll)
}

X = matrix(X)
X = rbind(X, X) # rbinding X more than twice improves the estimate
n = nrow(X)
Y = Y + rnorm(n, sd=0.02) # how to choose this value of sd???
D = distance(X)

counter = 0
g = optimise(nlg, interval=c(eps,var(Y)), D1=D, Y1=Y)$minimum

K = exp(-D) + diag(g, n)
Ki = solve(K)
tau2hat = drop(t(Y) %*% Ki %*% Y / n)

DX = distance(XX,X)
KX = exp(-DX)
DXX = distance(XX)
KXX = exp(-DXX) + diag(g, nrow(DXX))

mup = KX %*% Ki %*% Y
Sigmap = tau2hat * (KXX - KX %*% Ki %*% t(KX))
q1 = mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 = mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))

Sigma.int = tau2hat * (exp(-DXX) + diag(eps, nrow(DXX)) - KX %*% Ki %*% t(KX))
YY = rmvnorm(100, mup, Sigma.int)

trueY = c()
for (b in XX) {
  out = ode(y=init, times=times, func=sir, parms=c(beta=b, gamma=1))
  out = as.data.frame(out)
  out$time = NULL
  trueY = c(trueY, out$I[6])
}

matplot(XX, t(YY), type="l", lty=1, col="gray", xlab="beta", ylab="Infectious at  t=5")
points(X, Y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, trueY, lwd=2, col="green")
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)

legend(5, -0.2, c("Emulator mean", "Emulations", "True values", "95% Uncertainty"),
       lty =c(1, 1, 1, 2), lwd = c(2, 1, 2, 2), col=c("black", "gray", "green", "red"),
       bty = "n")
