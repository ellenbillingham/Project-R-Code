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
Sigma = exp(-D) + diag(eps, ncol(D))

XX = matrix(seq(0, 11, length=100), ncol=1)
DXX = distance(XX)
SXX = exp(-DXX) + diag(eps, ncol(DXX))
DX = distance(XX, X)
SX = exp(-DX)

Si = solve(Sigma)
mup = SX %*% Si %*% Y
Sigmap = SXX - SX %*% Si %*% t(SX)

YY = rmvnorm(100, mup, Sigmap)
q1 = mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 = mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))

matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="beta", ylab="Infectious at t=5")
points(X, Y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)

trueX = c(seq(0,11,length=100))
trueY = c()
for (b in trueX) {
  out = ode(y=init, times=times, func=sir, parms=c(beta=b, gamma=1))
  out = as.data.frame(out)
  out$time = NULL
  trueY = c(trueY, out$I[6])
}
lines(trueX, trueY, lwd=2, col="green")

legend(6, -1, c("Emulator mean", "Emulations", "True values", "95% Uncertainty"),
       lty =c(1, 1, 1, 2), lwd = c(2, 1, 2, 2), col=c("black", "gray", "green", "red"),
       bty = "n")
