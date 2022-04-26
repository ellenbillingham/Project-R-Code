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

#-------- 2D EMULATION - TEST ON VALUES OF BETA AND GAMMA ------------------
# follow similar approach to section 5.1.2 in surrogates book
library(plgp)

xb = c(seq(1, 10, length=8))
xg = c(seq(0.1, 2, length=8))
X = expand.grid(xb, xg)
Y = c()

# generate the outputs --> I values at time 5 for all values of beta and gamma

for (i in 1:64) {
  out = ode(y=init, times=times, func=sir, parms=c(beta=X[i,1], gamma=X[i,2]))
  out = as.data.frame(out)
  out$time = NULL
  Y = c(Y, out$I[6])
}

xxb = seq(0.5, 10.5, length = 40)
xxg = seq(0.05, 2.05, length = 40)
XX = expand.grid(xxb, xxg)

D = distance(X)
Sigma = exp(-D)

eps = sqrt(.Machine$double.eps)

DXX = distance(XX)
SXX = exp(-DXX) + diag(eps, ncol(DXX))
DX = distance(XX, X)
SX = exp(-DX)

Si = solve(Sigma)
mup = SX %*% Si %*% Y
Sigmap = SXX - SX %*% Si %*% t(SX)

sdp = sqrt(diag(Sigmap))

# create heat plots of the results
# left panel - image plot of mean over regularly gridded inputs XX
# right panel - image plot of sd over regularly gridded inputs XX
par(mfrow=c(1,2))
cols = heat.colors(128)
image(xxb, xxg, matrix(mup, ncol=length(xxb)), xlab="beta", ylab="gamma", col=cols)
points(X[,1], X[,2])
image(xxb, xxg, matrix(sdp, ncol=length(xxb)), xlab="beta", ylab="gamma", col=cols)
points(X[,1], X[,2])

