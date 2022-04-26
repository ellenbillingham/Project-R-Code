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

par(mfrow=c(1,1))
persp(xxb, xxg, matrix(mup, ncol=length(xxb)), theta=200, phi=15, xlab="beta", ylab="gamma", zlab="Infectious at t=5", ticktype="detailed")
res = persp(xxb, xxg, matrix(mup, ncol=length(xxb)), theta=200, phi=15, xlab="beta", ylab="gamma", zlab="Infectious at t=5", ticktype="detailed")
datapoints = trans3d(x=X[,1], y=X[,2], z=Y, pmat=res)
points(datapoints, pch=20, cex=2)
# sort of "ridge" around beta=5ish
# ie infections at t=5 seem to be highest for beta close to 5 at all gammas?

# with fewer testing points the surface looks very different
# lots of ups and downs