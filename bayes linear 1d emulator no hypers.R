#specify the beta values to be used as inputs:
n = 10
xj = c(seq(1,10,length=n))

#get the model output values from SIR model with inputs beta = xj:
library(deSolve)

sir = function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS = -beta * S * I
    dI =  beta * S * I - gamma * I
    dR =                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

infections = function(XX, g){
  trueY = c()
  init = c(S=1-1e-6, I=1e-6, R=0.0)
  times = seq(0, 15, by=1)
  for (b in XX) {
    out = ode(y=init, times=times, func=sir, parms=c(beta=b, gamma=g))
    out = as.data.frame(out)
    out$time = NULL
    trueY = c(trueY, out$I[6])
  }
  return(trueY)
}

D = infections(xj, 1)


#aim to approximate f(x) = proportion of population infected at t=5

#choose the polynomial terms to represent a simple linear model, b0 + b1*x
#also for now assume nugget is 0

#so emulator is given by f(x) = b0 + b1*x + u(x)
#treat b0 as known, choose prior expectation E(b0)=b0=0.5 (might need to change this)
#priorfit = lm(D ~ poly(xj,3,raw=TRUE))

priorfit = lm(D ~ 1)

#choose sigma_u and theta (starting with sigma_u = 1.5 and theta=0.14 same as paper)
#probably will need to change these
sigma_u = 1 #0.15
theta = 1 #1.25

E_fx = function(x) {
  return(predict(priorfit, data.frame(xj=x)))
}

Var_fx = sigma_u ^ 2

E_D = function(x) {
  return(c(seq(E_fx(x),E_fx(x),length=n))) 
}

Cov_fx_D = function(x) {
  cov = c()
  for (j in 1:n) {
    cov = c(cov, sigma_u^2 * exp(-((abs(x-xj[j])^2)/theta^2)))
  }
  return(cov)
}

Var_D = matrix(data = 0, nrow = n, ncol = n)
for (j in 1:n) {
  for (k in 1:n) {
    Var_D[j,k] = sigma_u^2 * exp(-((abs(xj[j]-xj[k])^2)/theta^2))
  }
  
}

# use these to calculate the adjusted expectation and variance:

ED_fx = function(x) {
  return(E_fx(x) + Cov_fx_D(x) %*% solve(Var_D) %*% (D - E_D(x)))
}

Var_D_fx = function(x) {
  if (x %in% xj) {
    return(0.0)
  }
  else {
    return(Var_fx - Cov_fx_D(x) %*% solve(Var_D) %*% Cov_fx_D(x)) #think the error is with the last bit in the multiplication
  }
}

# plot the findings:
XX = matrix(seq(0, 11, length=100), ncol=1)
YY = matrix(seq(0, 0, length=100), ncol=1)
for (i in 1:100) {
  YY[i,1] = ED_fx(XX[i,1])
}
trueY = infections(XX, 1)

Var_credible = matrix(seq(0, 0, length=100), ncol=1)
for (i in 1:100) {
  Var_credible[i,1] = Var_D_fx(XX[i,1])
}

c1 = YY - 3 * sqrt(Var_credible)
c2 = YY + 3 * sqrt(Var_credible)

matplot(XX, YY, type="l", lty=1, col="black", lwd=3, xlab="beta", ylab="Infectious at t=5", ylim=c(-2.5,2.5))
points(xj, D, pch=20, cex=2)
#lines(XX, predict(priorfit, data.frame(xj=XX)), col="black")
lines(XX, trueY, lwd=2, col="green")
lines(XX, c1, lwd=1, col="red", lty=2)
lines(XX, c2, lwd=1, col="red", lty=2)

legend("bottomright",legend=c("Updated expectation", "True values", "Credible region"), 
       col=c("black", "green", "red"),lty=c(1,1,2),lwd=c(3,2,1))

# can see from this plot that choosing a cubic prior gives a better fit than constant but still not great

#SOMETHING WRONG WITH VARIANCE CALCULATION --> CREDIBLE REGIONS INCORRECT?
