# generate an approximation of the SIR infectious curve for
# beta=5.5, gamma=1

# from the emulators, get estimate and uncertainty intervals for
# time and max infections

yest = 0.469036 # estimate of max infections
yci = c(0.4459551, 0.5314377) # 95% confidence interval for max infections

xest = exp(1.292612) # estimate of time of max
xci = exp(c(1.028836,1.450535)) # 95% confidence interval for time

# make a graph with the box formed by the uncertainty intervals
# and the point with the estimates as coordinates

plot(xest,yest,pch=20,cex=2,xlim=c(0,15),ylim=c(0,0.6), xlab="Time",
     ylab="Proportion of population infected")


# then draw a straight line from (0,0) to a random point in that box
# and symmetrically decline after this peak, then just be 0 when it
# reaches y=0

x = runif(5, xci[1], xci[2])
y = runif(5, yci[1], yci[2])

for (i in 1:length(x)) {
  lines(c(0,x[i]), c(0,y[i]), col="grey")
  lines(c(x[i],2*x[i]), c(y[i],0), col="grey")
  lines(c(2*x[i],15), c(0,0), col="grey")
}

# repeat a few times to get a few different approximations

# plot the true infectious curve from SIR to see how close it is

library(deSolve)

sir = function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS = -beta * S * I
    dI =  beta * S * I - gamma * I
    dR =                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

init = c(S=1-1e-6, I=1e-6, R=0.0)
parameters = c(beta=5.5, gamma=1)
times = seq(0, 15, by=0.1)
out = ode(y=init, times=times, func=sir, parms=parameters)
out = as.data.frame(out)
out$time = NULL
lines(times,out$I, lwd=2)

rect(xci[1],yci[1],xci[2],yci[2])

legend(10,0.4, c("95% Uncertainty","Approximations", "True Values"),
       lty=1, lwd=c(1,1,2), col=c(1,"grey", 1))