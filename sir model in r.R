# basic SIR model - nothing added

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
parameters = c(beta=2, gamma=1)
times = seq(0, 20, by=0.1)

out = ode(y=init, times=times, func=sir, parms=parameters)
out = as.data.frame(out)
out$time = NULL
#head(out, 10)

matplot(x=times, y=out$I, type="l", xlab="Time", ylab="Proportion of population infected", 
      lwd=2, lty=1, bty="l", col=2, ylim=c(0,0.6))

for (b in 2:6) {
  parameters = c(beta=b, gamma=1)
  out = ode(y=init, times=times, func=sir, parms=parameters)
  out = as.data.frame(out)
  out$time = NULL
  lines(x=times,y=out$I,type="l",lwd=2,bty="l",col=b)
}

legend(15, 0.4, sapply(2:6, function(x) as.expression(substitute(beta == B,
                                                                 list(B = as.name(x))))), 
       col=2:6, lty=1, lwd=2, bty="n")
#for (b in seq(1, 10, length=8)) {
#  out = ode(y=init, times=times, func=sir, parms=c(beta=b,gamma=1.8))
#  out = as.data.frame(out)
#  out$time = NULL
#  lines(x=times,y=out$I)
#}

# SIR with vaccination

library(deSolve)

sir = function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS = -beta * S * I             - vac * S
    dI =  beta * S * I - gamma * I
    dR =                 gamma * I + vac * S
    
    return(list(c(dS, dI, dR)))
  })
}

init = c(S=1-1e-6, I=1e-6, R=0.0)
parameters = c(beta=1.4247, gamma=0.14286, vac = 0.02)
times = seq(0, 50, by=1)

out = ode(y=init, times=times, func=sir, parms=parameters)
out = as.data.frame(out)
out$time = NULL
head(out, 10)

matplot(x=times, y=out, type="l", xlab="Time", ylab="Susceptible and Recovered", 
        main="SIR Model", lwd=1, lty=1, bty="l", col=2:4)

legend(40, 0.7, c("Susceptible", "Infected", "Recovered"), pch=1, 
       col=2:4, bty="n")

# basically to add something in, just adjust the differential equations in sir()
# then add in any extra parameter values into parameters

# e.g. SEIR:

seir = function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS = -beta * S * I
    dE =  beta * S * I - eta * E
    dI =                 eta * E - gamma * I
    dR =                           gamma * I
    
    return(list(c(dS, dE, dI, dR)))
  })
}

init = c(S=1-1e-6, E=0.0, I=1e-6, R=0.0)
parameters = c(beta=1.4, eta=0.2, gamma=0.2)
times = seq(0, 70, by=1)

out = ode(y=init, times=times, func=seir, parms=parameters)
out = as.data.frame(out)
out$time = NULL

matplot(x=times, y=out, type="l", xlab="Time", ylab="Proportion of Population", 
        main="SEIR Model", lwd=1, lty=1, bty="l", col=c(2,7,3,4))

legend(5, 0.6, c("Susceptible", "Exposed", "Infected", "Recovered"), lty=1, 
       col=c(2,7,3,4), bty="l", lwd=1)
