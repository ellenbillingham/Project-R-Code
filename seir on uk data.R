ukdata = read.csv("C:/Users/elkam/Dropbox/Ellie/university/project/overview_2022-01-18.csv")

ukdata$date = as.Date(ukdata$date, "%d/%m/%Y")
ukdata$DayNumber = seq.int(nrow(ukdata))
#ukdata = ukdata[order(ukdata$date),]

#plot(ukdata$date, ukdata$newCasesBySpecimenDateRollingSum)

# divide all the case numbers by total population (approx 67 million)
#ukdata$casesProp = ukdata$newCasesBySpecimenDateRollingSum / 67000000

#plot(ukdata$date, ukdata$casesProp)

library(deSolve)
library(minpack.lm)
library(reshape2)
library(Rcpp)

seir = function(time, state, parameters) {
  beta = as.numeric(parameters$beta)
  eta = as.numeric(parameters$eta)
  gamma = as.numeric(parameters$gamma)
  
  with(as.list(c(state, parameters)), {
    
    dS = -beta * S * I / N
    dE = beta * S * I / N - eta * E
    dI = eta * E - gamma * I
    dR = gamma * I
    
    return(list(c(dS, dE, dI, dR)))
  })
}

N = 67000000
I0 = ukdata$newCasesBySpecimenDateRollingSum[1]
R0 = 0.0
S0 = N - I0 - R0

init = c(S=S0, E=0.0, I=I0, R=R0)
times = 1:nrow(ukdata)

ssq = function(par) {
  beta = par$beta
  eta = par$eta
  gamma = par$gamma
  out = ode(y=init, times=times, func=seir, parms=list(beta=beta, eta=eta, gamma=gamma))
  outuk = data.frame(out)
  preduk = melt(outuk, id.var="time", measure.vars = "I", value.name="number")
  expuk = melt(ukdata, id.var="DayNumber", measure.vars = "newCasesBySpecimenDateRollingSum", value.name="number")
  ssqres = preduk$number - expuk$number
  
  return(ssqres)
}

parameters = list(beta=2, eta=1, gamma=1)
fitval = nls.lm(par=parameters, fn=ssq)
summary(fitval)

parameter_estimate = as.list(coef(summary(fitval))[,"Estimate"])

out_new = ode(y=init, times=times, func=seir, parms=parameter_estimate)
out_new = as.data.frame(out_new)

#plot(x=times, y=out_new$I, type="l", xlab="Time", ylab="Number of people infected", 
#     main="SEIR Model on UK Data", lwd=1, lty=1, bty="l", col = 2, ylim = c(0,1500000))
#points(x=ukdata$DayNumber, y=ukdata$newCasesBySpecimenDateRollingSum)


# Try on just the first 6 months to see if it models the start of the pandemic well:

uk1 = ukdata[26:117,]

N = 67000000
I0 = uk1$newCasesBySpecimenDateRollingSum[1]
R0 = 0.0
S0 = N - I0 - R0

init = c(S=S0, E=0.0, I=I0, R=R0)
times = 26:117

ssq1 = function(par) {
  beta = par$beta
  eta = par$eta
  gamma = par$gamma
  out = ode(y=init, times=times, func=seir, parms=list(beta=beta, eta=eta, gamma=gamma))
  outuk = data.frame(out)
  preduk = melt(outuk, id.var="time", measure.vars = "I", value.name="number")
  expuk = melt(uk1, id.var="DayNumber", measure.vars = "newCasesBySpecimenDateRollingSum", value.name="number")
  ssqres = preduk$number - expuk$number
  
  return(ssqres)
}

parameters = list(beta=2, eta=0.2, gamma=1)
fitval = nls.lm(par=parameters, fn=ssq1)
summary(fitval)

beta_est = coef(summary(fitval))["beta","Estimate"]
eta_est = coef(summary(fitval))["eta","Estimate"]
gamma_est = coef(summary(fitval))["gamma","Estimate"]

parameter_estimate = list(beta = beta_est, eta = eta_est, gamma = gamma_est)

out_new = ode(y=init, times=times, func=seir, parms=parameter_estimate)
out_new = as.data.frame(out_new)

plot(x=times, y=out_new$I, type="l", xlab="Time", ylab="Number of people infected", 
     main="SEIR Model on UK Data", lwd=1, lty=1, bty="l", col = 2)
points(x=uk1$DayNumber, y=uk1$newCasesBySpecimenDateRollingSum)

# not a great fit but does show a peak of infections at approximately
# the same time as the data