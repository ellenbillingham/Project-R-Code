# include libraries
library(deSolve)
library(minpack.lm)
library(reshape2)
library(Rcpp)

# read in the uk data
ukdata = read.csv("C:/Users/elkam/Dropbox/Ellie/university/project/overview_2022-01-18.csv")

# format the dates and day numbers
ukdata$date = as.Date(ukdata$date, "%d/%m/%Y")
ukdata$DayNumber = seq.int(nrow(ukdata))

# take data of just days 30 to 180
uk1 = ukdata[30:180,]

# set initial constants
N = 67000000
I0 = uk1$newCasesBySpecimenDateRollingSum[1]
R0 = 0.0
S0 = N - I0 - R0
init = c(S=S0, I=I0, Tr=0.0, R=R0)
times = uk1$DayNumber

# sitr model function
sitr = function(time, state, parameters, alpha, delta) {
  beta = as.numeric(parameters$beta) #infection rate from I to S
  #delta = as.numeric(parameters$delta) #people in T can infect people in S at reduced rate delta*beta
  #alpha = as.numeric(parameters$alpha) # fraction per unit time of infected people selected for treatment
  gamma = as.numeric(parameters$gamma) #recovery rate without treatment
  
  if (alpha==0.0) {
    nu = 0.0 #recovery rate with treatment
  }
  else {
    nu = as.numeric(parameters$nu)
  }
  
  with(as.list(c(state, parameters)), {
    
    dS = -beta * S * I / N - delta * beta * S * Tr / N
    dI = beta * S * I / N + delta * beta * S * Tr / N - alpha * I - gamma * I
    dT = alpha * I - nu * Tr
    dR = gamma * I + nu * Tr
    
    return(list(c(dS, dI, dT, dR)))
  })
}

# function to calculate residuals of model with current parameter values compared to the data
ssq = function(par, dat, alpha, delta) {
  out = ode(y=init, times=times, func=sitr, parms=par, alpha=alpha, delta=delta)
  outuk = data.frame(out)
  preduk = melt(outuk, id.var="time", measure.vars = "I", value.name="number")
  expuk = melt(dat, id.var="DayNumber", measure.vars = "newCasesBySpecimenDateRollingSum", value.name="number")
  ssqres = preduk$number - expuk$number
  
  return(ssqres)
}

# function that plots sitr model with fixed value of delta (probability of infected individual being given treatment)
plot_sitr = function(alpha, delta, casedata) {
  
  if (alpha==0.0) {
    parameters = list(beta=2, gamma=1)
  }
  else {
    parameters = list(beta=2, gamma=1, nu=1)
  }
  
  fitval = nls.lm(par=parameters, fn=ssq, dat=casedata, alpha=alpha, delta=delta)
  parameter_estimate = as.list(coef(summary(fitval))[,"Estimate"])
  #print(parameter_estimate)
  
  out_new = ode(y=init, times=times, func=sitr, parms=parameter_estimate, alpha=alpha, delta=delta)
  out_new = as.data.frame(out_new)
  
  plot(x=times, y=out_new$I, type="l", main=paste("alpha = ", alpha, ", delta = ", delta), lwd=1, lty=1, bty="l", col = 2, ylim=c(0,45000),
       ylab="", xlab="")
  points(x=times, y=casedata$newCasesBySpecimenDateRollingSum)
}

# try it:
#plot_sitr(0.5, 0.5, uk1)

# (FIXED!!!) problem - the model should be exactly SIR when alpha=0
# but the graph is not showing the same as the SIR in this case


# problem when alpha>0 and delta=0 - not sure what the issue is
# logically this should work - just means that the treated patients
# don't infect the susceptible people at all

par(mfrow=c(3,3))
for (a in seq(0.1,0.9,length=3)) {
  for (d in seq(0.1,0.9,length=3)) {
    plot_sitr(a,d,uk1)
  }
}
