library(deSolve)
library(reshape2)
library(ggplot2)

## Differential equation for SIR-R model of COVID-19 ##
ode_SIRS <- function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    #region = parameters$region
    InfectionParam = parameters$InfectionParam
    #VirusTransParam = parameters$VirusTransParam
    # Disaggregate fixed parameters
    beta0 = InfectionParam[1]
    g = InfectionParam[2]
    z = InfectionParam[3]
    # daily beta multiplier
    w_idx = ceiling(t/7) # find week index
    m.beta.d <- m.beta[w_idx]
    print(c(t,w_idx))
    #print(m.beta.d)
    m.beta.d=1
    
    beta.d = beta0 * m.beta.d
    # ODE equations
    #S=state[1]
    #I=state[2]
    #R=state[3]
    dS = -beta.d*S*I + z*R
    dI = beta.d*S*I - g*I
    dR = g*I - z*R
    dC = beta.d*S*I
    
    
    list(c(dS,dI,dR,dC))
  })
}

## A function to run SIR-R model of COVID-19 ##
run_SIRS <- function(region){
  with(parameters, {
    timestep=1
    # time span
    tspan = seq(T0, T0+Tend, by=timestep) 
    # Initial distribution by region
    if (region == 'Green'){ # Northeast
      # initial distribution
      IC = c(S=1-1000/71000000,I=1000/71000000,R=0,C=0)
      # region-specific oscillator
    }else if(region == 'Blue'){ # South
      IC = c(S=1-1000/67000000,I=1000/67000000,R=0,C=0)
    }else if (region == 'Red'){ # Midwest
      IC = c(S=1-1000/187000000,I=1000/187000000,R=0,C=0)
    }
  
    # run ODE solver
    odetime <- proc.time()
    #prob = rk(y=IC, times=tspan, func = ode_SIRS, parms = parameters, method='rk45dp7')
    prob = rk(y = IC, times = tspan, func = ode_SIRS, parms = parameters, method='rk45dp7',atol = 10^-8, rtol = 10^-10)  # TODO: constrain domain for x to be 0 or positive (x>0)
    proc.time() - odetime
    return(prob)
  })
} 

## A function to calculate the new cases ##
cal_case <- function(prob){
  new_case <-prob[-1,5]- prob[-nrow(prob),5]
  return(new_case)
}

## A function to calculate error ##
cal_cost <- function(region, state, parameters){
  
  # model outcome
  prob <- run_SIRS(region, parameters)
  I.m <- cal_case(prob)
  # observed case
  I.o <- covid_case %>% 
    filter(state == state)
  # sum of squared errors

  cost <- sum((I.m - I.o)^2)
  
  return(cost)
}

## Main ##
# Model set up parameters
T0 = 0.0
Tend = 730
n_days = as.integer(Tend-T0)

# Fixed infection parameters
beta0=0.185
gamma= 0.16
z = 0.005
# varying parameters (weekly beta multiplier)
m.beta = rep(1,as.integer(n_days/7))

# observational data


## Fitting daily beta by region
# parameters to run SIRS model
parameters = list(region = region, 
                  T0 = T0,
                  Tend = Tend,
                  n_days = n_days,
                  InfectionParam=c(beta0,gamma,z),
                  m.beta = m.beta)

optim(m.beta,fn = cal_cost)
# check beta.d over time
t<-seq(0,Tend,by=1)




# Test run with IC = c(1,0.01,0)
I.all <- data.frame(time=seq(0,Tend),green=prob.gr[,3],blue=prob.bl[,3],red=prob.rd[,3])
I.all.t <- melt(I.all,id.vars='time')
ggplot(data=I.all.t,aes(time,value, color=variable))+
  geom_line()+
  scale_color_manual(values=c("green","blue","red"))

