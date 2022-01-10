library(deSolve)
library(reshape2)
library(ggplot2)


## Differential equation for SIR-R model of COVID-19 ##
ode_SIRS <- function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    #region = parameters$region
    InfectionParam = parameters$InfectionParam
    NAOParam = parameters$NAOParam
    RSOParam = parameters$RSOParam
    VirusTransParam = parameters$VirusTransParam
    # Disaggregate fixed parameters
    beta0 = InfectionParam[1]
    g = InfectionParam[2]
    z = InfectionParam[3]
    # Disaggregate parameters for NAO
    c.nao = NAOParam[1]
    a.nao = NAOParam[2]
    p.nao = NAOParam[3]
    f.nao = NAOParam[4]
    p.rso = RSOParam[region,5]
    f.rso = RSOParam[region,6]
    
    # Caculate daily NAO multiplier
    m.nao.d = a.nao * sin((t-p.nao)/(1/f.nao))+c.nao 
    
    # Disaggregate parameters for changing viral transmissibility
    L = VirusTransParam[1]
    K = VirusTransParam[2]
    
    # Disaggregate parameters for RSO
    if (region %in% c("Green","Red")){
      if ((t%/%183)%%2 == 0){ # alternating amplitude for every half-year
        c.rso = RSOParam[region,1]
        a.rso = RSOParam[region,3]
      }else{
        c.rso = RSOParam[region,2]
        a.rso = RSOParam[region,4]}
      # Calculate daily RSO multiplier
      m.rso.d = a.rso * sin((t-p.rso)/(1/f.rso))+c.rso 
    } else{ # Blue states
      if (t<=91 | ((t+91)%/%183)%%2 == 0){ # alternating amplitude for every half-year
        c.rso = rso.matrix[region,1]
        a.rso = rso.matrix[region,3]
      }else{
        c.rso = rso.matrix[region,2]
        a.rso = rso.matrix[region,4]}
      m.rso.d= a.rso * sin((t+p.rso)/(1/f.rso))+c.rso
    }

    # Calculate daily viral transmissibility multiplier
    j.d = L + K*t
    # Calculate daily beta
    beta.d = beta0 * ((m.nao.d+m.rso.d)/2*j.d)
    # ODE equations
    #S=state[1]
    #I=state[2]
    #R=state[3]
    dS = -beta.d*S*I + z*R
    dI = beta.d*S*I - g*I
    dR = g*I - z*R
  

  list(c(dS,dI,dR))
  })
}

## A function to run SIR-R model of COVID-19 ##
run_SIRS <- function(region){
  T0 = 0.0
  Tend = 730
  timestep=1
  # time span
  tspan = seq(T0, T0+Tend, by=timestep) 
  # National average oscillator
  if (region == 'Green'){ # Northeast
    # initial distribution
    IC = c(S=1,I=1000/71000000,R=0)
    # region-specific oscillator
  }else if(region == 'Blue'){ # South
    IC = c(S=1,I=1000/67000000,R=0)
  }else if (region == 'Red'){ # Midwest
    IC = c(S=1,I=1000/187000000,R=0)
  }

  # Fixed parameters
  beta0=0.185
  gamma= 0.16
  z = 0.005
  c.nao= 1
  a.nao= 0.4
  p.nao = 273
  f.nao = 1/58.1
  L = 1
  K = 0.0018
  # Varying parameters
  c.rso1.green = 1
  c.rso2.green = 0.9
  a.rso1.green = 0.2
  a.rso2.green = 0.1
  c.rso1.blue = 1.1
  c.rso2.blue = 1
  a.rso1.blue = 0.4
  a.rso2.blue = 0.3  
  c.rso1.red = 0.9
  c.rso2.red = 1
  a.rso1.red = 0.1
  a.rso2.red = 0.2
  p.rso = 45
  f.rso=1/29.05
  rso.matrix = matrix(c(c.rso1.green,c.rso2.green,a.rso1.green,a.rso2.green,p.rso,f.rso,
                        c.rso1.blue,c.rso2.blue,a.rso1.blue,a.rso2.blue,p.rso,f.rso,
                        c.rso1.red,c.rso2.red,a.rso1.red,a.rso2.red, p.rso,f.rso), byrow = TRUE, nrow=3, dimnames = list(c("Green","Blue","Red"),c("c1","c2","a1","a2","p","f")))
  # parameters to run SIRS model
  parameters = list(region = region, 
                    InfectionParam=c(beta0,gamma,z), 
                    NAOParam=c(c.nao,a.nao,p.nao,f.nao), 
                    RSOParam=rso.matrix,
                    VirusTransParam=c(L,K))
  # run ODE solver
  odetime <- proc.time()
  prob = rk(y=IC, times=tspan, func = ode_SIRS, parms = parameters, method='rk45dp7')
  proc.time() - odetime
  return(prob)
  } 
## A function to calculate error ##
cal_cost <- function(){
  
  prob <- run_SIRS()
  I.m <- # calculate 
  return(cost)
}

## Main ##
# check beta.d over time
t<-seq(0,Tend,by=1)
m.nao.1 = a.nao * sin((t-p.nao)/(1/f.nao))+c.nao 
j.d = L + K*t

cal_rso <- function(region){
  m.rso = rep(0,Tend+1)
  for (t in 0:Tend){
    if(region %in% c("Red","Green")){
      if ((t%/%183)%%2 == 0){ # alternating amplitude for every half-year
        c.rso = rso.matrix[region,1]
        a.rso = rso.matrix[region,3]
      }else{
        c.rso = rso.matrix[region,2]
        a.rso = rso.matrix[region,4]
      }
      m.rso[t]= a.rso * sin((t-p.rso)/(1/f.rso))+c.rso
  }else{
    if (t<=91 | ((t+91)%/%183)%%2 == 0){ # alternating amplitude for every half-year
      c.rso = rso.matrix[region,1]
      a.rso = rso.matrix[region,3]
    }else{
      c.rso = rso.matrix[region,2]
      a.rso = rso.matrix[region,4]
    }
    m.rso[t]= a.rso * sin((t+p.rso)/(1/f.rso))+c.rso}    
  }
  
  plot(m.rso)
    return(m.rso)
}
m.rso.gr <- cal_rso('Green')
m.rso.bl <- cal_rso('Blue')
m.rso.rd <- cal_rso('Red')

beta.d.gr = beta0 * ((m.nao.1+m.rso.gr)/2*j.d)
beta.d.bl = beta0 * ((m.nao.1+m.rso.bl)/2*j.d)
beta.d.rd = beta0 * ((m.nao.1+m.rso.rd)/2*j.d)
beta.all <- data.frame(time=seq(0,Tend),green=beta.d.gr,blue=beta.d.bl,red=beta.d.rd)
beta.all.t <- melt(beta.all,id.vars='time')
ggplot(data=beta.all.t,aes(time,value, color=variable))+
  geom_line()


# Test run with IC = c(1,0.01,0)
I.all <- data.frame(time=seq(0,Tend),green=prob.gr[,3],blue=prob.bl[,3],red=prob.rd[,3])
I.all.t <- melt(I.all,id.vars='time')
ggplot(data=I.all.t,aes(time,value, color=variable))+
  geom_line()+
  scale_color_manual(values=c("green","blue","red"))


# find optimal parameter sets 




