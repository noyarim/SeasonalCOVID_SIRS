library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)


##########################
###### Functions #########
##########################

## Differential equation for SIR-R model of COVID-19 ##
ode_SIRS <- function(t, state, parameters){
  with(parameters,{
    #region = parameters$region
    InfectionParam = parameters$InfectionParam
    #VirusTransParam = parameters$VirusTransParam
    # Disaggregate fixed parameters
    beta0 = InfectionParam[1]
    g = InfectionParam[2]
    z = InfectionParam[3]
    # weekly beta multiplier
    m.beta = parameters$m.beta
    # calculate daily beta multiplier
    if(t ==0){
      w_idx = 1
    }else{
      w_idx = floor(t/n_days_per_unit)
      }
    # find week index
    m.beta.d <- m.beta[w_idx]

    beta.d = beta0 * m.beta.d
    # ODE equations
    S=state[1]
    I=state[2]
    R=state[3]
    C=state[4]
    dS = -beta.d*S*I + z*R
    dI = beta.d*S*I - g*I
    dR = g*I - z*R
    dC = beta.d*S*I
    
    
    return(list(c(dS,dI,dR,dC)))
  })
}

## A function to run SIR-R model of COVID-19 ##
run_SIRS <- function(parameters){
  with(parameters, {
    timestep=1
    # time span
    tspan = seq(T0, T0+Tend, by=timestep) 
    # # Initial distribution by region
    # if (region == 'Green'){ # Northeast
    #   # initial distribution
    #   IC = c(S=1-1000/71000000,I=1000/71000000,R=0,C=0)
    #   # region-specific oscillator
    # }else if(region == 'Blue'){ # South
    #   IC = c(S=1-1000/67000000,I=1000/67000000,R=0,C=0)
    # }else if (region == 'Red'){ # Midwest
    #   IC = c(S=1-1000/187000000,I=1000/187000000,R=0,C=0)
    # }

    # run ODE solver
    odetime <- proc.time()
    prob = rk(y=IC, times=tspan, func = ode_SIRS, parms = parameters, method='rk45dp7')
    #prob = rk(y = IC, times = tspan, func = ode_SIRS, parms = parameters, method='rk45dp7',atol = 10^-8, rtol = 10^-10)  # TODO: constrain domain for x to be 0 or positive (x>0)
    proc.time() - odetime
    return(prob)
  })
} 

## A function to calculate the number of new cases ##
cal_case <- function(prob){
  new_case <-(prob[-1,5]- prob[-nrow(prob),5])
  new_case[1] = 0

  return(new_case)
  

}

## A function to calculate case rate per 100,000 ##
cal_caserate <- function(prob, scale){
  new_case <- cal_case(prob)
  case_rate <- new_case / prob[-1,1] * 10^5
  
  if (scale %in% c('log','Log')){
    case_rate_log <- log(new_case)
    case_rate_log[is.infinite(case_rate_log)] = 0
    
    return(new_case_log)
  }else{
    return(case_rate)
}
}

## A function to calculate error ##
cal_cost <- function(m.beta, parameters, data, grouping, group){
  with(parameters, {
    parameters$m.beta = m.beta
    # model outcome
    prob <- run_SIRS(parameters)
    I.m <- cal_caserate(prob, 'natural')
    # observed case 
    temp <- data %>%
      filter(!!sym(grouping) == group) %>%
      group_by(date)%>%
      summarize(
        w.case_rate = weighted.mean(case_rate_sm, pop2020)
      ) %>%
      mutate(
        w.case_rate_log = log(w.case_rate),
        w.case_rate_log = ifelse(is.infinite(w.case_rate_log),0,w.case_rate_log)
      )
    
    I.o <- temp$w.case_rate
    #print(c(length(m.beta[m.beta<0]),length(I.m[is.na(I.m)]), length(I.o[is.na(I.o)])))
    
    # sum of squared errors
    cost <- sum((I.m - I.o)^2)
    #print(c(m.beta,cost))
    if (is.nan(cost)){
      cost = 10^10
    }
    #print(cost)
    return(cost)
  })
}

## A function to run SIR with the calibrated betas ##
run_SIRS_clbr <- function(m.beta, parameters){
  parameters$m.beta = m.beta
  prob = run_SIRS(parameters)
  I.m <- cal_caserate(prob,'natural')
  
  return(list(prob,I.m))
}



#####################
###### Data #########
#####################
# covid data by state
covid_by_state <- read.csv("data/covid_by_state_cln.csv")
# Try different grouping of states
covid_by_state2 <- covid_by_state %>%
  filter(date <= as.Date("2021-03-01"))%>%
  mutate(
    # group 1: divide the states into 4: NE, MW, S, W
    grp1 = ifelse(state_ab %in% c("VT","CT","PA","ME","MA","NH","NJ","NY","RI"),'NE',
                  ifelse(state_ab %in% c("MN","ND","SD","IA","NE","WI","MO","IN","IL","MI","OH","KS"),'MW',
                         ifelse(state_ab %in% c("MD","DE","WV","VA","TX","TN","SC","OK","NC","MS","LA","KY","GA","FL","AR","AL"),'S',
                                ifelse(state_ab %in% c("ID", "MT","CO","WY","UT","NV","AZ","CA","WA","OR","NM"),'W',NA)))),
    # group 2: divide the states into 2: N, S
    grp2 = ifelse(grp1 == 'S', 'S', 'N' )
  )

test <- covid_by_state2 %>%
  group_by(grp1, date) %>%
  summarize(sum_case_sm = sum(cases_sm),
            sum_pop = sum(pop2020)) %>%
  filter(date == "2020-03-01")%>%
  mutate(IC = sum_case_sm/sum_pop) 

############################
###### Calibration #########
###########################

# Model set up parameters
T0 = 0.0
Tend = length(levels(as.factor(covid_by_state2$date)))
n_days = as.integer(Tend-T0)
n_days_per_unit = 30.5
# Initial state distribution
IC = c(S=1-1000/7100000,I=1000/7100000,R=0,C=0) # Northeast
#IC = c(S=1-0.00000000767,I=0.00000000767,R=0,C=0) # Northeast

IC = c(S=1-1000/67000000,I=1000/67000000,R=0,C=0) # South
#IC = c(S=1-1000/187000000,I=1000/187000000,R=0,C=0) # Midwestt
# Fixed infection parameters
beta0=1#0.185
gamma= 0.16
z = 0.005

# varying parameters (weekly beta multiplier)
m.beta = rep(1,as.integer(n_days/n_days_per_unit))

# parameters to run SIRS model
parameters = list(
                  T0 = T0,
                  Tend = Tend,
                  n_days = n_days,
                  n_days_per_unit = n_days_per_unit,
                  IC = IC,
                  InfectionParam=c(beta0,gamma,z)
                  )

# Directed search of weekly beta values by group  
grp1_NE<-optim(m.beta, fn = cal_cost, parameters = parameters, data=covid_by_state2, grouping = 'grp1', group='NE', method="BFGS")
write.csv(grp1_NE$par, "results/grp1_NE_betas(monthly)2.csv")
grp1_S<-optim(m.beta, fn = cal_cost, parameters = parameters, data=covid_by_state2, grouping = 'grp1', group='S', method="BFGS")
write.csv(grp1_S$par, "results/grp1_S_betas(monthly)2.csv")
grp1_W<-optim(m.beta, fn = cal_cost, parameters = parameters, data=covid_by_state2, grouping = 'grp1', group='W', method="BFGS")
write.csv(grp1_W$par, "results/grp1_W_betas(monthly)2.csv")
grp1_MW<-optim(m.beta, fn = cal_cost, parameters = parameters, data=covid_by_state2, grouping = 'grp1', group='MW', method='BFGS')
write.csv(grp1_MW$par, "results/grp1_MW_betas(monthly)2.csv")


#########################
####### Plotting ########
#########################
# observed case rate
dt_grp1_o <- covid_by_state2 %>%
  group_by(grp1, date) %>%
  summarize(w_case_rate = weighted.mean(case_rate_sm, pop2020))%>%
  mutate(type='observed')
# calibrated betas
grp1_NE_betas <- read.csv("results/grp1_NE_betas.csv")$x
grp1_S_betas <- read.csv("results/grp1_S_betas.csv")$x
grp1_W_betas <- read.csv("results/grp1_W_betas.csv")$x
grp1_MW_betas <- read.csv("results/grp1_MW_betas.csv")$x
# modeled case rate
case_rate_grp1_NE = run_SIRS_clbr(grp1_NE_betas,parameters)[[2]]
case_rate_grp1_S = run_SIRS_clbr(grp1_S_betas,parameters)[[2]]
case_rate_grp1_W = run_SIRS_clbr(grp1_W_betas,parameters)[[2]]
case_rate_grp1_MW = run_SIRS_clbr(grp1_MW_betas,parameters)[[2]]
dt_grp1_m <- data.frame(
  date = dt_grp1_o$date,
  MW = case_rate_grp1_MW/100,
  NE = case_rate_grp1_NE/100,
  S = case_rate_grp1_S/100,
  W = case_rate_grp1_W/100
)
dt_grp1_m_t <- melt(dt_grp1_m, id.vars="date")
dt_grp1_m_t <- dt_grp1_m_t %>%
  rename(grp1 = variable, w_case_rate=value)%>%
  select(grp1,date,w_case_rate)%>%
  mutate(type = 'modeled')
dt_grp1_m_t$type = 'modeled'
dt_grp1_all <- rbind(dt_grp1_o, dt_grp1_m_t)
dt_grp1_all$grp1 = as.factor(dt_grp1_all$grp1)
# Plot the observed vs. modeled case rate
ggplot(dt_grp1_all, aes(date,w_case_rate,group=type,color=type)) +
  geom_line()+
  facet_wrap(~grp1)+
  theme_bw()
ggsave("graphs/calibration_fit(monthly).pdf")

# Test run with IC = c(1,0.01,0)
I.all <- data.frame(time=seq(0,Tend),green=prob.gr[,3],blue=prob.bl[,3],red=prob.rd[,3])
I.all.t <- melt(I.all,id.vars='time')
ggplot(data=I.all.t,aes(time,value, color=variable))+
  geom_line()+
  scale_color_manual(values=c("green","blue","red"))

