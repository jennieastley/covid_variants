# Setup ####

library(deSolve)
library(tidyverse)
library(gridExtra)

# Parameters V1 ####
R0 <- 0.3*5.08 # Reproductive number of V1
gamma1 <- 1/10 # Latent rate of V1
mu1 <- 0.0057 # Proportion of infected individuals that die from V1
beta1 <- R0 * gamma1 * (1 - mu1) # Transmission rate of V1
rho_init <- 0.01 # proportion vaccinated against V1 before V1 wave
rho <- 0.05 # Proportion of population vaccinated against V1 after V1 wave
i2 <- 0.001 # Proportion of population infected with V2 when it is first detected
v2_day <-365 # Day on which V2 first detected
phia <- 1.3 # Degree to which V2 is more transmissible than V1
phib <- 1 # Degree to which latent rate of V2 is higher than V1
phic <- 0.325 # Degree to which V2 is a V1 vaccine escape
phid <- 0.3 # Degree to which V2 is more virulent than V1
phie <- 0.335 # What proportion of individuals recovered from V1 are immune to V2
time_stop <- 730 # Total days in model

col_infections <- "violet"
col_deaths <- "turquoise"

parameters <- c(
  beta1 <- beta1,
  gamma1 <- gamma1,
  phia <- phia,
  phib <- phib,
  phic <- phic,
  phid <- phid,
  phie <- phie,
  mu1 <- mu1,
  rho_init <- rho_init,
  rho <- rho,
  i2 <- i2,
  v2_day <- v2_day,
  time_stop <- time_stop
)

# Times V1 ####

time_start1 <- 0
time_stop1 <- time_stop
deltat1 <- 1
tps1 <- seq(time_start1, time_stop1, by=deltat1)


# Times V2 ####

time_start2 <- v2_day
time_stop2 <- time_stop
deltat2 <- 1
tps2 <- seq(time_start2, time_stop2, by=deltat2)

# Initial conditions V1 ####

S1_0 <- 1 - rho_init
I1_0 <- 0.0001
D1_0 <- 0
R1_0 <- rho_init
CInc1_0 <- 0
init1 <- c(S1 = S1_0, I1 = I1_0, D1 = D1_0, R1 = R1_0, CInc1 = CInc1_0)

# Define model V1 ####

covid_V1 <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    dS1 <- - beta1 * S1 * I1
    dI1 <- beta1 * S1 * I1 - gamma1 * (1 - mu1) * I1 - gamma1 * mu1 * I1
    dD1 <- mu1 * gamma1 * I1
    dR1 <- (1 - mu1) * gamma1 * I1
    dCInc1 <- beta1 * S1 * I1
    
    return(list(c(dS1,dI1,dD1,dR1,dCInc1)))
  }
  )
}



# Define model V2 ####

covid_V2 <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    
    beta2 <- phia * beta1
    gamma2 <- phib * gamma1
    mu2 <- phid * mu1
    
    dS2 <- - beta2 * S2 * I2
    dI2 <- beta2 * S2 * I2 - gamma2 * I2
    dD2 <- mu2 * gamma2 * I2
    dR2 <- (1 - mu2) * gamma2 * I2
    dCInc2 <- beta2 * S2 * I2
    
    return(list(c(dS2,dI2,dD2,dR2,dCInc2)))
  }
  )
}

# Run model 1 ####

out1 <- as.data.frame(ode(y=init1, times=tps1, func = covid_V1, parms = parameters))

# Initial conditions V2 ####

S2_0 <- 1 - rho*(1 - phic) - i2 - out1[v2_day,4] - phie*out1[v2_day,5]
I2_0 <- i2
D2_0 <- 0
R2_0 <- rho*(1 - phic) + phie*out1[v2_day,5]
CInc2_0 <- 0
init2 <- c(S2 = S2_0, I2 = I2_0, D2 = D2_0, R2 = R2_0, CInc2 = CInc2_0)


# Run model 2 and manipulate data ####

out2_growth <- ode(y=init2, times=tps2, func = covid_V2, parms = parameters)

out2_initial1 <- matrix(0:(v2_day-1), nrow=v2_day, ncol=1)
colnames(out2_initial1) <- c("time")

out2_initial2 <- matrix(rep((1 - rho*(1 - phic) - i2 - out1[time_stop,4] - phie*out1[time_stop,5]),v2_day), nrow=v2_day, ncol=1)
colnames(out2_initial2) <- c("S2")

out2_initial3 <- matrix(rep(0,v2_day), nrow = v2_day, ncol=1)
colnames(out2_initial3) <- c("I2")

out2_initial4 <- matrix(rep(0,v2_day), nrow = v2_day, ncol=1)
colnames(out2_initial4) <- c("D2")

out2_initial5 <- matrix(rep((rho*(1 - phic) + i2 + phie*out1[time_stop,5]),v2_day), nrow=v2_day, ncol=1)
colnames(out2_initial5) <- c("R2")

out2_initial6 <- matrix(rep(0,v2_day), nrow = v2_day, ncol=1)
colnames(out2_initial6) <- c("CInc2")

out2_initial <- cbind(out2_initial1,out2_initial2,out2_initial3,out2_initial4,out2_initial5,out2_initial6)
tail(out2_initial)

out2_combined <- rbind(out2_initial,out2_growth)

out2 <- as.data.frame(out2_combined)

out <- merge(out1,out2,by="time")

output <- as_tibble(out) %>% 
  mutate(P1 = S1+I1+D1+R1,
         P2 = S2+I2+D2+R2,
         Inc1 = c(0, diff(CInc1)),
         Inc2 = c(0, diff(CInc2)),
         Death1 = c(0, diff(D1)),
         Death2 = c(0, diff(D2))) %>% 
  pivot_longer(names_to = "variable", cols = !1) %>% 
  mutate(variant = ifelse(str_ends(variable, "2"), "V2", "V1")
  )

# V1 plots ####

plot1 <- output %>% 
  filter(variable %in% c("Inc1"), time<v2_day,) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)),colour=col_infections)+
  theme_minimal() +
  labs(title = "Daily Incidence of V1 Wave", x=('Time (days)'), y =("Proportion of population"), colour="Compartment") +
  facet_wrap(~variable)

plot2 <- output %>% 
  filter(variable %in% c("Death1"), time<v2_day,) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)),colour=col_deaths)+
  theme_minimal() +
  labs(title = "Daily Deaths of V1 Wave", x=('Time (days)'), y =("Proportion of population"), colour="Compartment") +
  facet_wrap(~variable)

grid.arrange(plot1, plot2, ncol=2)

# Other compartment plots ####

output %>% 
  filter(variable %in% c("S1", "I1", "D1", "R1", "S2", "I2", "D2", "R2")) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "SIR Compartments", x=('Time (days)'), y =("Proportion of population"), colour="Compartment") +
  facet_wrap(~variant)

output %>% 
  filter(variable %in% c("Inc1", "Inc2")) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Daily Incidence", x=('Time (days)'), y =("Proportion of population"), colour="Compartment") +
  facet_wrap(~variant)

output %>% 
  filter(variable %in% c("Death1", "Death2")) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Daily Deaths", x=('Time (days)'), y =("Proportion of population"), colour="Compartment") +
  facet_wrap(~variant)

output %>% 
  filter(variable %in% c("Inc1", "Inc2")) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
    labs(title = "Daily Incidence", x=('Time (days)'), y =("Proportion of population"), colour="Compartment")# +

output %>% 
  filter(variable %in% c("Death1", "Death2")) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Daily Deaths", x=('Time (days)'), y =("Proportion of population"), colour="Compartment")


# Gathering results ####

output %>% 
  filter(variable %in% c("I2")) %>% 
  group_by(variable) %>%
  ggplot()+
  xlim(v2_day,time_stop) +
  geom_vline(xintercept=v2_day+107) +
  geom_hline(yintercept=0.0507) +
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Infections", x=('Time (days)'), y =("Proportion of population"), colour="Compartment")

output %>% 
  filter(variable %in% c("Inc2")) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_vline(xintercept=v2_day+98) +
  geom_hline(yintercept=0.00533) +
  xlim(v2_day,time_stop) +
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Daily Incidence", x=('Time (days)'), y =("Proportion of population"), colour="Compartment")

tail(out2_growth)

# *BEWARE - CAN BE LONG TO RUN* Sensitivity analysis - all five parameters on INFECTIONS ####

phia_vector<-seq(0.75,3.4,by=0.6625) 
phib_vector<-seq(1,1,by=1) 
phic_vector<-seq(0.5,0.875,by=0.09375) 
phid_vector<-seq(0.1,2,by=0.475) 
phie_vector<-seq(0.17,0.5,by=0.0825) 

result3<-matrix(0,nrow = 1,ncol = 7)
colnames(result3)<-c("time","I","phic","phia","phid","phib","phie")
result3<-as.data.frame(result3)
for (m in 1:length(phie_vector)){
  for(l in 1:length(phib_vector)){
    for (k in 1:length(phid_vector)){
      for (j in 1:length(phia_vector)){
        for (i in 1:length(phic_vector)){
          parameters["phic"]<-phic_vector[i]
          parameters["phia"]<-phia_vector[j]
          parameters["phid"]<-phid_vector[k]
          parameters["phib"]<-phib_vector[l]
          parameters["phie"]<-phie_vector[m]
          
          S2_0 <- 1 - rho*(1 - phic_vector[i]) - i2 - out1[v2_day,4] - phie_vector[m]*out1[v2_day,5]
          I2_0 <- i2
          D2_0 <- 0
          R2_0 <- rho*(1 - phic_vector[i]) + phie_vector[m]*out1[v2_day,5]
          CInc2_0 <- 0
          init2 <- c(S2 = S2_0, I2 = I2_0, D2 = D2_0, R2 = R2_0, CInc2 = CInc2_0)
          
          out_aux <- ode(y = init2, times = tps2, func = covid_V2, parms = parameters) 
          
          aux_mat<-cbind(out_aux[,c(1,3)],rep(phic_vector[i],length(tps2)), rep(phic_vector[j],length(tps2)), rep(phic_vector[k],length(tps2)), rep(phic_vector[l],length(tps2)), rep(phic_vector[m],length(tps2)))

          colnames(aux_mat)<-c("time","I","phic","phia","phid","phib","phie")
          result3<-rbind(result3,aux_mat)
        }
      }
    }
  }
}

result3 <- pivot_longer(result3, names_to="phi", cols=3:7 )



p13 <- out1 %>% 
  ggplot() +
  geom_line(aes(x=time,y=I1),colour=col_infections) +
  xlim(0,v2_day) + 
  ylim(0,0.4) +
  theme_minimal() +
  labs(title = "V1 Infections", x=('Time (days)'), y="Proportion of population")

p14 <- result3 %>% 
  ggplot() +
  geom_line(aes(x=time,y=I),colour=col_infections) +
  xlim(v2_day,time_stop) +
  ylim(0,0.4) +
  theme_minimal() +
  labs(title = "V2 Infections", x=('Time (days)'), y =("Proportion of population"))

grid.arrange(p13, p14,  nrow=2)

# *BEWARE - CAN BE LONG TO RUN* Sensitivity analysis - all five parameters on CUMULATIVE DEATHS ####

phia_vector<-seq(0.75,3.4,by=0.6625) 
phib_vector<-seq(1,1,by=1) 
phic_vector<-seq(0.5,0.875,by=0.09375) 
phid_vector<-seq(0.1,2,by=0.475) 
phie_vector<-seq(0.17,0.5,by=0.0825) 

result4<-matrix(0,nrow = 1,ncol = 7)
colnames(result4)<-c("time","D","phic","phia","phid","phib","phie")
result4<-as.data.frame(result4)
for (m in 1:length(phie_vector)){
  for(l in 1:length(phib_vector)){
    for (k in 1:length(phid_vector)){
      for (j in 1:length(phia_vector)){
        for (i in 1:length(phic_vector)){
          parameters["phic"]<-phic_vector[i]
          parameters["phia"]<-phia_vector[j]
          parameters["phid"]<-phid_vector[k]
          parameters["phib"]<-phib_vector[l]
          parameters["phie"]<-phie_vector[m]
          
          S2_0 <- 1 - rho*(1 - phic_vector[i]) - i2 - out1[v2_day,4] - phie_vector[m]*out1[v2_day,5]
          I2_0 <- i2
          D2_0 <- 0
          R2_0 <- rho*(1 - phic_vector[i]) + phie_vector[m]*out1[v2_day,5]
          CInc2_0 <- 0
          init2 <- c(S2 = S2_0, I2 = I2_0, D2 = D2_0, R2 = R2_0, CInc2 = CInc2_0)
          
          out_aux <- ode(y = init2, times = tps2, func = covid_V2, parms = parameters) 
          
          aux_mat<-cbind(out_aux[,c(1,4)],rep(phic_vector[i],length(tps2)), rep(phic_vector[j],length(tps2)), rep(phic_vector[k],length(tps2)), rep(phic_vector[l],length(tps2)), rep(phic_vector[m],length(tps2)))
          
          colnames(aux_mat)<-c("time","D","phic","phia","phid","phib","phie")
          result4<-rbind(result4,aux_mat)
        }
      }
    }
  }
}

result4 <- pivot_longer(result4, names_to="phi", cols=3:7 )

p15 <- out1 %>% 
  ggplot() +
  geom_line(aes(x=time,y=D1),colour=col_deaths) +
  xlim(0,v2_day) + 
  ylim(0,0.0125) +
  theme_minimal() +
  labs(title = "V1 Cumulative Deaths", x=('Time (days)'), y="Proportion of population")

p16 <- result4 %>% 
  ggplot() +
  geom_line(aes(x=time,y=D),colour=col_deaths) +
  xlim(v2_day,time_stop) +
  ylim(0,0.0125) +
  theme_minimal() +
  labs(title = "V2 Cumulative Deaths", x=('Time (days)'), y =("Proportion of population"))

grid.arrange(p15, p16,  nrow=2)

# Sensitivity analysis  - vaccine escape and cross immunity on cumulative INFECTIONS ####

parameters["phia"] <- 1.25 # Transmission rate
parameters["phib"] <- 1 # Latent rate
phic_vector<-seq(0.5,0.875,by=0.09375) # Vaccine escape
parameters["phid"] <- 0.25 # Virulence
phie_vector<-seq(0.17,0.5,by=0.0825) # Cross immunity

cumulative_infections_result<-matrix(0,nrow = length(phic_vector),length(phie_vector)) 
for (i in 1:length(phic_vector)){
  for (j in 1:length(phie_vector)){
    parameters["phic"]<-phic_vector[i] # Virulence mu
    parameters["phie"]<-phie_vector[j] # Transmissibility beta
    
    S2_0 <- 1 - rho*(1 - phic_vector[i]) - i2 - out1[v2_day,4] - phie_vector[j]*out1[v2_day,5]
    I2_0 <- i2
    D2_0 <- 0
    R2_0 <- rho*(1 - phic_vector[i]) + phie_vector[j]*out1[v2_day,5]
    CInc2_0 <- 0
    init2 <- c(S2 = S2_0, I2 = I2_0, D2 = D2_0, R2 = R2_0, CInc2 = CInc2_0)
    
    out <- ode(y = init2, times = tps2, func = covid_V2, parms = parameters) 
    
    cumulative_infections_result[i,j]<-tail(out[,6],1) 
  }
}

#contour plot of a matrix. don't want to draw the x or y axis. give labels. 10 levels:
contour(cumulative_infections_result,xaxt = "n", yaxt = "n",ylab="Cross Immunity phi_e",xlab="Vaccine Escape phi_c",nlevels = 10,labcex=1.3,col = hcl.colors(11, "Temps"),lwd =2)
axis(1, at=1:length(phic_vector)/length(phic_vector), labels=phic_vector)
axis(2, at=1:length(phie_vector)/length(phie_vector), labels=phie_vector)
title("Cumulative Infections (proportion of population)")


# Sensitivity analysis - worst case with beta on INFECTIONS ####

phia_vector<-seq(0.75,3.4,by=0.53) 
parameters["phib"] <- 1 # Latent rate
parameters["phic"] <- 0.875 # Vaccine escape
parameters["phid"] <- 2 # Virulence
parameters["phie"] <- 0.17 # Cross immunity

result_phia<-matrix(0,nrow = 1,ncol = 3)
colnames(result_phia)<-c("time","I","phia")
result_phia<-as.data.frame(result_phia)
for (i in 1:length(phia_vector)){
  parameters["phia"]<-phia_vector[i]
  out_aux <- ode(y = init2, times = tps2, func = covid_V2, parms = parameters) 
  
  aux_mat<-cbind(out_aux[,c(1,3)],rep(phia_vector[i],length(tps2)))
  colnames(aux_mat)<-c("time","I","phia")
  result_phia<-rbind(result_phia,aux_mat)
}

p1 <- out1 %>% 
  ggplot() +
  geom_line(aes(x=time,y=I1)) +
  xlim(0,v2_day) + 
  ylim(0,0.25) +
  theme_minimal() +
  labs(title = "V1 Infections", x=('Time (days)'), y="Proportion of population")

p2 <- result_phia %>% 
  ggplot() +
  geom_line(aes(x=time,y=I,colour=as.factor(phia))) +
  xlim(v2_day,time_stop) +
  ylim(0,0.25) +
  theme_minimal() +
  labs(title = "V2 Infections", x=('Time (days)'), y =(" "), colour="Beta Multiplying Factor")

grid.arrange(p1, p2, widths=c(0.38, 0.62), ncol=2)

# Sensitivity analysis - cross immunity on INFECTIONS ####

parameters["phia"] <- 1.25 # Transmissibility
parameters["phib"] <- 1 # Latent rate
parameters["phic"] <- 0.6 # Vaccine escape
parameters["phid"] <- 0.25 # Virulence
phie_vector<-seq(0.17,0.5,by=0.0825) 

result_phie<-matrix(0,nrow = 1,ncol = 3)
colnames(result_phie)<-c("time","I","phie")
result_phie<-as.data.frame(result_phie)
for (i in 1:length(phie_vector)){
  parameters["phie"]<-phie_vector[i]
  
  out1 <- as.data.frame(ode(y=init1, times=tps1, func = covid_V1, parms = parameters))
  
  S2_0 <- 1 - rho*(1 - phic) - i2 - out1[v2_day,4] - phie_vector[i]*out1[v2_day,5]
  I2_0 <- i2
  D2_0 <- 0
  R2_0 <- rho*(1 - phic) + phie_vector[i]*out1[v2_day,5]
  CInc2_0 <- 0
  init2 <- c(S2 = S2_0, I2 = I2_0, D2 = D2_0, R2 = R2_0, CInc2 = CInc2_0)
  
  out_aux <- ode(y = init2, times = tps2, func = covid_V2, parms = parameters) 
  
  aux_mat<-cbind(out_aux[,c(1,3)],rep(phie_vector[i],length(tps2)))
  colnames(aux_mat)<-c("time","I","phie")
  result_phie<-rbind(result_phie,aux_mat)
}

p7 <- out1 %>% 
  ggplot() +
  geom_line(aes(x=time,y=I1)) +
  xlim(0,v2_day) + 
  ylim(0,0.1) +
  theme_minimal() +
  labs(title = "V1 Infections", x=('Time (days)'), y="Proportion of population")

p8 <- result_phie %>% 
  ggplot() +
  geom_line(aes(x=time,y=I,colour=as.factor(phie))) +
  xlim(v2_day,time_stop) +
  ylim(0,0.1) +
  theme_minimal() +
  labs(title = "V2 Infections", x=('Time (days)'), y =(" "), colour="Cross Immunity")

grid.arrange(p7, p8, widths=c(0.4,0.6), ncol=2)


# 1 MONTH LATER - Sensitivity analysis - all three parameters on INFECTIONS ####

phia_vector<-seq(1.1,1.5,by=0.1) # Beta
parameters["phib"] <- 1 # Latent rate
phic_vector<-seq(0.25,0.4,by=0.05) # Vaccine escape
phid_vector<-seq(0.1,0.5,by=0.1) # Virulence
parameters["phie"] <- 0.3 # Cross immunity


result1<-matrix(0,nrow = 1,ncol = 5)
colnames(result1)<-c("time","I","phic","phia","phid")
result1<-as.data.frame(result1)
for (k in 1:length(phid_vector)){
 for (j in 1:length(phia_vector)){
  for (i in 1:length(phic_vector)){
    parameters["phic"]<-phic_vector[i]
    parameters["phia"]<-phia_vector[j]
    parameters["phid"]<-phid_vector[k]
    
    
    S2_0 <- 1 - rho*(1 - phic_vector[i]) - i2 - out1[v2_day,4] - phie*out1[v2_day,5]
    I2_0 <- i2
    D2_0 <- 0
    R2_0 <- rho*(1 - phic_vector[i]) + phie*out1[v2_day,5]
    CInc2_0 <- 0
    init2 <- c(S2 = S2_0, I2 = I2_0, D2 = D2_0, R2 = R2_0, CInc2 = CInc2_0)
    
    out_aux <- ode(y = init2, times = tps2, func = covid_V2, parms = parameters) 

    aux_mat<-cbind(out_aux[,c(1,3)],rep(phic_vector[i],length(tps2)), rep(phic_vector[j],length(tps2)), rep(phic_vector[k],length(tps2)))

    colnames(aux_mat)<-c("time","I","phic","phia","phid")
    result1<-rbind(result1,aux_mat)
  }
 }
}
result1 <- pivot_longer(result1, names_to="phi", cols=3:5 )


p9 <- out1 %>% 
  ggplot() +
  geom_line(aes(x=time,y=I1),colour=col_infections) +
  xlim(0,v2_day) + 
  ylim(0,0.07) +
  theme_minimal() +
  labs(title = "V1 Infections", x=('Time (days)'), y="Proportion of population")

p10 <- result1 %>% 
  ggplot() +
  geom_line(aes(x=time,y=I),colour=col_infections) +
  xlim(v2_day,time_stop) +
  ylim(0,0.07) +
  theme_minimal() +
  labs(title = "V2 Infections", x=('Time (days)'), y =("Proportion of population"))

grid.arrange(p9, p10, nrow=2)


# 1 MONTH LATER - Sensitivity analysis - all three parameters on CUMULATIVE DEATHS ####

phia_vector<-seq(1.1,1.5,by=0.1) # Beta
parameters["phib"] <- 1 # Latent rate
phic_vector<-seq(0.25,0.4,by=0.05) # Vaccine escape
phid_vector<-seq(0.1,0.5,by=0.1) # Virulence
parameters["phie"] <- 0.3 # Cross immunity


result2<-matrix(0,nrow = 1,ncol = 5)
colnames(result2)<-c("time","D","phic","phia","phid")
result2<-as.data.frame(result2)
for (k in 1:length(phid_vector)){
  for (j in 1:length(phia_vector)){
    for (i in 1:length(phic_vector)){
      parameters["phic"]<-phic_vector[i]
      parameters["phia"]<-phia_vector[j]
      parameters["phid"]<-phid_vector[k]
     
      S2_0 <- 1 - rho*(1 - phic_vector[i]) - i2 - out1[v2_day,4] - phie*out1[v2_day,5]
      I2_0 <- i2
      D2_0 <- 0
      R2_0 <- rho*(1 - phic_vector[i]) + phie*out1[v2_day,5]
      CInc2_0 <- 0
      init2 <- c(S2 = S2_0, I2 = I2_0, D2 = D2_0, R2 = R2_0, CInc2 = CInc2_0)
      
      out_aux <- ode(y = init2, times = tps2, func = covid_V2, parms = parameters) 
      
      aux_mat<-cbind(out_aux[,c(1,4)],rep(phic_vector[i],length(tps2)), rep(phic_vector[j],length(tps2)), rep(phic_vector[k],length(tps2)))

      colnames(aux_mat)<-c("time","D","phic","phia","phid")
      result2<-rbind(result2,aux_mat)
    }
  }
}
result2 <- pivot_longer(result2, names_to="phi", cols=3:5 )


p11 <- out1 %>% 
  ggplot() +
  geom_line(aes(x=time,y=D1),colour=col_deaths) +
  xlim(0,v2_day) + 
  ylim(0,0.004) +
  theme_minimal() +
  labs(title = "V1 Cumulative Deaths", x=('Time (days)'), y="Proportion of population")

p12 <- result2 %>% 
  ggplot() +
  geom_line(aes(x=time,y=D),colour=col_deaths) +
  xlim(v2_day,time_stop) +
  ylim(0,0.004) +
  theme_minimal() +
  labs(title = "V2 Cumulative Deaths", x=('Time (days)'), y =("Proportion of population"))

grid.arrange(p11, p12, nrow=2)

# 1 MONTH LATER - Sensitivity analysis - cross immunity on INFECTIONS ####

parameters["phia"] <- 1.3 # Transmissibility
parameters["phib"] <- 1 # Latent rate
parameters["phic"] <- 0.325 # Vaccine escape
parameters["phid"] <- 0.2 # Virulence
phie_vector<-seq(0.17,0.5,by=0.0825) 

result_phie<-matrix(0,nrow = 1,ncol = 3)
colnames(result_phie)<-c("time","I","phie")
result_phie<-as.data.frame(result_phie)
for (i in 1:length(phie_vector)){
  parameters["phie"]<-phie_vector[i]
  
  out1 <- as.data.frame(ode(y=init1, times=tps1, func = covid_V1, parms = parameters))
  
  S2_0 <- 1 - rho*(1 - phic) - i2 - out1[v2_day,4] - phie_vector[i]*out1[v2_day,5]
  I2_0 <- i2
  D2_0 <- 0
  R2_0 <- rho*(1 - phic) + phie_vector[i]*out1[v2_day,5]
  CInc2_0 <- 0
  init2 <- c(S2 = S2_0, I2 = I2_0, D2 = D2_0, R2 = R2_0, CInc2 = CInc2_0)
  
  out_aux <- ode(y = init2, times = tps2, func = covid_V2, parms = parameters) 

  aux_mat<-cbind(out_aux[,c(1,3)],rep(phie_vector[i],length(tps2)))
  colnames(aux_mat)<-c("time","I","phie")
  result_phie<-rbind(result_phie,aux_mat)
}

p17 <- out1 %>% 
  ggplot() +
  geom_line(aes(x=time,y=I1)) +
  xlim(0,v2_day) + 
  ylim(0,0.1) +
  theme_minimal() +
  labs(title = "V1 Infections", x=('Time (days)'), y="Proportion of population")

p18 <- result_phie %>% 
  ggplot() +
  geom_line(aes(x=time,y=I,colour=as.factor(phie))) +
  xlim(v2_day,time_stop) +
  ylim(0,0.1) +
  theme_minimal() +
  labs(title = "V2 Infections", x=('Time (days)'), y =(" "), colour="Cross Immunity")

grid.arrange(p17, p18, widths=c(0.4,0.6), ncol=2)