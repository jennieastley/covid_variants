# Setup ####
library(pacman)
p_load(deSolve,
       readxl,
       tidyverse,
       shiny,
       shinycssloaders,
       shinyBS)

# Parameters ####
R0 <- 0.3*5.08 # Reproductive number of V1
gamma1 <- 1/10 # Latent rate of V1
mu1 <- 0.0057 # Proportion of infected individuals that die from V1
beta1 <- R0 * gamma1 * (1 - mu1) # Transmission rate of V1
rho_init <- 0.01 # proportion vaccinated against V1 before V1 wave
rho <- 0.05 # Proportion of population vaccinated against V1 after V1 wave
i2 <- 0.001 # Proportion of population infected with V2 when it is first detected
v2_day <-365 # Day on which V2 first detected
phia <- 1.25 # Degree to which V2 is more transmissible than V1
phib <- 1 # Degree to which latent rate of V2 is higher than V1
phic <- 0.3 # Degree to which V2 is a V1 vaccine escape
phid <- 0.25 # Degree to which V2 is more virulent than V1
phie <- 0.5 # What proportion of individuals recovered from V1 are immune to V2
time_stop <- 730 # Total days in model

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

# Start app ####

# Define UI for application 
ui <- fluidPage(
  
  # Title
  titlePanel("SARS-CoV-2 new variant"),
  
  # Sidebar with slider inputs 
  sidebarLayout(
    sidebarPanel(
      actionButton("go", "Go"),
      
      bsCollapse(
        bsCollapsePanel("Multiplying factors: parameters of V2 compared to V1",
                        sliderInput(inputId="phia", label = "Transmission rate", value = 1, min=0.25, max=4,step=0.25),
                        sliderInput(inputId="phic", label = "Vaccine escape proportion", value = 0.25, min=0, max=1,step=0.05),
                        sliderInput(inputId="phid", label = "Virulence/Case Fatality Rate", value = 0.3, min=0.1, max=2.1,step=0.1),
                        sliderInput(inputId="phib", label = "Latent rate", value = 1, min=0.5, max=1.5,step=0.25),
                        sliderInput(inputId="phie", label = "Cross immunity", value = 0.1, min=0, max=1,step=0.05)
                        
        )
      ),
       
      bsCollapse(
        bsCollapsePanel("V1 transmission parameters",
                        sliderInput(inputId="rho_init", label = "Proportion vaccinated against V1 before initial wave", value = 0.01, min=0, max=1,step=0.05),
                        sliderInput(inputId="rho", label = "Proportion vaccinated against V1 after initial wave", value = 0.05, min=0, max=1,step=0.05),
                        sliderInput(inputId="mu1", label = "Case fatality rate of V1", value = 0.0057, min=0.01, max=0.1,step=0.001),
                        sliderInput(inputId="gamma1", label = "Latent rate of V1", value = 0.1, min=0.1, max=0.3,step=0.1),
                        sliderInput(inputId="beta1", label = "Transmission rate of V1", value = 0.152, min=0.1, max=0.8,step=0.05)
                        
        )
      )
    ),
    
    # Tabs for incidence and costs
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Daily incidence", plotOutput("incPlot") %>% withSpinner()),
                  tabPanel("Infections", plotOutput("infPlot") %>% withSpinner()),
                  tabPanel("Daily deaths", plotOutput("deathPlot") %>% withSpinner())
      )
    )
  )
)

# Server
server <- function(input, output) {
  modelOut <- eventReactive(input$go, {
    parameters["phia"] <- input$phia
    parameters["phib"] <- input$phib
    parameters["phic"] <- input$phic
    parameters["phid"] <- input$phid
    parameters["phie"] <- input$phie
    parameters["rho"] <- input$rho
    parameters["beta1"] <- input$beta1
    parameters["gamma1"] <- input$gamma1
    parameters["mu1"] <- input$mu1
    parameters["rho_init"] <- input$rho_init

    
    # Run the model
    
    # Initial conditions V1 
    
    S1_0 <- 1 - input$rho_init
    I1_0 <- 0.0001
    D1_0 <- 0
    R1_0 <- input$rho_init
    CInc1_0 <- 0
    init1 <- c(S1 = S1_0, I1 = I1_0, D1 = D1_0, R1 = R1_0, CInc1 = CInc1_0)
    
    # Define model V1 
    
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
    
    # Run model 1 
    
    out1 <- as.data.frame(ode(y=init1, times=tps1, func = covid_V1, parms = parameters))
    
    # Initial conditions V2 
    
    S2_0 <- 1 - input$rho*(1 - input$phic) - i2 - out1[v2_day,4] - input$phie*out1[v2_day,5]
    I2_0 <- i2
    D2_0 <- 0
    R2_0 <- input$rho*(1 - input$phic) + input$phie*out1[v2_day,5]
    CInc2_0 <- 0
    init2 <- c(S2 = S2_0, I2 = I2_0, D2 = D2_0, R2 = R2_0, CInc2 = CInc2_0)
    
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
    
  })
  
  #Define Plots
  
  output$deathPlot <- renderPlot({
    modelOut() %>% 
      filter(variable %in% c("Death1", "Death2")) %>% 
      group_by(variable) %>%
      ggplot()+
      geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
      theme_minimal() +
      labs(title = "Daily Deaths", x=('Time (days)'), y =("Proportion of population"), colour="Compartment")
  })
  
  output$incPlot <- renderPlot({
    modelOut() %>% 
      filter(variable %in% c("Inc1", "Inc2")) %>% 
      group_by(variable) %>%
      ggplot()+
      geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
      theme_minimal() +
      labs(title = "Daily Incidence", x=('Time (days)'), y =("Proportion of population"), colour="Compartment")# +
  })
  
  output$infPlot <- renderPlot({
    modelOut() %>% 
      filter(variable %in% c("I1", "I2")) %>% 
      group_by(variable) %>%
      ggplot()+
      geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
      theme_minimal() +
      labs(title = "Infections", x=('Time (days)'), y =("Proportion of population"), colour="Compartment")# +
  })
  

  
}

# Run the application 
shinyApp(ui = ui, server = server)
