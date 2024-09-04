library(shiny)
library(ggplot2)
library(tidyverse)


max_value <- 3 # Max value for contribution sliders.

# Define UI
ui <- fluidPage(
  # Application title
  titlePanel("Tripartite Model Simulation"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    sidebarPanel(
      sliderInput("tmax", "Time (tmax):", min = 50, max = 500, value = 100, step = 10),
      # sliderInput("k", "k:", min = 1, max = 20, value = 10, step = .1),
      # sliderInput("b", "b:", min = 50, max = 200, value = 100, step = 1),
      sliderInput("dBV", "Bacterial absorbtion rate of viruses (dBV):", min = 0.1, max = 1, value = 0.2, step = 0.1),
      sliderInput("dHB", "Host absorbtion rate of bacteria (dHB):", min = 0.1, max = 1, value = 0.3, step = 0.1),
      sliderInput("initial_H", "Initial Hosts:", min = 0, max = 100, value = 50, step = 5),
      sliderInput("initial_B", "Initial Bacteria:", min = 0, max = 100, value = 50, step = 5),
      sliderInput("initial_V", "Initial Virus:", min = 0, max = 100, value = 50, step = 5),
    ),
    
    mainPanel(
      plotOutput("distPlot"),
      
      # Sliders for FH, FB, FV below the plot
      h3("Contribution to Pools"),
      fluidRow(
        column(3, sliderInput("FH1", "FH[1]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FH2", "FH[2]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FH3", "FH[3]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FH4", "FH[4]:", min = 0, max = max_value, value = 1, step = 0.1))
      ),
      fluidRow(
        column(3, sliderInput("FB1", "FB[1]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FB2", "FB[2]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FB3", "FB[3]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FB4", "FB[4]:", min = 0, max = max_value, value = 1, step = 0.1))
      ),
      fluidRow(
        column(3, sliderInput("FV1", "FV[1]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FV2", "FV[2]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FV3", "FV[3]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FV4", "FV[4]:", min = 0, max = max_value, value = 1, step = 0.1))
      ),
      
      # Randomize buttons
      fluidRow(
        column(3, actionButton("randomize_FH", "Randomize FH")),
        column(3, actionButton("randomize_FB", "Randomize FB")),
        column(3, actionButton("randomize_FV", "Randomize FV")),
        column(3, actionButton("randomize_all", "Randomize All"))
      )
    )
  )
)




# Define server logic required to draw the plot
server <- function(input, output, session) {
  
  observeEvent(input$randomize_FH, {
    updateSliderInput(session, "FH1", value = runif(1, 0, max_value))
    updateSliderInput(session, "FH2", value = runif(1, 0, max_value))
    updateSliderInput(session, "FH3", value = runif(1, 0, max_value))
    updateSliderInput(session, "FH4", value = runif(1, 0, max_value))
  })
  
  observeEvent(input$randomize_FB, {
    updateSliderInput(session, "FB1", value = runif(1, 0, max_value))
    updateSliderInput(session, "FB2", value = runif(1, 0, max_value))
    updateSliderInput(session, "FB3", value = runif(1, 0, max_value))
    updateSliderInput(session, "FB4", value = runif(1, 0, max_value))
  })
  
  observeEvent(input$randomize_FV, {
    updateSliderInput(session, "FV1", value = runif(1, 0, max_value))
    updateSliderInput(session, "FV2", value = runif(1, 0, max_value))
    updateSliderInput(session, "FV3", value = runif(1, 0, max_value))
    updateSliderInput(session, "FV4", value = runif(1, 0, max_value))
  })
  
  observeEvent(input$randomize_all, {
    updateSliderInput(session, "FH1", value = runif(1, 0, max_value))
    updateSliderInput(session, "FH2", value = runif(1, 0, max_value))
    updateSliderInput(session, "FH3", value = runif(1, 0, max_value))
    updateSliderInput(session, "FH4", value = runif(1, 0, max_value))
    
    updateSliderInput(session, "FB1", value = runif(1, 0, max_value))
    updateSliderInput(session, "FB2", value = runif(1, 0, max_value))
    updateSliderInput(session, "FB3", value = runif(1, 0, max_value))
    updateSliderInput(session, "FB4", value = runif(1, 0, max_value))
    
    updateSliderInput(session, "FV1", value = runif(1, 0, max_value))
    updateSliderInput(session, "FV2", value = runif(1, 0, max_value))
    updateSliderInput(session, "FV3", value = runif(1, 0, max_value))
    updateSliderInput(session, "FV4", value = runif(1, 0, max_value))
  })
  
  output$distPlot <- renderPlot({
    # -------- PARAMETERS ----------------------------
    # Setup the parameters.
    tmax <- input$tmax
    
    k <- input$k
    b <- input$b
    FH <- c(input$FH1, input$FH2, input$FH3, input$FH4)  # Contribution to host pool
    FB <- c(input$FB1, input$FB2, input$FB3, input$FB4)  # Contribution to bacterial pool
    FV <- c(input$FV1, input$FV2, input$FV3, input$FV4)  # Contribution to viral pool
    
    # Absorption rates
    dBV <- input$dBV
    dHB <- input$dHB
    
    # Initial conditions from sliders
    initial_H <- input$initial_H
    initial_B <- input$initial_B
    initial_V <- input$initial_V
    
    # -------- DATA SETUP -------------------------
    # Preallocate timeseries data arrays based on simulation length, tmax.
    H <- rep(0, tmax) # Host
    B <- rep(0, tmax) # Bacteria
    V <- rep(0, tmax) # Virus
    B0 <- rep(0, tmax) # Bacteria with no virus
    B1 <- rep(0, tmax) # Bacteria with at least one virus
    H00 <- rep(0, tmax) # Hosts with 0 empty bact and 0 infected bact
    H01 <- rep(0, tmax) # Hosts with 0 empty bact and 1 infected bact
    H10 <- rep(0, tmax) # Hosts with at least 1 empty bact and 0 infected bact
    H11 <- rep(0, tmax) # Hosts with at least one of both
    
    # Set initial conditions
    H[1] <- initial_H
    B[1] <- initial_B
    V[1] <- initial_V
    
    # -------- FUNCTIONS ------------------------
    
    # Return the density parameter of bacteria and virus 
    muBV <- function(V, B) aBV(V / B) * V / B
    
    # Return the density parameter of hosts and bacteria
    muHB <- function(H, B) aHB(B / H) * B / H
    
    # Return the absorption rate of viruses by bacteria
    aBV <- function(r) dBV
    
    # Return the absorption rate of viruses by bacteria
    aHB <- function(r) dHB
    
    # Probability of getting exactly zero occurrences in a Poisson distribution
    P0 <- function(d) P(0, d)
    
    # Probability of getting at least one occurrence in a Poisson distribution
    P1 <- function(d) 1 - P0(d)
    
    # Compute probabilities in a Poisson distribution
    P <- function(m, d) {
      if (is.na(d)) {
        d <- 1e20  # If density goes to inf treat it as a huge number
      }
      
      if (d == 0) {
        return(ifelse(m == 0, 1, 0))
      }
      
      (exp(-d) * d^m) / factorial(m)
    }
    
    # Update the abundances of infected and uninfected bacteria at time t
    infect_bacteria <- function(t, B, V) {
      muBV_val <- muBV(V[t], B[t])
      B0 <- P0(muBV_val) * B[t] # Bacteria with no virus
      B1 <- P1(muBV_val) * B[t] # Bacteria with at least one virus
      list(B0 = B0, B1 = B1)
    }
    
    # Update the abundances of tripartite holobionts
    infected_bacteria_enter_hosts <- function(t, H, B0, B1) {
      muHB_B0 <- muHB(H[t], B0[t])
      muHB_B1 <- muHB(H[t], B1[t])
      H00 <- P0(muHB_B0) * P0(muHB_B1) * H[t]
      H01 <- P0(muHB_B0) * P1(muHB_B1) * H[t]
      H10 <- P1(muHB_B0) * P0(muHB_B1) * H[t]
      H11 <- P1(muHB_B0) * P1(muHB_B1) * H[t]
      list(H00 = H00, H01 = H01, H10 = H10, H11 = H11)
    }
    
    # Update the abundances of individuals for time t+1
    update_individuals <- function(t, H, B, V, H00, H01, H10, H11) {
      H[t+1] <- FH %*% c(H00, H01, H10, H11)
      B[t+1] <- FB %*% c(H00, H01, H10, H11)
      V[t+1] <- FV %*% c(H00, H01, H10, H11)
      list(H = H, B = B, V = V)
    }
    
    # =============== MAIN SIMULATION ==========================
    
    # Iterate the simulation
    for (t in 1:tmax) {
      bacteria <- infect_bacteria(t, B, V)
      B0[t] <- bacteria$B0
      B1[t] <- bacteria$B1
      
      hosts <- infected_bacteria_enter_hosts(t, H, B0, B1)
      H00[t] <- hosts$H00
      H01[t] <- hosts$H01
      H10[t] <- hosts$H10
      H11[t] <- hosts$H11
      
      if (t < tmax) {
        individuals <- update_individuals(t, H, B, V, H00[t], H01[t], H10[t], H11[t])
        H <- individuals$H
        B <- individuals$B
        V <- individuals$V
      }
    }
    
    # Make tibble of data and melt it
    data <- tibble(t = seq(1, tmax), H = H, B = B, V = V, B0 = B0, B1 = B1, 
                   H00 = H00, H01 = H01, H10 = H10, H11 = H11)
    
    melted_data <- data %>%
      pivot_longer(cols = -t, names_to = "variable", values_to = "value")
    
    # ---------------- PLOT ------------------------------
    theme_set(theme_gray())
    
    p <- ggplot(
      data = melted_data,
      mapping = aes(x = t, y = value, color = variable, linetype = variable)
    ) +
      geom_line(size = 1.5, alpha = 0.7) +  # Adjust alpha here
      scale_linetype_manual(values = c(
        H = "solid", B = "solid", V = "solid",
        H00 = "dotted", H01 = "dotted", H10 = "dotted", H11 = "dotted"
      )) +
      scale_y_log10() + 
      labs(x = "Time", y = "Abundance (Log)") +
      theme(legend.title = element_blank(), legend.position = "bottom")
    
    
    # Print the plot
    print(p)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
