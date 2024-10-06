library(shiny)
library(ggplot2)
library(tidyverse)

# Maximum value for contribution sliders
max_value <- 3

# Define UI
ui <- fluidPage(
  titlePanel("Tripartite Model Simulation"),
  
  sidebarLayout(
    sidebarPanel(
      # Simulation parameters with help text for guidance
      sliderInput("tmax", "Simulation Time (tmax):",
                  min = 50, max = 500, value = 100, step = 10),
      helpText("Total number of time steps in the simulation."),
      
      # Uncomment if 'k' and 'b' are needed with their help texts
      # sliderInput("k", "Parameter k:",
      #             min = 1, max = 20, value = 10, step = 0.1),
      # helpText("Description of parameter k."),
      # sliderInput("b", "Parameter b:",
      #             min = 50, max = 200, value = 100, step = 1),
      # helpText("Description of parameter b."),
      
      # Absorption rates with explanations
      sliderInput("dBV", "Bacterial Absorption Rate of Viruses (dBV):",
                  min = 0.1, max = 1, value = 0.2, step = 0.1),
      helpText("Rate at which bacteria absorb viruses."),
      sliderInput("dHB", "Host Absorption Rate of Bacteria (dHB):",
                  min = 0.1, max = 1, value = 0.3, step = 0.1),
      helpText("Rate at which hosts absorb bacteria."),
      
      # Initial conditions with explanations
      sliderInput("initial_H", "Initial Hosts:",
                  min = 0, max = 100, value = 50, step = 5),
      helpText("Initial number of hosts in the simulation."),
      sliderInput("initial_B", "Initial Bacteria:",
                  min = 0, max = 100, value = 50, step = 5),
      helpText("Initial number of bacteria in the simulation."),
      sliderInput("initial_V", "Initial Viruses:",
                  min = 0, max = 100, value = 50, step = 5),
      helpText("Initial number of viruses in the simulation.")
    ),
    
    mainPanel(
      plotOutput("distPlot"),
      
      # Contribution sliders with explanations
      helpText("Contribution to the host next gen from holobionts H00, H01, H10, and H11."),
      fluidRow(
        column(3, sliderInput("FH1", "FH[1]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FH2", "FH[2]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FH3", "FH[3]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FH4", "FH[4]:", min = 0, max = max_value, value = 1, step = 0.1))
      ),
      
      helpText("Contribution to the bacteria next gen from holobionts H00, H01, H10, and H11."),
      fluidRow(
        column(3, sliderInput("FB1", "FB[1]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FB2", "FB[2]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FB3", "FB[3]:", min = 0, max = max_value, value = 1, step = 0.1)),
        column(3, sliderInput("FB4", "FB[4]:", min = 0, max = max_value, value = 1, step = 0.1))
      ),
      
      helpText("Contribution to the virus next gen from holobionts H00, H01, H10, and H11."),
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

# Define server logic
server <- function(input, output, session) {
  
  # Helper function to randomize slider inputs
  randomize_sliders <- function(prefix) {
    for (i in 1:4) {
      updateSliderInput(session, paste0(prefix, i), value = runif(1, 0, max_value))
    }
  }
  
  # Observers for randomize buttons
  observeEvent(input$randomize_FH, {
    randomize_sliders("FH")
  })
  
  observeEvent(input$randomize_FB, {
    randomize_sliders("FB")
  })
  
  observeEvent(input$randomize_FV, {
    randomize_sliders("FV")
  })
  
  observeEvent(input$randomize_all, {
    randomize_sliders("FH")
    randomize_sliders("FB")
    randomize_sliders("FV")
  })
  
  # Render the plot
  output$distPlot <- renderPlot({
    # -------- ERROR HANDLING ----------------------------
    # Validate required inputs
    req(input$tmax > 0)
    req(input$dBV > 0)
    req(input$dHB > 0)
    req(input$initial_H >= 0)
    req(input$initial_B >= 0)
    req(input$initial_V >= 0)
    
    # Check that contribution vectors are numeric and non-negative
    contribution_inputs <- c(
      input$FH1, input$FH2, input$FH3, input$FH4,
      input$FB1, input$FB2, input$FB3, input$FB4,
      input$FV1, input$FV2, input$FV3, input$FV4
    )
    if (any(is.na(contribution_inputs)) || any(contribution_inputs < 0)) {
      showNotification("Contribution values must be non-negative numbers.", type = "error")
      return(NULL)
    }
    
    # -------- PARAMETERS ----------------------------
    tmax <- input$tmax
    
    # Uncomment if 'k' and 'b' are needed
    # k <- input$k
    # b <- input$b
    
    # Contribution vectors
    FH <- c(input$FH1, input$FH2, input$FH3, input$FH4)
    FB <- c(input$FB1, input$FB2, input$FB3, input$FB4)
    FV <- c(input$FV1, input$FV2, input$FV3, input$FV4)
    
    # Absorption rates
    dBV <- input$dBV
    dHB <- input$dHB
    
    # Initial conditions
    H <- numeric(tmax)
    B <- numeric(tmax)
    V <- numeric(tmax)
    H[1] <- input$initial_H
    B[1] <- input$initial_B
    V[1] <- input$initial_V
    
    # Preallocate additional variables
    B0 <- numeric(tmax)  # Bacteria with no virus
    B1 <- numeric(tmax)  # Bacteria with at least one virus
    H00 <- numeric(tmax) # Hosts with 0 empty bact and 0 infected bact
    H01 <- numeric(tmax) # Hosts with 0 empty bact and 1 infected bact
    H10 <- numeric(tmax) # Hosts with at least 1 empty bact and 0 infected bact
    H11 <- numeric(tmax) # Hosts with at least one of both
    
    # -------- FUNCTIONS ------------------------
    
    # Return the absorption rate of viruses by bacteria (placeholder)
    aBV <- function(r) {
      dBV  # Currently returns a constant, can be extended later
    }
    
    # Return the absorption rate of bacteria by hosts (placeholder)
    aHB <- function(r) {
      dHB  # Currently returns a constant, can be extended later
    }
    
    # Density parameters
    muBV <- function(V_t, B_t) {
      ifelse(B_t <= 0, 0, aBV(B_t) * (V_t / B_t))
    }
    
    muHB <- function(H_t, B_t) {
      ifelse(H_t <= 0, 0, aHB(H_t) * (B_t / H_t))
    }
    
    # Probability functions
    P0 <- function(d) dpois(0, d)
    P1 <- function(d) 1 - P0(d)
    
    # =============== MAIN SIMULATION ==========================
    # Try-catch block for error handling
    tryCatch({
      for (t in 1:tmax) {
        # Update bacteria infection status
        muBV_val <- muBV(V[t], B[t])
        B0[t] <- P0(muBV_val) * B[t]
        B1[t] <- P1(muBV_val) * B[t]
        
        # Update host infection status
        muHB_B0 <- muHB(H[t], B0[t])
        muHB_B1 <- muHB(H[t], B1[t])
        
        H00[t] <- P0(muHB_B0) * P0(muHB_B1) * H[t]
        H01[t] <- P0(muHB_B0) * P1(muHB_B1) * H[t]
        H10[t] <- P1(muHB_B0) * P0(muHB_B1) * H[t]
        H11[t] <- P1(muHB_B0) * P1(muHB_B1) * H[t]
        
        # Ensure no negative values
        H00[t] <- max(0, H00[t])
        H01[t] <- max(0, H01[t])
        H10[t] <- max(0, H10[t])
        H11[t] <- max(0, H11[t])
        
        if (t < tmax) {
          # Update populations for next time step
          H[t + 1] <- sum(FH * c(H00[t], H01[t], H10[t], H11[t]))
          B[t + 1] <- sum(FB * c(H00[t], H01[t], H10[t], H11[t]))
          V[t + 1] <- sum(FV * c(H00[t], H01[t], H10[t], H11[t]))
          
          # Ensure no negative values
          H[t + 1] <- max(0, H[t + 1])
          B[t + 1] <- max(0, B[t + 1])
          V[t + 1] <- max(0, V[t + 1])
        }
      }
    }, error = function(e) {
      showNotification("An error occurred during simulation. Please check your inputs.", type = "error")
      return(NULL)
    })
    
    # -------- DATA PREPARATION -------------------------
    # Create a data frame with simulation results
    data <- tibble(
      t = 1:tmax,
      H = H, B = B, V = V,
      B0 = B0, B1 = B1,
      H00 = H00, H01 = H01, H10 = H10, H11 = H11
    )
    
    # Remove non-finite or NA values
    data <- data %>% filter_all(all_vars(is.finite(.)))
    
    # Transform data for plotting
    melted_data <- data %>%
      pivot_longer(cols = -t, names_to = "variable", values_to = "value") %>%
      filter(value > 0)  # Remove zero or negative values for log scale
    
    # Check if there is data to plot
    validate(
      need(nrow(melted_data) > 0, "No data to plot. Please adjust your inputs.")
    )
    
    # -------- PLOT ------------------------------
    ggplot(
      data = melted_data,
      mapping = aes(x = t, y = value, color = variable, linetype = variable)
    ) +
      geom_line(size = 1.2, alpha = 0.8) +
      scale_linetype_manual(values = c(
        H = "solid", B = "solid", V = "solid",
        B0 = "dashed", B1 = "dashed",
        H00 = "dotted", H01 = "dotted", H10 = "dotted", H11 = "dotted"
      )) +
      scale_y_log10() +
      labs(x = "Time", y = "Abundance (Log Scale)", color = "Variables", linetype = "Variables") +
      theme_minimal() +
      theme(
        legend.position = "bottom"
      )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
