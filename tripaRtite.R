# Load Required Libraries
library(shiny)
library(ggplot2)
library(tidyverse)
library(shinyWidgets)  # For knob inputs

# Define Maximum Value and Knob Size
max_value <- 3
knob_size <- "80px"  # Adjust this value to make knobs smaller or larger

# Define Server Logic First
server <- function(input, output, session) {
  
  # Helper Function to Randomize Knob Inputs
  randomize_sliders <- function(prefix) {
    for (i in 1:4) {
      new_val <- runif(1, 0, max_value)
      updateKnobInput(session, paste0(prefix, i), value = new_val)
    }
  }
  
  # Observers for Randomize Buttons
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
  
  # Render the Plot
  output$distPlot <- renderPlot({
    # -------- ERROR HANDLING ----------------------------
    # Validate Required Inputs
    req(input$tmax > 0)
    req(input$dBV > 0)
    req(input$dHB > 0)
    req(input$initial_H >= 0)
    req(input$initial_B >= 0)
    req(input$initial_V >= 0)
    
    # Check that Contribution Vectors are Numeric and Non-negative
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
    
    # Contribution Vectors
    FH <- c(input$FH1, input$FH2, input$FH3, input$FH4)
    FB <- c(input$FB1, input$FB2, input$FB3, input$FB4)
    FV <- c(input$FV1, input$FV2, input$FV3, input$FV4)
    
    # Absorption Rates
    dBV <- input$dBV
    dHB <- input$dHB
    
    # Initial Conditions
    H <- numeric(tmax)
    B <- numeric(tmax)
    V <- numeric(tmax)
    H[1] <- input$initial_H
    B[1] <- input$initial_B
    V[1] <- input$initial_V
    
    # Preallocate Additional Variables
    B0 <- numeric(tmax)  # Bacteria with no virus
    B1 <- numeric(tmax)  # Bacteria with at least one virus
    H00 <- numeric(tmax) # Hosts with 0 empty bact and 0 infected bact
    H01 <- numeric(tmax) # Hosts with 0 empty bact and 1 infected bact
    H10 <- numeric(tmax) # Hosts with at least 1 empty bact and 0 infected bact
    H11 <- numeric(tmax) # Hosts with at least one of both
    
    # -------- FUNCTIONS ------------------------
    
    # Absorption Rate of Viruses by Bacteria (Placeholder)
    aBV <- function(r) {
      dBV  # Currently returns a constant, can be extended later
    }
    
    # Absorption Rate of Bacteria by Hosts (Placeholder)
    aHB <- function(r) {
      dHB  # Currently returns a constant, can be extended later
    }
    
    # Density Parameters
    muBV <- function(V_t, B_t) {
      ifelse(B_t <= 0, 0, aBV(B_t) * (V_t / B_t))
    }
    
    muHB <- function(H_t, B_t) {
      ifelse(H_t <= 0, 0, aHB(H_t) * (B_t / H_t))
    }
    
    # Probability Functions
    P0 <- function(d) dpois(0, d)
    P1 <- function(d) 1 - P0(d)
    
    # =============== MAIN SIMULATION ==========================
    # Try-Catch Block for Error Handling
    tryCatch({
      for (t in 1:tmax) {
        # Update Bacteria Infection Status
        muBV_val <- muBV(V[t], B[t])
        B0[t] <- P0(muBV_val) * B[t]
        B1[t] <- P1(muBV_val) * B[t]
        
        # Update Host Infection Status
        muHB_B0 <- muHB(H[t], B0[t])
        muHB_B1 <- muHB(H[t], B1[t])
        
        H00[t] <- P0(muHB_B0) * P0(muHB_B1) * H[t]
        H01[t] <- P0(muHB_B0) * P1(muHB_B1) * H[t]
        H10[t] <- P1(muHB_B0) * P0(muHB_B1) * H[t]
        H11[t] <- P1(muHB_B0) * P1(muHB_B1) * H[t]
        
        # Ensure No Negative Values
        H00[t] <- max(0, H00[t])
        H01[t] <- max(0, H01[t])
        H10[t] <- max(0, H10[t])
        H11[t] <- max(0, H11[t])
        
        if (t < tmax) {
          # Update Populations for Next Time Step
          H[t + 1] <- sum(FH * c(H00[t], H01[t], H10[t], H11[t]))
          B[t + 1] <- sum(FB * c(H00[t], H01[t], H10[t], H11[t]))
          V[t + 1] <- sum(FV * c(H00[t], H01[t], H10[t], H11[t]))
          
          # Ensure No Negative Values
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
    # Create a Data Frame with Simulation Results
    data <- tibble(
      t = 1:tmax,
      H = H, B = B, V = V,
      B0 = B0, B1 = B1,
      H00 = H00, H01 = H01, H10 = H10, H11 = H11
    )
    
    # Remove Non-finite or NA Values
    rdata <- data %>% filter_all(all_vars(is.finite(.)))
    
    # Transform Data for Plotting
    melted_data <- data %>%
      pivot_longer(cols = -t, names_to = "variable", values_to = "value") %>%
      filter(value > 0)  # Remove zero or negative values for log scale
    
    # Check if There is Data to Plot
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

# Define UI Below the Server Logic
ui <- fluidPage(
  # Custom CSS for Enhanced Aesthetics
  tags$head(
    tags$style(HTML("
      /* Style for the contribution matrix table */
      table {
        width: 100%;
        border-collapse: collapse;
        margin-bottom: 20px;
      }
      th, td {
        padding: 10px;
        text-align: center;
        border: 1px solid #ddd;
      }
      th {
        background-color: #f2f2f2;
        font-weight: bold;
      }
      /* Optional: Hover effect for table cells */
      td:hover {
        background-color: #f9f9f9;
      }
    "))
  ),
  
  titlePanel("Tripartite Model Simulation"),
  
  sidebarLayout(
    sidebarPanel(
      # Simulation Parameters with Help Text
      sliderInput("tmax", "Simulation Time (tmax):",
                  min = 50, max = 500, value = 100, step = 10),
      
      # Absorption Rates with Explanations
      sliderInput("dBV", "Bacterial Absorption Rate of Viruses (dBV):",
                  min = 0.0, max = 1, value = 0.2, step = 0.001),
      
      sliderInput("dHB", "Host Absorption Rate of Bacteria (dHB):",
                  min = 0.0, max = 1, value = 0.3, step = 0.001),
      
      # Initial Conditions with Explanations
      sliderInput("initial_H", "Initial Hosts:",
                  min = 0, max = 10000, value = 50, step = 1),
      sliderInput("initial_B", "Initial Bacteria:",
                  min = 0, max = 10000, value = 50, step = 1),
      sliderInput("initial_V", "Initial Viruses:",
                  min = 0, max = 10000, value = 50, step = 1),

    ),
    
    mainPanel(
      plotOutput("distPlot"),
      
      # Explanatory Text Above the Contribution Matrix
      helpText("Contribution Matrix: how many new individuals (Host, Bacteria, or Virus) each holobiont (H00, H01, H10, or H11) creates."),
      
      # Contribution Matrix Grid
      tags$table(
        # Header Row with Column Names
        tags$tr(
          tags$th(NULL),  # Empty top-left corner
          lapply(c("H00", "H01", "H10", "H11"), function(col_label) {
            tags$th(col_label)
          })
        ),
        # Hosts Row
        tags$tr(
          tags$td("Hosts"),
          lapply(c("FH1", "FH2", "FH3", "FH4"), function(input_id) {
            tags$td(
              knobInput(
                inputId = input_id,
                label = NULL,  # No label inside the cell
                value = 1,
                min = 0,
                max = max_value,
                step = 0.1,
                displayPrevious = FALSE,
                fgColor = "#66C2A5",  # Host-related color
                inputColor = "#1F78B4",
                width = knob_size  # Use the knob_size variable
              )
            )
          })
        ),
        # Bacteria Row
        tags$tr(
          tags$td("Bacteria"),
          lapply(c("FB1", "FB2", "FB3", "FB4"), function(input_id) {
            tags$td(
              knobInput(
                inputId = input_id,
                label = NULL,
                value = 1,
                min = 0,
                max = max_value,
                step = 0.1,
                displayPrevious = FALSE,
                fgColor = "#FC8D62",  # Bacteria-related color
                inputColor = "#E78AC3",
                width = knob_size
              )
            )
          })
        ),
        # Viruses Row
        tags$tr(
          tags$td("Viruses"),
          lapply(c("FV1", "FV2", "FV3", "FV4"), function(input_id) {
            tags$td(
              knobInput(
                inputId = input_id,
                label = NULL,
                value = 1,
                min = 0,
                max = max_value,
                step = 0.1,
                displayPrevious = FALSE,
                fgColor = "#8DA0CB",  # Virus-related color
                inputColor = "#E5C494",
                width = knob_size
              )
            )
          })
        )
      ),
      
      # Randomize Buttons Below the Grid
      fluidRow(
        column(3, actionButton("randomize_FH", "Randomize FH")),
        column(3, actionButton("randomize_FB", "Randomize FB")),
        column(3, actionButton("randomize_FV", "Randomize FV")),
        column(3, actionButton("randomize_all", "Randomize All"))
      )
    )
  )
)

# Run the Application
shinyApp(ui = ui, server = server)
