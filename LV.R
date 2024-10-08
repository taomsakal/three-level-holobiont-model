# Load necessary packages
suppressPackageStartupMessages({
  if(!require(shiny)) {install.packages("shiny"); library(shiny)}
  if(!require(mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}  # For pmvnorm function
  if(!require(gtools)) {install.packages("gtools"); library(gtools)}    # Permutation with repeats                                               # Load feasibility analysis package
})

# Function that computes the average feasibility of a system 
# governed by the interaction matrix A
# Inputs: 
#   A = interaction matrix
#   num = number of iterations to repeat the calculation
# Output: 
#   out = the average feasibility of the system
Omega <- function(A, num = 30) {
  A <- as.matrix(A)
  S <- nrow(A)
  
  # Omega function (for a single iteration)
  omega <- function(S, Sigma) {
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    return(d[1])  # Community-level feasibility
  }
  
  # Rule out errors by checking if the matrix is solvable
  f <- function(m) class(try(solve(t(m) %*% m), silent = TRUE)) == "matrix"
  
  # Initialize output variable
  out <- 0
  
  # If matrix is not solvable, return 0
  if (all(f(A) == FALSE)) {
    return(0)
  } else {
    # Perform the calculation 'num' times and accumulate the results
    for (i in 1:num) {
      Sigma <- solve(t(A) %*% A)
      out <- out + omega(S, Sigma)
    }
    
    # Compute the average feasibility over the 'num' iterations
    out <- out / num
    
    return(out)
  }
}

# Function that computes the feasibility domain of a subset of species 
# in a large S(>2)-species community analytically
# Note that this is for just one region
# Inputs: A = S by S interaction matrix of the entire community, diag(A) = -1 
#         species = index vector of the subset of species
# Output: omega_comm_analytical = the normalized feasibility of the subset 
#         of species while others are transient
Omega_comm_analytical <- function(A, species){
  S <- nrow(A)
  if(is.null(species)){
    omega_comm_analytical <- Omega(-diag(S))
  }else{
    species <- sort(species)
    S <- nrow(A)
    B <- matrix(0, nrow = S, ncol = S)
    diag(B) <- 1
    
    other <- c(1:S)[-species]
    A_species <- cbind(A[,species], B[,other])
    omega_comm_analytical <- Omega(A_species)
  }
  omega_comm_analytical
}

# Function that computes the feasibility domain of a subset of species 
# in a large S(>2)-species community analytically
# Note that this is for a combination of regions
# Inputs: A = S by S interaction matrix of the entire community, diag(A) = -1 
#         species = index vector of the subset of species
# Output: omega_comm_analytical_all = the normalized feasibility of the subset 
#         of species across all possible comms
Omega_comm_analytical_all <- function(A, species){
  species <- sort(species)
  S <- nrow(A)
  S_sub <- length(species)
  S_other <- S - S_sub
  
  omega_comm_analytical_all <- Omega_comm_analytical(A, species)
  for(s in 1:S_other){
    other_comb <- combinations(S_other, s, c(1:S)[-species])
    n_other_comb <- nrow(other_comb)
    for(n in 1:n_other_comb){
      other <- other_comb[n,]
      comm <- sort(c(species, other))
      omega_comm_analytical_all <- omega_comm_analytical_all + Omega_comm_analytical(A, comm)
    }
  }
  omega_comm_analytical_all
}

# Shiny UI
ui <- fluidPage(
  titlePanel("Host-Bacteria-Phage Coexistence"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Set Interaction Parameters"),
      numericInput("A12", "Interaction (Bacteria on Host)", value = 0.3, step = 0.1),
      numericInput("A21", "Interaction (Host on Bacteria)", value = 0.3, step = 0.1),
      numericInput("A23", "Interaction (Phage on Bacteria)", value = 0.3, step = 0.1),
      numericInput("A32", "Interaction (Bacteria on Phage)", value = 0.3, step = 0.1),
      
      actionButton("calculate", "Calculate"),
      
      p("Note: aij < 0 indicates competitive effects, aij > 0 indicates facilitative effects. aii = -1."),
      p("Phage and Host do not interact, so their interaction values are set to zero.")
    ),
    
    mainPanel(
      h3("Results"),
      verbatimTextOutput("matrixOutput"),
      verbatimTextOutput("coexist_HBP"),
      verbatimTextOutput("hostContribution")
    )
  )
)

# Shiny Server
server <- function(input, output) {
  
  observeEvent(input$calculate, {
    # Define the interaction matrix A based on user inputs, with diagonals set to -1
    # Set the interaction between Phage (species 3) and Host (species 1) to 0
    A <- matrix(c(-1, input$A12, 0,
                  input$A21, -1, input$A23,
                  0, input$A32, -1), nrow = 3, byrow = TRUE)
    
    # Display the interaction matrix
    output$matrixOutput <- renderText({
      matrix_string <- paste(capture.output(print(A)), collapse = "\n")
      paste("Interaction matrix:\n", matrix_string)
    })
    
    # Calculate the probability of host-bacteria-phage coexistence
    coexist_HBP <- Omega(A)
    output$coexist_HBP <- renderPrint({
      paste("Probability of host-bacteria-phage coexistence: ", round(coexist_HBP, 3))
    })
    
    # Calculate the contribution of the host to bacteria-phage coexistence
    coexist_BP <- Omega(A[c(2, 3), c(2, 3)])
    coexist_BP_WithHost <- Omega_comm_analytical_all(A, species = c(2, 3))
    hostContribution <- coexist_BP_WithHost - coexist_BP
    output$hostContribution <- renderPrint({
      paste("Host contribution to bacteria-phage coexistence: ", round(hostContribution, 3))
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
