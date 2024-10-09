# ---------------------------#
# Virus Growth Rate Analysis
# ---------------------------#

# Load Required Libraries
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

# Libraries for parallel processing
library(furrr)
library(future)
library(progressr)

# ---------------------------#
# Define Biological Scenarios
# ---------------------------#
# Each scenario represents a different biological interaction between hosts, bacteria, and viruses.
# Scenarios include how bacteria affect hosts and how viruses behave (e.g., lytic).
#
# A line like FH <- c(1,2,3,4) means that each holobiont H00 creates one new
# host, H01 creates two, H10 creates three, and H11 creates four. FB and FV
# calculate the next generation of bacteria and viruses, respectively.
#
# Recall that
# H00 = host with no bacteria at all
# H01 = host with at least one infected bacteria, no uninfected bacteria
# H10 = host with no infected bacteria, at least one uninfected
# H11 = host with at least one of both kinds of bacteria
#
# Note that we have defined the senarios as a function to allow for easy 
# modification and reuse.

define_scenarios <- function() {
  scenarios <- list(
    bacteria_help_slightly_lytic = list(
      description = "Bacteria help host, viruses are slightly lytic",
      compute_contributions = function(k, phi, b) {
        FH <- c(1, 2, 2, 2)
        FB <- c(0, k / phi, k, k / phi)
        FV <- c(0, phi * b, 0, phi * b)
        list(FH = FH, FB = FB, FV = FV)
      },
      # Fixed parameters for this scenario
      fixed_params = list(k = 100, phi = 2)
    ),
    
    commensal_bacteria_slightly_lytic = list(
      description = "Bacteria is commensal, viruses are slightly lytic",
      compute_contributions = function(k, phi, b) {
        FH <- c(2, 2, 2, 2)
        FB <- c(0, k / phi, k, k / phi)
        FV <- c(0, phi * b, 0, phi * b)
        list(FH = FH, FB = FB, FV = FV)
      },
      # Fixed parameters for this scenario
      fixed_params = list(k = 100, phi = 3)
    ),
    
    parasitic_bacteria_slightly_lytic = list(
      description = "Bacteria is parasitic, viruses are slightly lytic",
      compute_contributions = function(k, phi, b) {
        FH <- c(2, 1.5, 1.5, 1.5)
        FB <- c(0, k / phi, k, k / phi)
        FV <- c(0, phi * b, 0, phi * b)
        list(FH = FH, FB = FB, FV = FV)
      },
      # Fixed parameters for this scenario
      fixed_params = list(k = 100, phi = 3)
    ),
    
    bacterial_mutualism_more_nuance = list(
      description = "Bacterial mutualism, more nuance",
      compute_contributions = function(k, phi, b) {
        FH <- c(0.1, 2.5, 4, 3)
        FB <- c(0, k / 4, k, k / 2)
        FV <- c(0, b, 0, 2 * b)
        list(FH = FH, FB = FB, FV = FV)
      },
      # Fixed parameters for this scenario
      fixed_params = list(k = 10, phi = 2)
    )
  )
  
  return(scenarios)
}

# -------- SIMULATION FUNCTION ----------------------------
simulate_dynamics <- function(scenario, k, phi, b, 
                              tmax = 200, 
                              dBV = 0.1, dHB = 0.05, 
                              initial_H = 50, initial_B = 100, initial_V = 200) {
  
  # Compute contribution vectors based on the scenario
  contributions <- scenario$compute_contributions(k, phi, b)
  FH <- contributions$FH
  FB <- contributions$FB
  FV <- contributions$FV
  
  # Density Parameters
  muBV <- function(V_t, B_t) {
    ifelse(B_t <= 0, 0, dBV * (V_t / B_t))
  }
  
  muHB <- function(H_t, B_t) {
    ifelse(H_t <= 0, 0, dHB * (B_t / H_t))
  }
  
  # Probability Functions
  P0 <- function(d) dpois(0, d)
  P1 <- function(d) 1 - P0(d)
  
  # Initialize Variables
  H <- numeric(tmax)
  B <- numeric(tmax)
  V <- numeric(tmax)
  H[1] <- initial_H
  B[1] <- initial_B
  V[1] <- initial_V
  
  # Preallocate data lists for timeseries
  B0 <- numeric(tmax)  # Bacteria with no virus
  B1 <- numeric(tmax)  # Bacteria with at least one virus
  H00 <- numeric(tmax) # Hosts with 0 empty bact and 0 infected bact
  H01 <- numeric(tmax) # Hosts with 0 empty bact and 1 infected bact
  H10 <- numeric(tmax) # Hosts with at least 1 empty bact and 0 infected bact
  H11 <- numeric(tmax) # Hosts with at least one of both
  
  # Main Simulation Loop
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
  
  # Create a Data Frame with Simulation Results
  data <- tibble(
    t = 1:tmax,
    H = H, 
    B = B, 
    V = V,
    B0 = B0, 
    B1 = B1,
    H00 = H00, 
    H01 = H01, 
    H10 = H10, 
    H11 = H11
  )
  
  return(data)
}


# -------- R0 CALCULATION FUNCTION ----------------------------
# Calculate the average growth rate of the given population over the
# last timesteps.
calculate_avg_R0 <- function(sim_data, population = "V", n = 50) {
  # Validate population input
  if (!population %in% c("H", "B", "V")) {
    stop("Population must be one of 'H', 'B', or 'V'.")
  }
  
  # Ensure there are enough timesteps to calculate R0
  if (nrow(sim_data) < (n + 1)) {
    warning(paste("Not enough timesteps to calculate R0 for", population, 
                  ". Required:", n + 1, 
                  "Found:", nrow(sim_data)))
    return(NA)
  }
  
  # Select the appropriate column based on population
  population_col <- sim_data[[population]]
  
  # Extract the last (n + 1) values to calculate growth rates
  pop_recent <- tail(population_col, n + 1)
  
  # Split into pop_prev (Pop(t)) and pop_next (Pop(t+1))
  pop_prev <- pop_recent[-length(pop_recent)] # Pop(t)
  pop_next <- pop_recent[-1]                   # Pop(t+1)
  
  # Calculate growth rates R0(t) = Pop(t+1) / Pop(t)
  growth_rates <- pop_next / pop_prev
  
  # Handle cases where Pop(t) is zero to avoid division by zero
  growth_rates[pop_prev == 0] <- NA
  growth_rates[is.infinite(growth_rates) | is.nan(growth_rates)] <- NA
  
  # Calculate average R0, excluding NA values
  R0_avg <- mean(growth_rates, na.rm = TRUE)
  
  return(R0_avg)
}



# -------- MAIN ANALYSIS WITH FUTURE AND FURRR PARALLEL PROCESSING ----------------------------
# Sweep through the various bust sizes. We do this in a parallel manner.
# (Note the parallelzation part is mostly written by GPT o1-mini. Seems to work
# the same as before in all my tests, just must faster.)

# Define Scenarios
scenarios <- define_scenarios()

# Define Burst Size (b) Values to Sweep Over
b_values <- seq(1, 800, by = 5)  # Burst sizes from 1 to 800 in steps of 5

# Create a Grid of All Scenario and b Combinations
simulation_grid <- expand.grid(
  scenario_name = names(scenarios),
  b = b_values,
  stringsAsFactors = FALSE
)

# Set Up Future Plan
# Use multisession for parallel processing on multiple cores
# You can also use multicore on Unix-like systems, but multisession is more portable
num_cores <- future::availableCores() - 1  # Reserve one core for the system
plan(multisession, workers = num_cores)

# Initialize Progress Handler
handlers(global = TRUE)
handlers("txtprogressbar")  # Use text-based progress bar

# Define a progress handler using progressr and furrr
progressr::with_progress({
  
  p <- progressr::progressor(along = 1:nrow(simulation_grid))
  
  # Define a function to run a single simulation
  run_simulation <- function(scenario_name, b) {
    # Retrieve the scenario details
    scenario <- scenarios[[scenario_name]]
    description <- scenario$description
    fixed_k <- scenario$fixed_params$k
    fixed_phi <- scenario$fixed_params$phi
    
    # Run Simulation with Error Handling
    sim_data <- tryCatch({
      simulate_dynamics(
        scenario = scenario,
        k = fixed_k,
        phi = fixed_phi,
        b = b,
        tmax = 200,         # Maximum time steps
        dBV = 0.2,          # Absorption rate dBV
        dHB = 0.1,          # Absorption rate dHB
        initial_H = 50,     # Initial Hosts
        initial_B = 100,    # Initial Bacteria
        initial_V = 200     # Initial Viruses
      )
    }, error = function(e) {
      # Log the error message
      warning(paste("Simulation failed for scenario", scenario_name, "with b =", b, ":", e$message))
      return(NULL)
    })
    
    # If simulation failed, return NULL
    if (is.null(sim_data)) {
      return(NULL)
    }
    
    # Define Number of Timesteps Over Which to Average R0
    average_timesteps <- 50  # Last 50 timesteps
    
    # Calculate average R0 for Hosts, Bacteria, and Viruses
    R0_H_avg <- calculate_avg_R0(sim_data, population = "H", average_timesteps)
    R0_B_avg <- calculate_avg_R0(sim_data, population = "B", average_timesteps)
    R0_V_avg <- calculate_avg_R0(sim_data, population = "V", average_timesteps)
    
    # Update progress
    p()
    
    # Return a tibble row if all R0_avg are not NA
    if (!is.na(R0_H_avg) && !is.na(R0_B_avg) && !is.na(R0_V_avg)) {
      return(tibble(
        scenario = scenario_name,
        description = description,
        R0_H = R0_H_avg,
        R0_B = R0_B_avg,
        R0_V = R0_V_avg
      ))
    } else {
      return(NULL)
    }
  }
  
  # Execute simulations in parallel with furrr::future_pmap
  R0_results <- simulation_grid %>%
    mutate(result = furrr::future_pmap(
      list(scenario_name, b),
      run_simulation
    )) %>%
    # Remove rows where result is NULL (failed simulations)
    filter(!sapply(result, is.null)) %>%
    # Unnest the result column to get a flat data frame
    unnest(result)
})

# Reset future plan to sequential
plan(sequential)

# Display the R0 Results
print(R0_results)

# -------- VISUALIZE RESULTS ----------------------------

# Reshape the R0_results data to long format
R0_long <- R0_results %>%
  pivot_longer(
    cols = starts_with("R0_"),
    names_to = "Population",
    values_to = "R0"
  ) %>%
  mutate(
    Population = recode(Population,
                        "R0_H" = "Host",
                        "R0_B" = "Bacteria",
                        "R0_V" = "Virus")
  )

# Define color and shape mappings for scenarios
scenario_colors <- c(
  "bacteria_help_slightly_lytic" = "#1b9e77",
  "commensal_bacteria_slightly_lytic" = "#d95f02",
  "parasitic_bacteria_slightly_lytic" = "#7570b3",
  "bacterial_mutualism_more_nuance" = "#e7298a"
)

scenario_shapes <- c(
  "bacteria_help_slightly_lytic" = 16,
  "commensal_bacteria_slightly_lytic" = 17,
  "parasitic_bacteria_slightly_lytic" = 15,
  "bacterial_mutualism_more_nuance" = 18
)

# Create the combined plot with facets
combined_plot <- ggplot(R0_long, aes(x = b, y = R0, color = scenario, shape = scenario)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_line(aes(group = scenario), size = 1, alpha = 0.7) +
  scale_color_manual(
    values = scenario_colors,
    labels = c(
      "bacteria_help_slightly_lytic" = "Bacteria Help Hosts",
      "commensal_bacteria_slightly_lytic" = "Commensal",
      "parasitic_bacteria_slightly_lytic" = "Parasitic",
      "bacterial_mutualism_more_nuance" = "Mutualistic"
    )
  ) +
  scale_shape_manual(
    values = scenario_shapes,
    labels = c(
      "bacteria_help_slightly_lytic" = "Bacteria Help Hosts",
      "commensal_bacteria_slightly_lytic" = "Commensal",
      "parasitic_bacteria_slightly_lytic" = "Parasitic",
      "bacterial_mutualism_more_nuance" = "Mutualistic"
    )
  ) +
  labs(
    title = "Average Growth Rates (R0) vs Burst Size (b) for Hosts, Bacteria, and Viruses",
    x = "Burst Size (b)",
    y = expression(R[0]),
    color = "Scenario",
    shape = "Scenario"
  ) +
  facet_wrap(~ Population, scales = "free_y") +  # Create separate panels for each population
  expand_limits(y = 0) +  # Ensure y-axis starts at zero
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 14)  # Increase facet label size
  )

# Display the combined plot
print(combined_plot)

# Optionally, save the plot
# ggsave("Combined_R0_vs_b.png", plot = combined_plot, width = 12, height = 8, dpi = 300)
