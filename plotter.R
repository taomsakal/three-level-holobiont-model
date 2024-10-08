# ---------------------------#
# Virus Growth Rate Analysis
# ---------------------------#

# Load Required Libraries
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

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
  
  # Validate contributions
  contribution_inputs <- c(FH, FB, FV)
  if (any(is.na(contribution_inputs)) || any(contribution_inputs < 0)) {
    stop("Contribution values must be non-negative numbers.")
  }
  
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
calculate_avg_R0 <- function(sim_data, average_timesteps = 50) {
  # Ensure there are enough timesteps to calculate R0
  if (nrow(sim_data) < (average_timesteps + 1)) {
    warning(paste("Not enough timesteps to calculate R0. Required:", average_timesteps + 1, 
                  "Found:", nrow(sim_data)))
    return(NA)
  }
  
  # Extract the last (average_timesteps + 1) V values to calculate growth rates
  V_recent <- tail(sim_data$V, average_timesteps + 1)
  
  # Split into V_prev (V(t)) and V_next (V(t+1))
  V_prev <- V_recent[-length(V_recent)] # V(t)
  V_next <- V_recent[-1]               # V(t+1)
  
  # Calculate growth rates R(t) = V(t+1) / V(t)
  growth_rates <- V_next / V_prev
  
  # Handle cases where V(t) is zero to avoid division by zero
  growth_rates[V_prev == 0] <- NA
  growth_rates[is.infinite(growth_rates) | is.nan(growth_rates)] <- NA
  
  # Calculate average R0, excluding NA values
  R0_avg <- mean(growth_rates, na.rm = TRUE)
  
  return(R0_avg)
}

# -------- MAIN ANALYSIS ----------------------------

# Define Scenarios
scenarios <- define_scenarios()

# Define Burst Size (b) Values to Sweep Over
b_values <- seq(1, 800, by = 5)  # Burst sizes from 1 to 1000 in steps of 5

# Initialize a Data Frame to Store R0 Results
R0_results <- tibble(
  scenario = character(),
  description = character(),
  b = numeric(),
  R0 = numeric()
)

# Define Number of Timesteps Over Which to Average R0
average_timesteps <- 50  # Last 50 timesteps

# Total number of simulations for progress tracking
total_simulations <- length(scenarios) * length(b_values)
current_simulation <- 1

# Loop Through Each Scenario and Each b Value
for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  description <- scenario$description
  fixed_k <- scenario$fixed_params$k
  fixed_phi <- scenario$fixed_params$phi
  
  for (b in b_values) {
    # Run Simulation
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
      warning(paste("Simulation failed for scenario", scenario_name, "with b =", b, ":", e$message))
      return(NULL)
    })
    
    # Skip to next iteration if simulation failed
    if (is.null(sim_data)) {
      next
    }

    # Calculate average R0 over the last n timesteps
    R0_avg <- calculate_avg_R0(sim_data, average_timesteps)
    
    # Append to R0_results only if R0_avg is not NA
    if (!is.na(R0_avg)) {
      R0_results <- R0_results %>%
        add_row(
          scenario = scenario_name,
          description = description,
          b = b,
          R0 = R0_avg
        )
    }
    
    # Progress update
    cat(sprintf("Completed %d of %d simulations\n", current_simulation, total_simulations))
    current_simulation <- current_simulation + 1
  }
}

# Remove any rows with NA R0 (if any) - redundant now but kept for safety
R0_results <- R0_results %>% filter(!is.na(R0))

# Display the R0 Results
print(R0_results)

# -------- VISUALIZE RESULTS ----------------------------

# Beautify the plot
ggplot(R0_results, aes(x = b, y = R0, color = scenario, shape = scenario)) +
  geom_point(size = 2) +
  geom_line(aes(group = scenario), size = 1) +
  scale_color_manual(
    values = c(
      "bacteria_help_slightly_lytic" = "#1b9e77",
      "commensal_bacteria_slightly_lytic" = "#d95f02",
      "parasitic_bacteria_slightly_lytic" = "#7570b3",
      "bacterial_mutualism_more_nuance" = "#e7298a"
    ),
    labels = c(
      "Bacteria Help Hosts",
      "Commensal",
      "Parasitic",
      "Mutualistic"
    )
  ) +
  scale_shape_manual(
    values = c(
      "bacteria_help_slightly_lytic" = 16,
      "commensal_bacteria_slightly_lytic" = 17,
      "parasitic_bacteria_slightly_lytic" = 15,
      "bacterial_mutualism_more_nuance" = 18
    ),
    labels = c(
      "Bacteria Help Hosts",
      "Commensal",
      "Parasitic",
      "Mutualistic"
    )
  ) +
  labs(
    title = expression(paste("Average Virus Growth Rate ", R[0], " vs Burst Size (", b, ")")),
    x = "Burst Size (b)",
    y = expression(R[0]),
    color = "Scenario",
    shape = "Scenario"
  ) +
  #scale_x_log10() +  # Use logarithmic scale for burst size if necessary
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )
