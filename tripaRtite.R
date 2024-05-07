library(ggplot2)
library(tidyverse)


# -------- DATA SETUP -------------------------
# Here we store the timeseries data of individuals and holobiont abundances.
# These arrays are preallocated based off simulation length, tmax.
tmax = 100
H = rep(0, tmax) # Host
B = rep(0, tmax) # Bacteria
V = rep(0, tmax) # Virus
B0 = rep(0, tmax) # Bacteria with no virus
B1 = rep(0, tmax) # Bacteria with at least one virus
H00 <- rep(0, tmax) # Hosts with 0 empty bact and 0 infected bact
H01 <- rep(0, tmax) # Hosts with 0 empty bact and 1 infected bact
H10 <- rep(0, tmax) # Hosts with at least 1 empty bact and 0 infected bact
H11 <- rep(0, tmax) # Hosts with at least one of both


# -------- PARAMETERS ----------------------------
# Here we setup the parameters.

# Bacterial mutualism case
# todo: figure out what k, phi, and b were and comment that in
k <- 10 
phi <- 2
b <- 100
FH <- c(.9, 1.1, 1.4, 1.2)  # Contribution to host pool
FB <- c(0, k/2, k, k/2)  # Contribution to bacterial pool
FV <- c(0, b, 0, b*2)  # Contribution to viral pool




# ------------ SIMULATION FUNCTIONS ------------------------
# All the functions that run the simulation



#' Return the density parameter of bacteria and virus 
muBV <- function(t) aBV(V[t]/B[t]) * V[t]/B[t]


#' Return the density parameter of hosts and bacteria
muHB <- function(H, B) aHB(B/H) * B/H  # Todo: raise error if this is NaN

#' Return the absorption rate(?) of viruses by bacteria
#' 
#' @param r virus/bacteria ratio
aBV <- function(r) dBV
dBV <- 0.2 # right now the rates are constant.

#' Return the absorption rate(?) of viruses by bacteria
#' @param virus/bacteria ratio
aHB <- function(r) dHB
dHB <- 0.3


#' Probability of getting exactly zero occurrences.
#' 
#' @param d Density parameter.
#' @return Probability of getting exactly zero occurrences.
P0 <- function(d) P(0, d)

#' Probability of getting at least one occurrence.
#' 
#' @param d Density parameter.
#' @return Probability of getting at least one occurrence.
P1 <- function(d) {
  return(1 - P0(d))
}

#' Functions to compute probabilities in a Poisson distribution.
#' 
#' @param m Number of occurrences.
#' @param d Density parameter.
#' @return Probability of getting exactly \code{m} occurrences.
P <- function(m, d) {
  if (is.na(d)) {
    d = 1e20  # If density goes to inf treat it as a huge number
  }
  
  if (d == 0) {
    if (m == 0) return(1) else return(0)
  }
  return((exp(-d) * d^m) / factorial(m))
}


#' Update the abundances of infected and uninfected bacteria at time t
#' 
#' This function uses the previous timestep to calculate the abundence of 
#' uninfected bacteria (B0) and infected bacteria (B1). It then updates the
#' data arrays with those values.
infect_bacteria = function(t){
  B0[t] <<- P0(muBV(t)) * B[t] # Bacteria with no virus
  B1[t] <<- P1(muBV(t)) * B[t] # Bacteria with at least one virus
}

#' Update the abundances of tripartite holobionts
infected_bacteria_enter_hosts = function(t){
  H00[t] <<- P0(muHB(H[t], B0[t])) * P0(muHB(H[t], B1[t])) * H[t]
  H01[t] <<- P0(muHB(H[t], B0[t])) * P1(muHB(H[t], B1[t])) * H[t]
  H10[t] <<- P1(muHB(H[t], B0[t])) * P0(muHB(H[t], B1[t])) * H[t]
  H11[t] <<- P1(muHB(H[t], B0[t])) * P1(muHB(H[t], B1[t])) * H[t]
}

#' Update the abundances of individuals
#' 
#' This splits up the holobionts at time t to update the individual
#' abundences for time t+1. This is the start of the next generation.
update_individuals = function(t){
  # Linear contribution functions
  H[t+1] <<- FH %*% c(H00[t], H01[t], H10[t], H11[t])
  B[t+1] <<- FB %*% c(H00[t], H01[t], H10[t], H11[t])
  V[t+1] <<- FV %*% c(H00[t], H01[t], H10[t], H11[t])
}


# =============== MAIN SIMULATION ==========================

# Random Initial Conditions
H[1] = runif(1)*100
B[1] = runif(1)*100
V[1] = runif(1)*100

# Iterate the simulation
for (t in 1:tmax) {
  infect_bacteria(t) 
  infected_bacteria_enter_hosts(t) #
  if (t < tmax){
  update_individuals(t) 
  }
}

# Make tibble of data and melt it
data <- tibble(t = seq(1,tmax), H = H, B = B, V = V, B0 = B0, B1 = B1, 
               H00 = H00, H01 = H01, H10 = H10, H11 = H11)

melted_data <- data %>%
  pivot_longer(cols = -t, names_to = "variable", values_to = "value")




# ---------------- PLOT ------------------------------
# Todo: Clean up the code, choose better colors, 
# clean up legend, maybe split into two graphs



theme_set(theme_gray())

p <- ggplot(
  data = melted_data,
  mapping = aes(x = t, y = value, color = variable, linetype = variable)
) +
  geom_line(size = 1.5, alpha = 0.7) +  # Adjust alpha here
  scale_linetype_manual(values = c(H = "solid", B = "solid", V = "solid", 
                                   H00 = "dotted", H01 = "dotted", H10 = "dotted", H11 = "dotted")) +
  labs(x = "Time", y = "Abundance") 

# Print the plot
print(p)

# Print the log plot
print(p + scale_y_log10())


