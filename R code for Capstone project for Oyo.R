install.packages("deSolve")
library(deSolve)


# Define the SIR-SI model differential equations
sir_si_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Human equations
    dS_H <- (Lambda_H + theta_H) - (beta_H * (1 - psi) * S_H * I_V)/N_H - (mu_H * S_H)
    dI_H <- (beta_H * (1 - psi) * S_H * I_V)/N_H - (gamma_H * (1 + sigma) + mu_H + delta_H) * I_H
    dR_H <- gamma_H * (1 + sigma) * I_H - (mu_H + theta_H) * R_H
    
    # Mosquito equations
    dS_V <- Lambda_V - (beta_V * (1 - psi) * S_V * I_V)/N_H - (mu_V * S_V)
    dI_V <- (beta_V * (1 - psi) * S_V * I_V)/N_H  - (mu_V * I_V)
    
    # Return the rate of change
    list(c(dS_H, dI_H, dR_H, dS_V, dI_V))
  })
}

# Initial state values
N_H <-  7358018  # Total human population
S_H0 <- 3985103  # Initial susceptible humans
I_H0 <- 1765924  # Initial infected humans
R_H0 <- 1606991  # Initial recovered humans
S_V0 <- 1000000  # Initial susceptible mosquitoes (assume 1 million)
I_V0 <- 10000  # Initial infected mosquitoes (assume 10,000)

initial_state <- c(S_H = S_H0, I_H = I_H0, R_H = R_H0, S_V = S_V0, I_V = I_V0)

# Define parameters
parameters <- c(
  Lambda_H = 26182, # human recruitment rate
  Lambda_V = 0, # Assume constant mosquito population for simplicity
  beta_H = 0.0001,  # Transmission rate from mosquitoes to humans
  beta_V = 0.0002,  # Transmission rate from humans to mosquitoes
  gamma_H = 0.0714,  # Recovery rate in humans (1/14 days)
  mu_H = 0.0009,  # Human mortality rate
  mu_V = 0.1,  # Mosquito mortality rate
  psi = 0.6,  # ITN effectiveness
  sigma = 0.5,  # Antimalarial drug effectiveness
  theta_H = 0.001, # waning of immunity 
  delta_H = 0.09) #death induced


# Time frame for simulation (e.g., 1 year with daily steps)
time <- seq(0, 72, by = 1)

# Solve the differential equations
output <- ode(y = initial_state, times = time, func = sir_si_model, parms = parameters)

print(output)
# Convert output to a data frame
output_df <- as.data.frame(output)


# Plot the results
plot(output_df$time, output_df$S_H, type = "l", col = "blue", ylim = c(0, max(output_df$S_H, output_df$I_H, output_df$R_H)), 
     xlab = "Time (days)", ylab = "Population", main = "SIR-SI Model: Oyo with 60% intervention")
lines(output_df$time, output_df$I_H, col = "red")
lines(output_df$time, output_df$R_H, col = "green")
lines(output_df$time, output_df$S_V, col = "orange")
legend("right", legend = c("Susceptible (Humans)", "Infected (Humans)", "Recovered (Humans)"), col = c("blue", "red", "green"), lty = 1)

# Plot Mosquito population
plot(output_df$time, output_df$S_V, type = "l", col = "orange", ylim = c(0, max(output_df$S_V, output_df$I_V)), 
     xlab = "Time (days)", ylab = "Mosquito Population", main = "SIR-SI Model: Kano Mosquito Population")
lines(output_df$time, output_df$I_V, col = "purple")
legend("right", legend = c("Susceptible (Mosquitoes)", "Infected (Mosquitoes)"), col = c("orange", "purple"), lty = 1)


############################################################################################
# with zero Interventions
install.packages("deSolve")
library(deSolve)


# Define the SIR-SI model differential equations
sir_si_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Human equations
    dS_H <- (Lambda_H + theta_H) - (beta_H * (1 - psi) * S_H * I_V)/N_H - (mu_H * S_H)
    dI_H <- (beta_H * (1 - psi) * S_H * I_V)/N_H - (gamma_H * (1 + sigma) + mu_H + delta_H) * I_H
    dR_H <- gamma_H * (1 + sigma) * I_H - (mu_H + theta_H) * R_H
    
    # Mosquito equations
    dS_V <- Lambda_V - (beta_V * (1 - psi) * S_V * I_V)/N_H - (mu_V * S_V)
    dI_V <- (beta_V * (1 - psi) * S_V * I_V)/N_H  - (mu_V * I_V)
    
    # Return the rate of change
    list(c(dS_H, dI_H, dR_H, dS_V, dI_V))
  })
}

# Initial state values
N_H <-  7358018  # Total human population
S_H0 <- 3985103  # Initial susceptible humans
I_H0 <- 1765924  # Initial infected humans
R_H0 <- 1606991  # Initial recovered humans
S_V0 <- 1000000  # Initial susceptible mosquitoes (assume 1 million)
I_V0 <- 10000  # Initial infected mosquitoes (assume 10,000)

initial_state <- c(S_H = S_H0, I_H = I_H0, R_H = R_H0, S_V = S_V0, I_V = I_V0)

# Define parameters
parameters <- c(
  Lambda_H = 26182, # human recruitment rate
  Lambda_V = 0, # Assume constant mosquito population for simplicity
  beta_H = 0.0001,  # Transmission rate from mosquitoes to humans
  beta_V = 0.0002,  # Transmission rate from humans to mosquitoes
  gamma_H = 0.0714,  # Recovery rate in humans (1/14 days)
  mu_H = 0.0009,  # Human mortality rate
  mu_V = 0.1,  # Mosquito mortality rate
  psi = 0,  # ITN effectiveness
  sigma = 0,  # Antimalarial drug effectiveness
  theta_H = 0.001, # waning of immunity 
  delta_H = 0.09) #death induced


# Time frame for simulation (e.g., 1 year with daily steps)
time <- seq(0, 72, by = 1)

# Solve the differential equations
output <- ode(y = initial_state, times = time, func = sir_si_model, parms = parameters)

print(output)
# Convert output to a data frame
output_df <- as.data.frame(output)


# Plot the results
plot(output_df$time, output_df$S_H, type = "l", col = "blue", ylim = c(0, max(output_df$S_H, output_df$I_H, output_df$R_H)), 
     xlab = "Time (days)", ylab = "Population", main = "SIR-SI Model: Oyo with zero intervention")
lines(output_df$time, output_df$I_H, col = "red")
lines(output_df$time, output_df$R_H, col = "green")
lines(output_df$time, output_df$S_V, col = "orange")
legend("right", legend = c("Susceptible (Humans)", "Infected (Humans)", "Recovered (Humans)"), col = c("blue", "red", "green"), lty = 1)

# Plot Mosquito population
plot(output_df$time, output_df$S_V, type = "l", col = "orange", ylim = c(0, max(output_df$S_V, output_df$I_V)), 
     xlab = "Time (days)", ylab = "Mosquito Population", main = "SIR-SI Model: Kano Mosquito Population")
lines(output_df$time, output_df$I_V, col = "purple")
legend("right", legend = c("Susceptible (Mosquitoes)", "Infected (Mosquitoes)"), col = c("orange", "purple"), lty = 1)


############################################################################################
 # with 90% intervention
install.packages("deSolve")
library(deSolve)


# Define the SIR-SI model differential equations
sir_si_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Human equations
    dS_H <- (Lambda_H + theta_H) - (beta_H * (1 - psi) * S_H * I_V)/N_H - (mu_H * S_H)
    dI_H <- (beta_H * (1 - psi) * S_H * I_V)/N_H - (gamma_H * (1 + sigma) + mu_H + delta_H) * I_H
    dR_H <- gamma_H * (1 + sigma) * I_H - (mu_H + theta_H) * R_H
    
    # Mosquito equations
    dS_V <- Lambda_V - (beta_V * (1 - psi) * S_V * I_V)/N_H - (mu_V * S_V)
    dI_V <- (beta_V * (1 - psi) * S_V * I_V)/N_H  - (mu_V * I_V)
    
    # Return the rate of change
    list(c(dS_H, dI_H, dR_H, dS_V, dI_V))
  })
}

# Initial state values
N_H <-  7358018  # Total human population
S_H0 <- 3985103  # Initial susceptible humans
I_H0 <- 1765924  # Initial infected humans
R_H0 <- 1606991  # Initial recovered humans
S_V0 <- 1000000  # Initial susceptible mosquitoes (assume 1 million)
I_V0 <- 10000  # Initial infected mosquitoes (assume 10,000)

initial_state <- c(S_H = S_H0, I_H = I_H0, R_H = R_H0, S_V = S_V0, I_V = I_V0)

# Define parameters
parameters <- c(
  Lambda_H = 26182, # human recruitment rate
  Lambda_V = 0, # Assume constant mosquito population for simplicity
  beta_H = 0.0001,  # Transmission rate from mosquitoes to humans
  beta_V = 0.0002,  # Transmission rate from humans to mosquitoes
  gamma_H = 0.0714,  # Recovery rate in humans (1/14 days)
  mu_H = 0.0009,  # Human mortality rate
  mu_V = 0.1,  # Mosquito mortality rate
  psi = 0.9,  # ITN effectiveness
  sigma = 0.8,  # Antimalarial drug effectiveness
  theta_H = 0.001, # waning of immunity 
  delta_H = 0.09) #death induced


# Time frame for simulation (e.g., 1 year with daily steps)
time <- seq(0, 72, by = 1)

# Solve the differential equations
output <- ode(y = initial_state, times = time, func = sir_si_model, parms = parameters)

print(output)
# Convert output to a data frame
output_df <- as.data.frame(output)


# Plot the results
plot(output_df$time, output_df$S_H, type = "l", col = "blue", ylim = c(0, max(output_df$S_H, output_df$I_H, output_df$R_H)), 
     xlab = "Time (days)", ylab = "Population", main = "SIR-SI Model: Oyo with 90% intervention")
lines(output_df$time, output_df$I_H, col = "red")
lines(output_df$time, output_df$R_H, col = "green")
lines(output_df$time, output_df$S_V, col = "orange")
legend("right", legend = c("Susceptible (Humans)", "Infected (Humans)", "Recovered (Humans)"), col = c("blue", "red", "green"), lty = 1)

# Plot Mosquito population
plot(output_df$time, output_df$S_V, type = "l", col = "orange", ylim = c(0, max(output_df$S_V, output_df$I_V)), 
     xlab = "Time (days)", ylab = "Mosquito Population", main = "SIR-SI Model: Kano Mosquito Population")
lines(output_df$time, output_df$I_V, col = "purple")
legend("right", legend = c("Susceptible (Mosquitoes)", "Infected (Mosquitoes)"), col = c("orange", "purple"), lty = 1)


############################################################################################