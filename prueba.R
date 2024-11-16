# Load required libraries upfront
required_packages <- c("ALA", "dplyr", "purrr", "hypergeo", "ggplot2", 
                      "parallel", "gtools", "GA", "future", "rstan", "furrr")
invisible(lapply(required_packages, library, character.only = TRUE))

# Set environment variables and options
Sys.setenv("R_MAX_NUM_DLL" = 9999999)
options(mc.cores = parallel::detectCores())

# Initialize parameters
seed <- 123
n_bootstrap <- 1000
N1 <- N2 <- 500

# Helper function to process data
process_school_data <- function(datos1, school_code) {
  datos1 %>%
    mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
    filter(school == school_code)
}

# Function to calculate summary statistics
get_summary_stats <- function(data, school_based, tv_based) {
  filtered_data <- data %>%
    filter(school.based == school_based & tv.based == tv_based)
  
  X <- filtered_data %>%
    group_by(stage) %>%
    summarise(Bin = sum(binTHKS))
  
  n <- filtered_data %>%
    group_by(stage) %>%
    summarise(n = n())
  
  list(X = X, n = n)
}

# Process all school combinations efficiently
school_combinations <- list(
  yy = list(code = "404", school = "yes", tv = "yes"),
  yn = list(code = "408", school = "yes", tv = "no"),
  ny = list(code = "508", school = "no", tv = "yes"),
  nn = list(code = "409", school = "no", tv = "no")
)

# Process all schools in parallel
school_results <- future_map(school_combinations, function(combo) {
  datos <- process_school_data(datos1, combo$code)
  get_summary_stats(datos, combo$school, combo$tv)
})

# Optimized density functions
densBB_H_fucntor <- function(x1, x2, n1, n2, u, l, a, a0, a1, a2, i, use_log = FALSE) {
  function(theta1) {
    log_f <- lgamma(n1 + a) + lgamma(n2 + a) -
      lgamma(x1[i] + a1) - lgamma(n1 - x1[i] + a - a1) -
      lgamma(x2[i] + a2) - lgamma(n2 - x2[i] + a - a2) -
      log(genhypergeo(U = u, L = l, check_mod = TRUE, z = 1)) +
      (a1 + x1[i] - 1) * log(theta1) +
      (a2 + a0 + (n1 - x1[i]) - 1) * log(1 - theta1) +
      (a2 + x2[i] - 1) * log(theta1) +
      (a1 + a0 + (n2 - x2[i]) - 1) * log(1 - theta1) -
      a * log(1 - theta1 * theta1)
    
    if (!use_log) exp(log_f) else log_f
  }
}

# Vectorized binomial function
f2 <- Vectorize(function(X1, X2, theta1) {
  exp(dbinom(X1, n1, theta1, log = TRUE) + 
      dbinom(X2, n2, theta1, log = TRUE))
})

# Set prior parameters
alphas_opt <- c(0.8373879, 0.8410984, 0.8053298)
a1 <- alphas_opt[2]
a2 <- alphas_opt[3]
a0 <- alphas_opt[1]
a <- sum(alphas_opt)

# Compile Stan model once
stan_model <- rstan::stan_model(file = 'BBpost3.stan', verbose = TRUE)

# Generate random variables for Monte Carlo integration
set.seed(seed)
Theta1 <- rbeta(N2, 1, 1)

