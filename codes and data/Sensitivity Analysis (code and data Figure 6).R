# Interplay between original variance and purely Laplace noise designs

## Section 5.2 (Figure 3)

library(ggplot2)
library(tidyr)

N = 1000*(1:3)
sigmaSq = 0.08^(c(1,1.5,2))
ep = 1 # variable
k = length(N)

eta = 200

Og.Var <- function(x,N,sigmaSq){
  sum((N^2/x)*sigmaSq / sum(N)^2)
}

og.alg <- function(ep){
  ogSolution = sqrt(sigmaSq)*N / sum(sqrt(sigmaSq)*N)*eta
  comboA_butk = as.list(as.data.frame(rbind(floor(ogSolution[-k]), ceiling(ogSolution[-k]))))
  comboA_butk = expand.grid(comboA_butk)
  
  comboA <- t(apply(comboA_butk, MARGIN = 1, function(y) {
    xk = eta - sum(y)
    x = c(y, xk)
    if (min(x) > 0 & max(x) < (eta-k+1)) {
      return(x)
    } else {
      return(NULL)  # Exclude this combination
    }
  }))
  
  adjust.output <- function(output) {
    if (is.list(output)) {
      output <- do.call(rbind, output)  
    }
    return(output)
  }
  comboA = adjust.output(comboA)
  
  var_value = apply(comboA, MARGIN = 1, function(n){
    Og.Var(n,N,sigmaSq)
  })
  ogSolution = comboA[which(var_value == min(var_value)),]
  var_ogSolution = Og.Var(ogSolution,N,sigmaSq)
  
  return(ogSolution)
}


Pure.Lap.Var <- function(x,N,sigmaSq,ep){
  sum((N^2/x)*2*(1/log(1+N*(exp(ep)-1)/x))^2 / sum(N)^2)
}

pureLap.alg <- function(ep){
  pureLapSolution = N / sum(N)*eta
  comboA_butk = as.list(as.data.frame(rbind(floor(pureLapSolution[-k]), ceiling(pureLapSolution[-k]))))
  comboA_butk = expand.grid(comboA_butk)
  
  comboA <- t(apply(comboA_butk, MARGIN = 1, function(y) {
    xk = eta - sum(y)
    x = c(y, xk)
    if (min(x) > 0 & max(x) < (eta-k+1)) {
      return(x)
    } else {
      return(NULL)  # Exclude this combination
    }
  }))
  
  adjust.output <- function(output) {
    if (is.list(output)) {
      output <- do.call(rbind, output)  
    }
    return(output)
  }
  comboA = adjust.output(comboA)
  
  var_value = apply(comboA, MARGIN = 1, function(n){
    Pure.Lap.Var(n,N,sigmaSq,ep)
  })
  pureLapSolution = comboA[which(var_value == min(var_value)),]
  var_pureLapSolution = Pure.Lap.Var(pureLapSolution,N,sigmaSq,ep)
  
  return(pureLapSolution)
}



library(nloptr)

fn.Var <- function(x,N,sigmaSq,ep){
  Og.Var(x,N,sigmaSq) + Pure.Lap.Var(x,N,sigmaSq,ep)
}

true.alg <- function(ep){
  # Define the objective function
  fun <- function(params) {
    x <- params
    # Define the original objective function
    #return(fn.Var(x,N,sigmaSq,ep))
    return(sum(N)^2*fn.Var(x,N,sigmaSq,ep)) # fn.Var times sum(N)^2 is easier for nloptr to run, avoiding rounding in decimals
  }
  
  # Define the equality constraint function
  constraint_fun <- function(params) {
    x <- params
    return(eta - sum(x))
  }
  
  # Define the initial guess for the parameters
  initial_guess <- rep(eta/k,k)
  lb = rep(1, k)
  ub = rep(eta, k)
  
  # Set up the optimization using nloptr
  result_nloptr <- nloptr(
    x0 = initial_guess,               # Starting point
    eval_f = fun,                     # Objective function
    lb = lb,                     # Lower bounds 
    ub = ub,                 # Upper bounds
    eval_g_eq = constraint_fun,     # Inequality constraints
    opts = list(
      "algorithm" = "NLOPT_LN_COBYLA",  # Algorithm that handles nonlinear constraints
      "xtol_rel" = 1.0e-20
    )
  )
  
  print(result_nloptr)
  xStar = result_nloptr$solution
  var_xStar = fn.Var(xStar, N, sigmaSq, ep)
  
  
  comboT_butk = as.list(as.data.frame(rbind(floor(xStar[-k]), ceiling(xStar[-k]))))
  comboT_butk = expand.grid(comboT_butk)
  
  comboT <- t(apply(comboT_butk, MARGIN = 1, function(y) {
    xk = eta - sum(y)
    x = c(y, xk)
    if (min(x) > 0 & max(x) < (eta-k+1)) {
      return(x)
    } else {
      return(NULL)  # Exclude this combination
    }
  }))
  
  adjust.output <- function(output) {
    if (is.list(output)) {
      output <- do.call(rbind, output)  
    }
    return(output)
  }
  comboT = adjust.output(comboT)
  
  ### Select n_init from set T
  var_value = apply(comboT, MARGIN = 1, function(n){
    fn.Var(n,N,sigmaSq,ep)
  })
  n_init = comboT[which(var_value == min(var_value)),]
  var_n_init = fn.Var(n_init,N,sigmaSq,ep)
  
  ### eigenvalues
  c = (exp(ep)-1)*N
  # x components with Laplace noise
  x = xStar
  H_f = diag( 2*N^2*( sigmaSq/x^3 + ((x+c)^2*(log(1+c/x))^2 - c*(4*x+3*c)*log(1+c/x) + 3*c^2)/ ((x^3*(x+c)^2)*(log(1+c/x))^4) ) / sum(N)^2 )
  H_f_info = eigen(H_f)
  H_f_eigenvalues = H_f_info$values #ifelse statement for negative eigenvalues
  max_eigenvalues = max(H_f_eigenvalues)
  min_eigenvalues = min(H_f_eigenvalues)
  
  
  ### radius 
  r = sqrt(2*(var_n_init - var_xStar)/min_eigenvalues)
  
  ### Select integer optimal from set S
  lowUpBound <- rbind(xStar - r, xStar + r)
  lowUpBound <- ifelse(lowUpBound < 0, 0, lowUpBound)
  comboS_butk <- lapply(1:(k-1), function(i){
    seq(ceiling(lowUpBound[1,i]), floor(lowUpBound[2,i]))
  })
  comboS_butk = expand.grid(comboS_butk)
  comboS <- t(apply(comboS_butk, MARGIN = 1, function(y) {
    xk = eta - sum(y)
    x = c(y, xk)
    if (min(x) > 0 & max(x) < (eta-k+1)) {
      return(x)
    } else {
      return(NULL)  # Exclude this combination
    }
  }))
  comboS = adjust.output(comboS)
  
  
  ### Select nStar from set S
  var_values = apply(comboS, MARGIN = 1, function(n){
    fn.Var(n,N,sigmaSq,ep)
  })
  nStar = comboS[which(var_values == min(var_values)),]
  var_nStar = fn.Var(nStar,N,sigmaSq,ep)
  
  return(nStar)
}

### 300% ###
sigmaSq = 0.08^(c(1,1.5,2))*3
temp3 = sapply(10^seq(-2,2,0.25), function(i){
  c(i, true.alg(i))
})

### 120% ###
sigmaSq = 0.08^(c(1,1.5,2))*1.2
temp1.2 = sapply(10^seq(-2,2,0.25), function(i){
  c(i, true.alg(i))
})

### 100% ###
sigmaSq = 0.08^(c(1,1.5,2))
temp1 = sapply(10^seq(-2,2,0.25), function(i){
  c(i, true.alg(i))
})

### 90% ###
sigmaSq = 0.08^(c(1,1.5,2))*0.9
temp0.9 = sapply(10^seq(-2,2,0.25), function(i){
  c(i, true.alg(i))
})

### 30% ###
sigmaSq = 0.08^(c(1,1.5,2))*0.3
temp0.3 = sapply(10^seq(-2,2,0.25), function(i){
  c(i, true.alg(i))
})






# Convert to data frames for easier manipulation
data1 <- data.frame(
  epsilon = temp1[1, ],
  solution1 = temp3[2, ],
  solution2 = temp1.2[2, ],
  solution3 = temp1[2, ],
  solution4 = temp0.9[2, ],
  solution5 = temp0.3[2, ]
)

data2 <- data.frame(
  epsilon = temp1[1, ],
  solution1 = temp3[3, ],
  solution2 = temp1.2[3, ],
  solution3 = temp1[3, ],
  solution4 = temp0.9[3, ],
  solution5 = temp0.3[3, ]
)
data3 <- data.frame(
  epsilon = temp1[1, ],
  solution1 = temp3[4, ],
  solution2 = temp1.2[4, ],
  solution3 = temp1[4, ],
  solution4 = temp0.9[4, ],
  solution5 = temp0.3[4, ]
)

# data2 <- data.frame(
#   epsilon = temp[1, ],
#   solution1 = tempPlus[3, ],
#   solution2 = tempNegative[3, ],
#   solution3 = temp[3, ]
# )
# 
# data3 <- data.frame(
#   epsilon = temp[1, ],
#   solution1 = tempPlus[4, ],
#   solution2 = tempNegative[4, ],
#   solution3 = temp[4, ]
# )

# Add a column to identify the groups
data1$group <- "Group 1"
data2$group <- "Group 2"
data3$group <- "Group 3"

# Combine the datasets
combined_data <- rbind(data1, data2, data3)

# Reshape the combined dataset
long_data <- pivot_longer(
  combined_data,
  cols = starts_with("solution"),
  names_to = "solution",
  values_to = "value"
)


# Plot the data
ggplot(long_data, aes(x = epsilon, y = value, color = solution)) +
  geom_line(size = 1) + # Add lines for each solution
  geom_point() +        # Add points to highlight data
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),  # Specify the breaks
                labels = c("0.01", "0.1", "1", "10", "100")) +     # Use logarithmic scale for x-axis
  scale_color_manual(
    values = c("solution1" = "#003FAD", "solution2" = "#2C6EFF", "solution3" = "#619CFF", "solution4" = "#88B4FF", "solution5" = "#AFCBFF" ),  # Keep colors
    labels = c("solution1" = "300%", "solution2" = "120%", "solution3" = "100%", "solution4" = "90%", "solution5" = "30%")  # Change legend labels
  ) +
  labs(
    title = "",
    x = expression(epsilon),
    y = "Subsampling Sizes",
    color = "Integer Solution"
  ) +
  facet_wrap(~group, scales = "free_x") + # Create separate panels for each group
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 41, hjust = 1, size = 13),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 17),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 17),  # facet title font size
  )





