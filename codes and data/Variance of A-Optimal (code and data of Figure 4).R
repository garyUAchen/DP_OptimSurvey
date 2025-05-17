# variance of A-optimal

## Section 5.1 (Figure 6 in appendix B)

library(ggplot2)
library(tidyr)

N = 1000*(7:10)
sigmaSq = 0.08^(1:4)
k = length(N)
eta = 200
ep=1

Og.Var <- function(x,N,sigmaSq){
  sum((1/x)*sigmaSq)
}

Pure.Lap.Var <- function(x,N,sigmaSq,ep){
  sum((1/x)*2*(1/log(1+N*(exp(ep)-1)/x))^2)
}

fn.Var <- function(x,N,sigmaSq,ep){
  Og.Var(x,N,sigmaSq) + Pure.Lap.Var(x,N,sigmaSq,ep)
}


naive.alg <- function(ep){
  naiveSolution = sqrt(sigmaSq)*N / sum(sqrt(sigmaSq)*N)*eta
  comboA_butk = as.list(as.data.frame(rbind(floor(naiveSolution[-k]), ceiling(naiveSolution[-k]))))
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
    fn.Var(n,N,sigmaSq,ep)
  })
  naiveSolution = comboA[which(var_value == min(var_value)),]
  var_naiveSolution = fn.Var(naiveSolution,N,sigmaSq,ep)
  
  return(var_naiveSolution)
}


library(nloptr)

true.alg <- function(ep){
  # Define the objective function
  fun <- function(params) {
    x <- params
    # Define the original objective function
    #return(fn.Var(x,N,sigmaSq,ep))
    return(fn.Var(x,N,sigmaSq,ep)) # fn.Var times sum(N)^2 is easier for nloptr to run, avoiding rounding in decimals
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
  H_f = diag( 2*( sigmaSq/x^3 + ((x+c)^2*(log(1+c/x))^2 - c*(4*x+3*c)*log(1+c/x) + 3*c^2)/ ((x^3*(x+c)^2)*(log(1+c/x))^4) ) )
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
  
  return(var_nStar)
}

df1 <- sapply(10^seq(-2,2,0.25), function(ep){
  
  c(ep, naive.alg(ep) / true.alg(ep))
  
})

df1 <- as.data.frame(t(df1))
colnames(df1) <- c("epsilon", "ratio")
ggplot(df1, aes(x = epsilon, y = ratio)) +
  geom_line(color = "blue") +  # Line plot
  geom_point(color = "blue") +  # Points for clarity
  scale_x_log10() +  # Log scale for x-axis (epsilon)
  labs(
    title = "Ratio vs. Epsilon on A-optimal with Laplace noise",
    x = expression(epsilon),
    y = "Ratio"
  ) +
  theme_minimal() 


####################################################################################################################################

Og.Var <- function(x,N,sigmaSq){
  sum((1/x)*sigmaSq)
}

Pure.TuLap.Var <- function(x,N,sigmaSq,ep){
  sum(2*(N*(exp(ep)-1)+x)/(exp(ep)-1)^2 ) + sum((1/x)*(1/12))
}

fn.Var <- function(x,N,sigmaSq,ep){
  Og.Var(x,N,sigmaSq) + Pure.TuLap.Var(x,N,sigmaSq,ep)
}

true.alg <- function(ep){
  # Define the objective function
  fun <- function(params) {
    x <- params
    # Define the original objective function
    #return(fn.Var(x,N,sigmaSq,ep))
    return(fn.Var(x,N,sigmaSq,ep)) # fn.Var times sum(N)^2 is easier for nloptr to run, avoiding rounding in decimals
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
  H_f = diag( 2*( (sigmaSq + 1/12) /x^3 ) )
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
  
  return(var_nStar)
}


df2 <- sapply(10^seq(-2,2,0.25), function(ep){
  
  c(ep, naive.alg(ep) / true.alg(ep))
  
})

df2 <- as.data.frame(t(df2))
colnames(df2) <- c("epsilon", "ratio")
ggplot(df2, aes(x = epsilon, y = ratio)) +
  geom_line(color = "blue") +  # Line plot
  geom_point(color = "blue") +  # Points for clarity
  scale_x_log10() +  # Log scale for x-axis (epsilon)
  labs(
    title = "Ratio vs. Epsilon on A-optimal with Tulap noise",
    x = expression(epsilon),
    y = "Ratio"
  ) +
  theme_minimal() 

####################################################################################################################################

Og.Var <- function(x,N,sigmaSq){
  sum((1/x)*sigmaSq)
}

Pure.DLap.Var <- function(x,N,sigmaSq,ep){
  sum( 2*(N*(exp(ep)-1)+x) / (exp(ep)-1)^2 ) 
}

fn.Var <- function(x,N,sigmaSq,ep){
  Og.Var(x,N,sigmaSq) + Pure.DLap.Var(x,N,sigmaSq,ep)
}

true.alg <- function(ep){
  # Define the objective function
  fun <- function(params) {
    x <- params
    # Define the original objective function
    #return(fn.Var(x,N,sigmaSq,ep))
    return(fn.Var(x,N,sigmaSq,ep)) # fn.Var times sum(N)^2 is easier for nloptr to run, avoiding rounding in decimals
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
  
  #print(result_nloptr)
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
  H_f = diag( 2*( sigmaSq/x^3 + 0.001 ) )
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
  
  return(var_nStar)
}



df3 = sapply(10^seq(-2,2,0.25), function(ep){
  print(ep)
  c(ep, naive.alg(ep) / true.alg(ep))
  
})

df3 <- as.data.frame(t(df3))
colnames(df3) <- c("epsilon", "ratio")
ggplot(df3, aes(x = epsilon, y = ratio)) +
  geom_line(color = "blue") +  # Line plot
  geom_point(color = "blue") +  # Points for clarity
  scale_x_log10() +  # Log scale for x-axis (epsilon)
  labs(
    title = "Ratio vs. Epsilon on A-optimal with DLap noise",
    x = expression(epsilon),
    y = "Ratio"
  ) +
  theme_minimal() 

####################################################################################################

df <- cbind(df1, df2[,2], df3[,2])
colnames(df) <- c("epsilon", "Laplace","TuLap","DLap")

# Reshape data into long format
df_long <- pivot_longer(
  df,
  cols = c("Laplace", "TuLap","DLap"),
  names_to = "Method",
  values_to = "Value"
)

library(ggplot2)


ggplot(df_long, aes(x = epsilon, y = Value, color = Method)) +
  geom_line(size = 1) +
  geom_point() +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),  # Specify the breaks
                labels = c("0.01", "0.1", "1", "10", "100")) +
  labs(
    title = "",
    x = expression(epsilon),
    y = "Variance Ratio",
    color = "Mechanism"
  ) +
  facet_wrap(~Method, scales = "free_y") +  # Separate panels for each method with independent y-scales
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 41, hjust = 1, size = 13),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 17),  # facet title font size
  )

