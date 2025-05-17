# Algorithm 1

## Take "Variance of population mean with local epsilon-DP Laplace" as an example 

library(tictoc)

# N = 1000*(1:5)
# sigmaSq = 0.08^(1:5)
# ep = 1
# k = length(N)
# eta = 200

N = 1000*(20:11)
sigmaSq = 0.08^((11:20)/10)
ep = 1
k = length(N)
eta = 100000


tic()
Og.Var <- function(x,N,sigmaSq){
  sum((N^2/x)*sigmaSq / sum(N)^2)
}

Pure.Lap.Var <- function(x,N,sigmaSq,ep){
  sum((N^2/x)*2*(1/log(1+N*(exp(ep)-1)/x))^2 / sum(N)^2)
}

fn.Var <- function(x,N,sigmaSq,ep){
  Og.Var(x,N,sigmaSq) + Pure.Lap.Var(x,N,sigmaSq,ep)
}

###
# library(nloptr)
# # Define the objective function
# fun <- function(params) {
#   x <- params
#   # Define the original objective function
#   #return(fn.Var(x,N,sigmaSq,ep))
#   return(sum(N)^2*fn.Var(x,N,sigmaSq,ep)) # fn.Var times sum(N)^2 is easier for nloptr to run, avoiding rounding in decimals
# }
# 
# # Define the equality constraint function
# constraint_fun <- function(params) {
#   x <- params
#   return(eta - sum(x))
# }
# 
# # Define the initial guess for the parameters
# initial_guess <- rep(eta/k,k)
# lb = rep(1, k)
# ub = rep(eta, k)
# 
# # Set up the optimization using nloptr
# result_nloptr <- nloptr(
#   x0 = initial_guess,               # Starting point
#   eval_f = fun,                     # Objective function
#   lb = lb,                     # Lower bounds 
#   ub = ub,                 # Upper bounds
#   eval_g_eq = constraint_fun,     # Inequality constraints
#   opts = list(
#     "algorithm" = "NLOPT_LN_COBYLA",  # Algorithm that handles nonlinear constraints
#     "xtol_rel" = 1.0e-20
#   )
# )
# 
# print(result_nloptr)
# xStar = result_nloptr$solution
# var_xStar = fn.Var(xStar, N, sigmaSq, ep)

library(alabama)

# Define the objective function
fun <- function(params) {
  x <- params
  return(sum(N)^3.2 * fn.Var(x, N, sigmaSq, ep))
}

# Define constraints (equality and bounds as inequality constraints)
constraint_fun <- function(params) {
  c(eta - sum(params),  # Equality constraint: sum(x) = eta
    params - lb,        # Lower bound: x >= lb
    ub - params)        # Upper bound: x <= ub
}

constraint_fun1 <- function(params) {
  c(params - lb,        # Lower bound: x >= lb
    ub - params)        # Upper bound: x <= ub
}

constraint_fun2 <- function(params) {
  eta - sum(params)
}

# Define initial guess, bounds, and constraints
initial_guess <- rep(eta / k, k)
lb <- rep(1, k)
ub <- rep(eta, k)

# Optimization using alabama package
result_alabama <- auglag(
  par = initial_guess,
  fn = fun,
  hin = constraint_fun1,  # Inequality constraints (includes bounds and equality)
  heq = constraint_fun2,
  control.outer = list(method = "BFGS", maxit = 1000, tol = 1e-10)
)

xStar <- result_alabama$par
var_xStar <- fn.Var(xStar, N, sigmaSq, ep)
sum(xStar)

#

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
#H_f = diag( 2*N^2*( sigmaSq/x^3 + ((x+c)^2*(log(1+c/x))^2 - c*(4*x+3*c)*log(1+c/x) + 3*c^2)/ ((x^3*(x+c)^2)*(log(1+c/x))^4) ) / sum(N)^2 )
H_f_scaled = diag( 2*N^2*( sigmaSq/x^3 + ((x+c)^2*(log(1+c/x))^2 - c*(4*x+3*c)*log(1+c/x) + 3*c^2)/ ((x^3*(x+c)^2)*(log(1+c/x))^4) ) )
H_f_info = eigen(H_f_scaled)
H_f_eigenvalues <- H_f_info$values / sum(N)^2  # Rescale back
max_eigenvalues = max(H_f_eigenvalues)
min_eigenvalues = min(H_f_eigenvalues)


### radius 
r = sqrt(2*(var_n_init - var_xStar)/min_eigenvalues)


### Select integer optimal from set S
lowUpBound <- rbind(xStar - r, xStar + r)
lowUpBound <- ifelse(lowUpBound < 0, 0, lowUpBound)
lowUpBound <- ifelse(lowUpBound > eta, eta, lowUpBound)
comboS_butk <- lapply(1:(k-1), function(i){
  seq(ceiling(lowUpBound[1,i]), floor(lowUpBound[2,i]), by = 1)
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
toc()


result_alabama$par


