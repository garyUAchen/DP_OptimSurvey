library(tictoc)
### Variance of population average with local epsilon-DP Laplace

# Algorithm

k = 14
eta = 100000

end = 10+k

N = 1000*(end:11)
sigmaSq = 0.08^((11:end)/10)
ep = 1



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

library(alabama)

# Define the objective function
fun <- function(params) {
  x <- params
  return(sum(N)^3 * fn.Var(x, N, sigmaSq, ep))
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
result_alabama <- alabama::auglag(
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

var_n_init/var_xStar - 1

toc()


