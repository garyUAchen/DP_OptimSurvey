library(tictoc, lib.loc="~/R/x86_64-pc-linux-gnu-library/4.2")

N = 1000*(20:11)
sigmaSq = 0.08^((11:20)/10)
ep = 1
k = length(N)

### total subsampling size: variable in simulation
eta <- as.numeric(Sys.getenv("eta"))
###

Og.Var <- function(x,N,sigmaSq){
  sum((N^2/x)*sigmaSq / sum(N)^2)
}

Pure.Lap.Var <- function(x,N,sigmaSq,ep){
  sum((N^2/x)*2*(1/log(1+N*(exp(ep)-1)/x))^2 / sum(N)^2)
}

fn.Var <- function(x,N,sigmaSq,ep){
  Og.Var(x,N,sigmaSq) + Pure.Lap.Var(x,N,sigmaSq,ep)
}

#################################################################################
tic()
find_min_variance_combination <- function(eta, k) {
  # Initialize variables to store the minimum variance and corresponding combination
  min_variance <- Inf
  best_combination <- NULL
  
  # Recursive function to generate combinations and calculate variance
  generate_combination <- function(sum_so_far, count, combination) {
    if (count == k - 1) {
      # The last element is determined by the remaining sum
      last_value <- eta - sum_so_far
      if (last_value > 0) {
        current_combination <- c(combination, last_value)
        
        # Calculate the variance for the current combination
        current_variance <- fn.Var(x = current_combination,N=N,sigmaSq=sigmaSq,ep=ep)
        
        # Update minimum variance and best combination if needed
        if (current_variance < min_variance) {
          min_variance <<- current_variance
          best_combination <<- current_combination
        }
      }
    } else {
      # Iterate over possible values for the current position
      for (i in 1:(eta - sum_so_far - (k - count - 1))) {
        generate_combination(sum_so_far + i, count + 1, c(combination, i))
      }
    }
  }
  
  # Start recursion with initial parameters
  generate_combination(0, 0, c())
  
  # Return the combination with the smallest variance
  list(best_combination = best_combination, min_variance = min_variance)
}

# Example usage

result <- find_min_variance_combination(eta=eta, k=k)
print(result$best_combination)
print(result$min_variance)
toc()


