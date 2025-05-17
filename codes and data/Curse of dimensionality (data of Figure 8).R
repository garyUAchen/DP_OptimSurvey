# Load necessary libraries
library(ggplot2)
library(tidyr)

# Input the data
eta_values <- c(1000, 10000, 100000)
k_values1 <- c(2, 4, 6, 8, 10, 12, 14)
k_values2 <- c(14, 16, 18, 20, 22, 24, 26)

# Table 1 (Radius)
radius_values <- matrix(c(
  0.414, 0.789, 0.963, 1.006, 1.245, 1.505, 1.93,
  0.118, 0.252, 1.165, 1.293, 1.263, 1.399, 1.391,
  0.647, 0.759, 0.828, 1.236, 1.605, 1.44, 1.625
), nrow = length(eta_values), byrow = TRUE)

data_radius <- expand.grid(eta = eta_values, k = k_values1)
data_radius$value <- as.vector(radius_values)
data_radius$type <- "Radius"

# Table 2 (Gap)
gap_values <- matrix(c(
  1.71E-05, 3.00E-05, 1.97E-05, 3.75E-05, 4.29E-05, 4.04E-05, 4.37E-05,
  6.22E-08, 1.76E-07, 1.92E-07, 2.22E-07, 3.51E-07, 3.78E-07, 4.67E-07,
  2.96E-10, 4.82E-10, 2.71E-10, 9.73E-10, 1.53E-09, 1.27E-09, 1.99E-09
), nrow = length(eta_values), byrow = TRUE)

data_gap <- expand.grid(eta = eta_values, k = k_values2)
data_gap$value <- as.vector(gap_values)
data_gap$type <- "Gap"

# Table 3 (Computation Time)
computation_values <- matrix(c(
  0.305, 0.58, 2.194, 8.599, 34.943, 154.982, 648.92,
  0.224, 0.608, 2.344, 8.213, 35, 162.2, 661.46,
  0.176, 0.505, 2.184, 8.585, 35.467, 164.407, 654.94
), nrow = length(eta_values), byrow = TRUE)

data_computation <- expand.grid(eta = eta_values, k = k_values2)
data_computation$value <- as.vector(computation_values)
data_computation$type <- "Computation Time"

# Combine all data frames
data_combined <- rbind(data_radius, data_gap, data_computation)

# Reorder factor levels for desired left-to-right order
levels_order <- c("Radius", "Computation Time", "Gap")
data_combined$type <- factor(data_combined$type, levels = levels_order)

# Create the combined plot
plot <- ggplot(data_combined, aes(x = k, y = value, color = factor(eta), group = eta)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ type, scales = "free", ncol = 3) +
  labs(
    x = "Number of Groups (k)",
    y = "Value (log scale)",
    color = expression("Sample size")
  ) +
  scale_y_log10() +
  scale_x_continuous(n.breaks = 8) + 
  theme_minimal()

# Create a named vector for custom titles
custom_titles <- c(
  "Radius" = "Radius (r)",
  "Computation Time" = "Runtime (sec.)",
  "Gap" = "Optimality Gap"
)

# Use the labeller argument in facet_wrap
plot <- ggplot(data_combined, aes(x = k, y = value, color = factor(eta), group = eta)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ type, scales = "free", ncol = 3, labeller = labeller(type = custom_titles)) +
  labs(
    x = "Number of Groups (k)",
    y = "Value (log-scaled y-axis)",
    color = expression("Sample Size")
  ) +
  scale_y_log10() +
  scale_x_continuous(n.breaks = 8) + 
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

# Display the plot
plot


