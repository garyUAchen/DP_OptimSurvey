library(tidyr)
library(dplyr)
library(ggplot2)

# Import data from Excel
data <- data.frame(
  eta = c(30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 100, 1000, 10000, 100000),
  grid = c(105.46, 256.68, 485.58, 948.27, 1615.97, 2150.22, 3554.8, 5567.57, 8715.83, 13387.384, 100000, 100000, 100000, 100000),
  algorithm = c(0.123, 0.16, 0.116, 0.08, 0.085, 0.086, 0.136, 0.225, 0.168, 0.147, 0.451, 0.22,	0.174,	0.473)
)

# Transform data into long format for ggplot
long_data <- pivot_longer(data, cols = c(grid, algorithm), names_to = "method", values_to = "value")


# Convert eta into a factor for custom x-axis, ensuring smaller values are on the left
long_data$eta_label <- factor(
  long_data$eta, 
  levels = c(as.character(sort(unique(long_data$eta[long_data$eta <= 48]))), "100", "1000", "10000", "100000")
)
long_data$eta_label[27] = "100000"
long_data$eta_label[28] = "100000"

ggplot(long_data, aes(x = eta_label, y = value, color = method, group = method)) +
  geom_line() +  # Draw lines between points
  geom_point() + # Add points for clarity
  scale_y_log10(
    limits = c(0.05, 15000),  # Set limits while preserving log scale
    breaks = c(1, 10, 100, 1000, 10000, 100000),  # Customize breaks
    labels = scales::label_comma()  # Format labels with commas
  ) +
  scale_color_manual(
    values = c("algorithm" = "#F8766D", "grid" = "#00BFC4"),
    labels = c("algorithm" = "Algorithm 1", "grid" = "Exhaustive")
  ) +
  labs(
    title = "",
    x = expression("Sampe Size"~eta),
    y = "Computation Time (sec.)",
    color = "Method"
  ) +
  geom_vline(xintercept = 10.5, linetype = "dashed", color = "black") + # Add a dashed line between 48 and 100
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 16),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
  )



