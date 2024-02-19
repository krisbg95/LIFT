# Load necessary libraries
library(ggplot2)
library(dplyr) 
library(ggpubr)

# Generate example data
set.seed(123)
x <- seq(0, 10, by = 0.5)
y <- 2*x^2 - 3*x + rnorm(length(x), mean = 0, sd = 5)

# Create a data frame
data <- data.frame(x = x, y = y)

# Define a function to plot the data and fit a model
plot_and_fit <- function(data, degree, title) {
  # Fit a polynomial regression model
  model <- lm(y ~ poly(x, degree), data = data)
  
  # Create a data frame for predictions
  new_data <- data.frame(x = seq(0, 10, length.out = 100))
  new_data$y_pred <- predict(model, newdata = new_data)
  
  # Plot the data and predictions
  p <- ggplot(data, aes(x = x, y = y)) + 
    ylim(-5, 170) + 
    xlim(0,10)+
    geom_point(color = "#00AFBB", size = 2) +
    geom_line(data = new_data, aes(x = x, y = y_pred), color = "#E7B800",size = 1) +
    ggtitle(title) +
    theme_minimal()
  
  return(p)
}

# Plot underfitting, overfitting, and regular fit
p_underfit <- plot_and_fit(data, degree = 1, title = "Underfit")
p_overfit <- plot_and_fit(data, degree = 10, title = "Overfit")
p_regular_fit <- plot_and_fit(data, degree = 2, title = "Regular Fit")
 
ggarrange(p_underfit,p_overfit,p_regular_fit,ncol=3) 
dev.print(svg, "Bias_Variance_tradeoff.svg",height = 4, width = 7)
