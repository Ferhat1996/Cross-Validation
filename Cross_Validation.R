
#**************************Nonparametric-Econometrics**************************#

#*******************************Ferhat-Kutluca*********************************#

#*******************************Problem-Set-3**********************************#

#*********************************Exercise_3***********************************#

# Install required packages
# install.packages(c("np","ICV"))

# Load required libraries
library(np)
library(ICV)

#******************************************************************************#

# Load the "faithful" dataset
data(faithful)

#****************Briefly describe the data and state the source****************#


#The data used in the provided R script comes from the "faithful" dataset, which is included in the base R installation.
#This dataset is a well-known example in statistical analysis and is often used for educational purposes.
#The "faithful" dataset contains information about the eruption duration and waiting time between eruptions of the Old 
#Faithful geyser in Yellowstone National Park. Each observation corresponds to a recorded eruption event, with two variables:

#eruptions: The duration of the eruption (in minutes).
#waiting: The waiting time until the next eruption (in minutes).

#This dataset is frequently used to illustrate concepts related to statistical analysis, including kernel density estimation,
#as demonstrated in the provided R script. The script utilizes the "faithful" dataset to compare different  methods for
#bandwidth selection in kernel density estimation, employing custom functions as well as packages like 'np' and 'ICV'.#

#**************************************(i)*************************************#

# (i) Custom method for optimal bandwidth estimation

# Define the Gaussian kernel function
kernel <- function(v) {
  (1 / sqrt(2 * pi)) * exp(-0.5 * v^2)
}

# Define the Gaussian overlined_kernel function which is the twofold convolution kernel derived from k(.)
overlined_kernel <- function(v) {
  exp(-v^2 / 4) / sqrt(4 * pi)
}

# Define the modified cross-validation criterion function
cv_criterion <- function(h, data) {
  n <- length(data)
  total_sum_1 <- 0
  total_sum_2 <- 0
  
  for (i in 1:n) {
    for (j in 1:n) {
      total_sum_1 <- total_sum_1 + overlined_kernel((data[i] - data[j]) / h)
    }
  }
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        total_sum_2 <- total_sum_2 + kernel((data[i] - data[j]) / h)
      }
    }
  }
  
  return(total_sum_1 / (n^2 * h) - 2 * total_sum_2 / (n * (n - 1) * h))
}

#******************************************************************************#

# Find the optimal bandwidth using the optimize function
initial_guess <- 2  # Adjust the initial guess accordingly
optimal_bandwidth <- optimize(cv_criterion, interval = c(0.001, 10), data = faithful$waiting)$minimum

#************************************(ii)**************************************#

# (ii) Applying the least squares cross-validation method to estimate the density
# Kernel density estimation using the optimal bandwidth
density_estimate_custom <- density(faithful$waiting, bw = optimal_bandwidth)

# Plotting the results
plot(density_estimate_custom, main = "Custom Method", col = "blue")

# Print results for (ii)
cat("Custom Method - Optimal Bandwidth:", optimal_bandwidth, "\n")

#***********************************(iii.1)************************************#

# (iii.1) Using the 'np' package for kernel density estimation with cross-validation

# Kernel density estimation using the 'np' package with cross-validation
# npudens Nonparametric Density Estimation
# npudensbw Nonparametric Density Bandwidth Selection
# Reference: "chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/np/vignettes/np.pdf"

density_estimate_np <- npudens(faithful$waiting)

# Plotting the result
plot(density_estimate_np, main = "NP Package", col = "red")
mtext(paste("Bandwidth: ", round(density_estimate_np$bw, 3)), side = 1, line = 2)

# Print results for (iii.1)
cat("np Package - Bandwidth:", density_estimate_np$bw, "\n")

#************************************(iii.2)***********************************#

# (iii.2) Using the 'ICV' package for kernel density estimation
# Kernel density estimation using the 'ICV' package
# KDE_ICV Computing the Gaussian density estimate based on h_ICV
# h_ICV Calculation of the ICV bandwidth for the Gaussian density estimator corresponding to expression (12) of Savchuk, Hart, and Sheather(2010).
# Reference: "chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/ICV/ICV.pdf"

# Plotting the result
density_estimate_icv <- KDE_ICV(faithful$waiting)

density_estimate_icv_Bandwidth <- h_ICV(faithful$waiting)

# Print results for (iii.2)
cat("ICV Package - Bandwidth:", density_estimate_icv_Bandwidth, "\n")

#******************************************************************************#

















