rm(list=ls())

# Make results reproducible
set.seed(12345)
#------------------------------------------------

# Part 1
# For this Part, load in the trees data set from the datasets package.
library(datasets)
head(trees)
library(mvtnorm)

# Define criterion to be minimised in Gaussian process regression
gp_criterion = function(p,x,y) {
  sig_sq = exp(p[1])
  rho_sq = exp(p[2])
  tau_sq = exp(p[3])
  Mu = rep(0, length(x)) #The starting point of the optimisation
  Sigma = sig_sq * exp( - rho_sq * outer(x, x, '-')^2 ) + tau_sq * diag(length(x))
  ll = dmvnorm(y, Mu, Sigma, log = TRUE)
  return(-ll)
}

regression_fit <- function(x_g, x, y, p= 1, method= "BFGS"){
  #an integer >=0, indicates the order of the polynomial regression (defaults to 1)
  #the algorithm used to optimise the parameters of the Gaussian process
  #regression (defaults to ‘BFGS’)

  #The function should fit a polynomial regression with no intercept term and powers of
  #the covariate up to p; also, it should fit a Gaussian process regression with parameters
  #optimised using the method indicated. The starting point of the optimisation should be
  #c(0,0,0) in the log scale.

  #defined variables
  x_rep <- matrix(rep(x, p+1), ncol = p+1, nrow = length(x))
  X <- sweep(x_rep, 2, 0:p, "^")
  fun <- lsfit(X, y, intercept = FALSE) #no intercept term
  X_g_rep <- matrix(rep(x_g, p+1), ncol = p+1, nrow = length(x_g))
  X_g <- sweep(X_g_rep, 2, 0:p, "^")
  pred1<-X_g%*%fun$coefficients
  X <- sweep(x_rep, 2, 0:p, "^")

  if(method == "NR"){
    # The NR method is implemented in R's function nlminb
    result <- nlminb(start = rep(0, 3), x= x,y= y, objective = gp_criterion)
    print(gp_criterion(result$par, x, y))}

  # 1) BFGS - a fancier Newton-type method
  else if(method == "BFGS"){
    result <- optim(rep(0, 3), x= x, y= y, gp_criterion, method = "BFGS")
    print(gp_criterion(result$par, x, y))}

  # 2) Nelder-Mead - a simplex/heuristic method
  else if(method == "NM"){
    result <- optim(rep(0, 3), x= x, y= y, gp_criterion, method = "Nelder-Mead")
    print(gp_criterion(result$par, x, y))}

  # 3) Simulated annealing - a cool method based on thermodynamics - requires no
  #                          derivatives
  else if(method == "SANN"){
    result <- optim(rep(0, 3), x= x, y= y, gp_criterion, method = 'SANN')
    print(gp_criterion(result$par,x,y))}
    sig_sq <- exp(result$par[1])
    rho_sq <- exp(result$par[2])
    tau_sq <- exp(result$par[3])

  # Create covariance matrices
  C <- sig_sq * exp( - rho_sq * outer(x_g, x, "-")^2 )
  Sigma <- sig_sq * exp( - rho_sq * outer(x, x, "-")^2 ) + tau_sq * diag(length(x))

  # Now create predictions
  pred2 <- C %*% solve(Sigma, y)
  output<-data.frame("Polynomial" = pred1,"Gaussian" = pred2)
  return(output)}

#The function returns a list of two vectors that contain the predicted values using the
#polynomial regression (pred1) and Gaussian process regression (pred2), respectively.

# Create the plot here.
x <- scale(trees$Girth)[,1] #a vector of scaled covariates
y <- scale(trees$Height)[,1] #a vector of scaled responses
x_g <- pretty(x, n = 100) #a vector of scaled predictors

#Use an order of polynomial 4 and using the BFGS optimisation method for the GP regression model.
fun<-regression_fit(x_g, x, y, p=4, method="BFGS")

#Create a scatterplot of the data along with the fitted models and ensure this is appropriately
#labelled and includes a legend. Save this as ‘regression.pdf’.

#open pdf

pdf("regression.pdf", width = 12, height = 8)

#create a plot
plot(x,y, main = "Polynomial vs Gaussian",ylab="Height (in feet)",xlab = "Girth (in inches)",)
lines(x_g, fun$Polynomial, col = 'coral', lty = 'longdash', lwd=2)
lines(x_g, fun$Gaussian, col = 'blue', lty = 'longdash', lwd=2)
legend("bottomright", c("Polynomial","Gaussian"),fill=c("coral","blue"))
#close the pdf file
dev.off()

#------------------------------------------------

# Part 2

# We will be using data from the nycflights13 package.
library(nycflights13)
head(flights)
library(magrittr)
library(ggplot2)
library(reshape2)
str(flights)

#------------------------------------------------
# P2.a)
#------------------------------------------------

# Create a new dataset 'flights_2' that contains only the flights from 'EWR' to 'LAX'.
# Recast the 'carrier' variable as a factor, with levels in the following order:
# 'UA', 'VX', 'AA'.

# Solution
flights_2 <-
  flights %>%
  subset(and(origin == 'EWR',dest== 'LAX'))
flights_2$carrier <- factor(flights_2$carrier, levels = c('UA', 'VX', 'AA'))

str(flights_2)

#------------------------------------------------
# P2.b)
#------------------------------------------------

# Create a barplot where the bars show the number of flights from 'EWR' to 'LAX' for
# each of the carriers.  Save the plot as 'plot_1.pdf".
# Solution

p <- ggplot(flights_2, aes(x=carrier, fill = carrier)) + geom_bar()

p <- p + theme_minimal() + # Nicer theme
  scale_y_continuous(breaks = seq(0, 4000, by = 500)) +
  ylab('Number of Flights') + # Proper axis label
  xlab('Carrier') + # Proper axis label
  ggtitle('Number of Flights from EWR to LAX for each Carrier') +
  theme(legend.position="none")  # Remove legend
p

ggsave(p, file = 'plot_1.pdf', width = 12, height = 8)

#------------------------------------------------
# P2.c)
#------------------------------------------------

# Calculate the average air time for each carrier for flights from 'EWR' to 'LAX'.

# Plot the estimated densities for each of the underlying empirical distributions
# (i.e. 1 figure with 3 continuous lines, each corresponding to a different carrier).
# Save the plot as "plot_2.pdf".

# Solution
my_mean <- flights_2 %>% aggregate(air_time ~ carrier, data = ., FUN = 'mean')

q <- ggplot(flights_2, aes(x= air_time, color= carrier)) +
  geom_density() + ggtitle('Estimated Densities for each of Carrier') +
  theme(legend.position= "bottom") +
  ylab('Density') + # Proper axis label
  xlab('Air Time') # Proper axis label

q

ggsave(q, file = 'plot_2.pdf', width = 12, height = 8)

#------------------------------------------------
# P2.d)
#------------------------------------------------

# When producing the plot for P2.c) the following warning message appears:
# "Removed 45 rows containing non-finite values (stat_density)."

# Why did we get this warning message?
# Answer:there were 45 rows with NA values

# What could be done to avoid this message?
# Answer:ignore NA use na.rm inside the geom_density function
#
clean <- ggplot(flights_2, aes(x= air_time, color= carrier)) +
  geom_density(na.rm = TRUE) + ggtitle('Estimated Densities for each of Carrier') +
  theme(legend.position= "bottom") +
  ylab('Density') + # Proper axis label
  xlab('Air Time') # Proper axis label

clean #yeay, no more warning message

#------------------------------------------------
# P2.e)
#------------------------------------------------

# Using the magrittr format, define a function called 'speed' that takes a flights
# data.frame and adds a new column with value equal to the average speed in miles
# per hour.
# Plot bloxplots for the speed by month, for all flights from 'EWR' to 'LAX'.
# Save the plot as "plot_3.pdf".

# Solution
# Now the magrittr way:
speed <- . %>%  transform(ave = (distance/air_time)*60) #adds a new column with value equal to the average speed in miles per hour.

d<- flights_2 %>% speed

r<-ggplot(d, aes(x=as.factor(month), y=ave)) +
  ggtitle('The Speed for All Flights from EWR to LAX') +
  geom_boxplot(na.rm = TRUE)+
  ylim(350,550) +
  ylab('Average Speed (miles/h)') + # Proper axis label
  xlab('Month') # Proper axis label

r

ggsave(r, file = 'plot_3.pdf', width = 12, height = 8)

#------------------------------------------------
# P2.f)
#------------------------------------------------

# Create multiple scatterplots to visually explore how delay at departure affects
# delay at arrival by carriers ('EWR' to 'LAX' only).
# The scatterplots share the same y-axis but have different x-axes and different points
# colours.
# Save the plot as "plot_4.pdf".

# Solution

s <- ggplot(flights_2, aes(x = dep_delay , y = arr_delay)) +
  geom_point(aes(shape = carrier, colour = carrier), na.rm = TRUE) +
  facet_wrap(~carrier) +
  ggtitle('The Relationship between Departure Delay and Arrival Delay by Carriers (EWR to LAX Route)') +
  ylab('Arrival Delay') + # Proper axis label
  xlab('Departure Delay') + # Proper axis label
  theme(legend.position="none")  # Remove legend

s

ggsave(s, file = 'plot_4.pdf', width = 12, height = 8)


# End -------------------------------------------

