## Conor McLachlan Sayer S1673058
library(ggplot2)
library(debug)

#mtrace(linmod)

##TODO Write unit tests for code (for example, converting inputs to correct types)

## Writing Functions
linmod <- function(formula, dat){

  ## Setting up data for calculating linear model 
  formula <- formula(formula) # Converting formula to class formula (incase not)
  yname <- all.vars(formula)[1] # y - obtaining name of response variable
  y <- model.frame(formula,dat)[[1]] # y - obtaining original response variables
  
  x <- model.matrix(formula, dat) # creating model matrix for ind variables
  
  n <- nrow(x)
  p <- ncol(x)
  
  ## Calculating linear model
  
  ## Estimating beta vector
  qrx <- qr(x) ## Get QR decomposition of model matrix
  Qty <- qr.qty(qrx, y)[1:p] # get Q(T) y
  beta <- backsolve(qr.R(qrx), Qty) # solve RB = Q(t)y for B
  names(beta) <- colnames(x)
  
  ## Estimating fitted valued mu
  mu <- qr.fitted(qrx, y)
  
  ## Estimating variance and sigma (sd) of the response & residuals
  RSS <- drop(t(y) %*% y - mu %*% y) # RSS = diff between y and fitted y
  var_y <- RSS/(n-p)
  sigma <- sqrt(var_y)
  
  ## Estimating covariance matrix of least squares estimators V
  
  rn <- nrow(qr.R(qrx))
  I <- `diag<-`(matrix(0,rn,rn),1)
  inverseRT <- forwardsolve(t(qr.R(qrx)),I) # First find R(-T)
  
  V <- backsolve(qr.R(qrx), inverseRT*var_y) # solve R*Var(Bhat) = R(-T) Sigma^2
  colnames(V) <- names(beta)
  rownames(V) <- names(beta)

  ## flev, containing a list of vector levels for each factor variable in dat
  flev <- 0
  ### TODO Make this work, so factors work with predict function 
  
  
  # Returning a list containing specified elements (of class "linmod")
  
  returnlist <- list(beta, V, mu, y, yname, formula, flev, sigma)
  
  names(returnlist) <- c("beta", "V", "mu", "y", "yname", 
                         "formula", "flev", "sigma")
  
  class(returnlist) <- "linmod"
  
  returnlist
}

print.linmod <- function(x,...){
  estimates <- cbind(x$beta, sqrt(diag(x$V)))
  colnames(estimates) <- c("Estimate","s.e.")
  print(x$formula)
  print(estimates)
}

plot.linmod <- function(x,...){
  res <- x$y - x$mu
  plot(x$mu, res, # plotting residuals against fitted values of y
       xlab="fitted values",
       ylab="residuals")
  abline(a=0, b=0, lty = 2)
}

predict.linmod <- function(x, newdata){

  newdata <- cbind(rep(0,nrow(newdata)), newdata) # adding dummy response data to newdata
  colnames(newdata) <- c(x$yname, colnames(newdata[-1]))
  
  ## TODO: Make work with factors
  newdata <- model.matrix(x$formula, newdata) # Converting new data to model.matrix
  
  beta <- matrix(x$beta) # converting estimated betas to model.matrix
  
  prediction <- newdata %*% beta
  
  colnames(prediction) <- x$yname
  
  prediction
}


## Testing Functions 

summary(lm(cty ~ cyl*hwy, mpg))
linmod(cty ~ cyl*hwy, mpg)


mpg1 <- data.frame(mpg) ## convert to regular DF
head(mpg1)
summary(mpg1)

mpg1$trans <- factor(gsub("\\(.*\\)","",mpg1$trans)) # strip out extra info on trans

formula <- formula(cty ~ trans * displ)

summary(lm(cty ~ trans * displ, mpg1))
test <- linmod(cty ~ trans * displ, mpg1)

plot.linmod(test)



newdata1 <- as.data.frame(cbind(matrix(rep(c(2,4,6), 20), 60, 1), matrix(floor(runif(n=60, 26, 40)))))
colnames(newdata1) <- c("cyl","hwy")


## Problem: factors which aren't of type 'factor' will break the predict function
newdata2 <- as.data.frame(cbind(matrix(rep(c("auto", "manual"), 30), 60, 1), matrix(round(runif(n=60, 1.6, 7),1))))
colnames(newdata2) <- c("trans","displ")
newdata2$trans <- as.factor(newdata2$trans)
newdata2$displ <- as.numeric(newdata2$trans)



predict.linmod(test, newdata2)



