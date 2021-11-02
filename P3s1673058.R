## Conor McLachlan Sayer S1673058
library(ggplot2)
library(debug)

mtrace(predict.linmod)
mtrace.off()

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
  
  ## Estimating fitted values mu
  mu <- qr.fitted(qrx, y)
  
  ## Estimating variance and sigma (sd) of the response & residuals
  RSS <- drop(t(y) %*% y - mu %*% y) # RSS = diff between y and fitted y
  var_y <- RSS/(n-p)
  sigma <- sqrt(var_y)
  
  ## Estimating covariance matrix of least squares estimators V
  
  rn <- nrow(qr.R(qrx))
  I <- `diag<-`(matrix(0,rn,rn),1) # Creating identity matrix, same size as R matrix
  inverseRT <- forwardsolve(t(qr.R(qrx)),I) # Find inverse of transpose of R matrix
  
  V <- backsolve(qr.R(qrx), inverseRT*var_y) # solve R*Var(Bhat) = R(-T) Sigma^2
  colnames(V) <- names(beta) # Add names corresponding to beta coefficients
  rownames(V) <- names(beta)
  
  ## flev: list of vector levels for each factor variable in dat
  
  factors <- which(sapply(dat, is.factor) == TRUE) # find factor variables in dat
  factors <- factors[names(factors) %in% all.vars(formula)] # only keep factors in formula
  
  if(length(factors) == 1){ # no need to use lapply if only one factor
    flev <- list(levels(dat[,factors]))
    names(flev) <- names(factors)
  }else if(length(factors) > 1){ # if > 1 factor, extract levels for every factor
    flev <- lapply(dat[,factors], levels)
    names(flev) <- names(factors)
  }else{ # else return empty list
    flev=list()
  }
  
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
  
  if(!x$yname %in% colnames(newdata)){ # if orig response variable not in newdata, add dummy data 
  newdata <- cbind(rep(0,nrow(newdata)), newdata) # adding dummy response data to newdata
  colnames(newdata) <- c(x$yname, colnames(newdata[-1]))
  }
  
  for (i in names(x$flev)){ #for every factor in newdata, ensure alignment with levels in original data  
    newdata[,i] <- factor(newdata[,i], levels = x$flev[[i]]) 
  }
  newdata <- model.matrix(x$formula, newdata) # Converting new data to model.matrix
  
  beta <- matrix(x$beta) # converting estimated betas to matrix
  
  prediction <- newdata %*% beta # as betas only px1, no need for QR of newdata 
  
  colnames(prediction) <- x$yname
  
  prediction
}

## Testing Functions 
#####

mpg1 <- data.frame(mpg) ## convert to regular DF
head(mpg1)
summary(mpg1)

mpg1$trans <- factor(gsub("\\(.*\\)","",mpg1$trans)) # strip out extra info on trans
mpg1$trans <- as.factor(mpg1$trans)
mpg1$class <- as.factor(mpg1$class)

## Normal test data

formula_mpg <- formula(cty ~ cyl + displ + year)

lm_test <- lm(formula_mpg, mpg1)
linmod_test <- linmod(formula_mpg, mpg1)

summary(lm_test)
linmod_test

plot.linmod(linmod_test)

newdata1 <- mpg1[sample(234,100),c("cyl", "displ", "year")]

sum(predict.linmod(linmod_test, newdata1) == predict.lm(lm_test, newdata1))

## Testing factors
formula_mpg2 <- formula(cty ~ cyl + trans + class)

lm_test2 <- lm(formula_mpg2, mpg1)
linmod_test2 <- linmod(formula_mpg2, mpg1)
summary(lm_test2)
linmod_test2

newdata2 <- mpg1[sample(234, 100),c("cyl","trans","class")]
newdata2$trans <- factor(newdata2$trans, levels=c("auto", "semi-auto", "manual"))
newdata2$class <- as.character(newdata2$class)

sum(predict.linmod(linmod_test2, newdata2)==predict.lm(lm_test2, newdata2))


#### IRIS data

formula_iris <- formula(Sepal.Length ~ Species * Sepal.Width + Petal.Length)

lm_test_iris <- lm(formula_iris, iris)
linmod_test_iris <- linmod(formula_iris, iris)

summary(lm_test_iris)
linmod_test_iris

newdata_iris <- iris[sample(150, 30),c("Species","Sepal.Width","Petal.Length")]
newdata_iris[,"Species"] <- factor(newdata_iris[,"Species"], levels=c("setosa", "versicolor", "virginica","fioih"))

sum(predict.linmod(linmod_test_iris, newdata_iris) == predict.lm(lm_test_iris, newdata_iris))

#### mtcars

formula_mtcars <- formula(mpg ~ disp + gear + wt)

lm_test_car <- lm(formula_mtcars, mtcars)
linmod_test_car <- linmod(formula_mtcars, mtcars)
summary(lm_test_car)
linmod_test_car

plot.linmod(linmod_test_car)

newdata_cars <- mtcars[sample(32,12),c("disp","gear","wt")]

sum(predict.linmod(linmod_test_car, newdata_cars) == predict.lm(lm_test_car, newdata_cars))


#### Toothgrowth
formula_toothgrowth <- formula(dose ~ supp * len)

lm_test_tooth <- lm(formula_toothgrowth, ToothGrowth)
linmod_test_tooth <- linmod(formula_toothgrowth, ToothGrowth)
summary(lm_test_tooth)
linmod_test_tooth


plot.linmod(linmod_test_tooth)

newdata_tooth <- ToothGrowth[sample(60, 15), c("supp","len")]
newdata_tooth$supp <- as.character(newdata_tooth$supp)

sum(predict.linmod(linmod_test_tooth, newdata_tooth)==predict.lm(lm_test_tooth, newdata_tooth))

#####