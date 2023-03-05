### STAT 341 A2 ###

# Problem 1

# 1b)

rho <- function(theta) {
  (2 * theta[1]^2) - (1.05 * theta[1]^4) + (theta[1]^6 /6) + 
    (theta[1] * theta[2]) + theta[2]^2
}

g <- function(theta) {
  c((4 * theta[1]) - (4.2 * theta[1]^3) + theta[1]^5 + theta[2], 
    theta[1] + (2 * theta[2]))
}

# 1c)

gradientDescent <- function(theta = 0, rhoFn, gradientFn, lineSearchFn, testConvergenceFn, 
                            maxIterations = 100, tolerance = 1e-06, relative = FALSE, lambdaStepsize = 0.01, 
                            lambdaMax = 0.5) {
  
  converged <- FALSE
  i <- 0
  
  while (!converged & i <= maxIterations) {
    g <- gradientFn(theta)  ## gradient
    glength <- sqrt(sum(g^2))  ## gradient direction
    if (glength > 0) 
      g <- g/glength
    
    lambda <- lineSearchFn(theta, rhoFn, g, lambdaStepsize = lambdaStepsize, 
                           lambdaMax = lambdaMax)
    
    thetaNew <- theta - lambda * g
    converged <- testConvergenceFn(thetaNew, theta, tolerance = tolerance, 
                                   relative = relative)
    theta <- thetaNew
    i <- i + 1
  }
  
  ## Return last value and whether converged or not
  list(theta = theta, converged = converged, iteration = i, fnValue = rhoFn(theta))
}

### line searching could be done as a simple grid search
gridLineSearch <- function(theta, rhoFn, g, lambdaStepsize = 0.01, lambdaMax = 1) { ## grid of lambda values to search
  lambdas <- seq(from = 0, by = lambdaStepsize, to = lambdaMax)
  ## line search
  rhoVals <- sapply(lambdas, function(lambda) {
    rhoFn(theta - lambda * g) 
  })
  ## Return the lambda that gave the minimum
  lambdas[which.min(rhoVals)] 
}

testConvergence <- function(thetaNew, thetaOld, tolerance = 1e-10, 
                            relative = FALSE) { 
  sum(abs(thetaNew - thetaOld)) < if (relative) 
    tolerance * sum(abs(thetaOld)) else tolerance 
}

# (−1.5, 1.5)
# converged to minima point A
gradientDescent(theta = c(-1.5, 1.5), rhoFn = rho, gradientFn = g, 
                lineSearchFn = gridLineSearch, testConvergenceFn = testConvergence, 
                maxIterations = 1000)

# (1.5, 1.5)
# converged to minima point C
gradientDescent(theta = c(1.5, 1.5), rhoFn = rho, gradientFn = g, 
                lineSearchFn = gridLineSearch, testConvergenceFn = testConvergence, 
                maxIterations = 1000)

# (1.5, −1.5)
# converged to minima point C
gradientDescent(theta = c(1.5, -1.5), rhoFn = rho, gradientFn = g, 
                lineSearchFn = gridLineSearch, testConvergenceFn = testConvergence, 
                maxIterations = 1000)

# (−1.5, −1.5)
# converged to minima point A
gradientDescent(theta = c(-1.5, -1.5), rhoFn = rho, gradientFn = g, 
                lineSearchFn = gridLineSearch, testConvergenceFn = testConvergence, 
                maxIterations = 1000)

# (−1,−1)
# converged to minima point B
gradientDescent(theta = c(-1, -1), rhoFn = rho, gradientFn = g, 
                lineSearchFn = gridLineSearch, testConvergenceFn = testConvergence, 
                maxIterations = 1000)

# 1d)

# values for theta1 and theta2 
xaxs <- yaxs <- seq(-2, 2, by=0.01)
# matrix to contain values of p(theta1,ptheta2) 
z <- matrix(NA, nrow=length(xaxs), ncol=length(yaxs))
for (i in 1:length(xaxs)) {
  for (j in 1:length(yaxs))
    z[i,j]=rho(c(xaxs[i],yaxs[j]))
}

# local minimas of three-hump camel function
xpnt <- c(-1.75,0,1.75)
ypnt <- c(0.875,0,-0.875)

# points from part c)
sxpnt <- c(-1.5,1.5,1.5,-1.5,-1)
sypnt <- c(1.5,1.5,-1.5,-1.5,-1)

# contour plot with local minimas and line segements
# connected with points from part c)
conlevels <- seq(0.5,10,by=0.5)
contour(xaxs, yaxs, z, levels=conlevels,
        xlab="Theta_1",
        ylab="Theta_2",
        main="Three-Hump Camel Function (2D)",
        cex.main=0.8,cex.lab=0.8)
points(xpnt,ypnt,pch=16,cex=0.6)
points(sxpnt,sypnt,pch=15,col="green",cex=0.6)
segments(-1.75,0.875,-1.5,1.5,col="green")
segments(1.75,-0.875,1.5,1.5,col="green")
segments(1.75,-0.875,1.5,-1.5,col="green")
segments(-1.75,0.875,-1.5,-1.5,col="green")
segments(0,0,-1,-1,col="green")
text(xpnt,ypnt,labels=c("A","B","C"),pos=c(1,1,1),cex=0.65)


### QUESTION 2 ###

gme <- read.csv("GME.csv")

# a)
adjc <- gme["Adj_Close"]
par(mfrow = c(1, 3))
# histogram
hist(adjc[1:127,],
     xlab = "Price in USD",
     main = "GME Closing prices \n Jul 28 2020 to Jan 27 2021",
     cex.main=0.9,breaks="Scott")
# boxplot
boxplot(adjc[1:127,], 
        main = "GME Closing prices \n Jul 28 2020 to Jan 27 2021",
        cex.main=0.9)
# quantile plot
qvals <- sort(adjc[1:127,])
pvals <- ppoints(127)
plot(pvals,qvals,pch=19,col=adjustcolor("grey", alpha = 0.5),
     xlim =c(0,1),
     xlab = "Proportion p",
     ylab = bquote("Q"["y"]~"(p)"),
     main = "GME Closing prices \n Jul 28 2020 to Jan 27 2021",
     cex.main=0.9)


# c)
adjc <- gme["Adj_Close"]
day <- gme["Day_Num"]
plot(day[1:127,],adjc[1:127,],
     xlab = "Day",
     ylab = "Closing Price in USD",
     main = "Closing price vs Day", pch=16)
abline(lm(adjc[1:127,] ~ day[1:127, ]), col = "red")
legend("topleft",legend="LS line",col="red",
       cex=0.75,bty = "n", lty = 1)

# d)

# get 3 outliers (index 125, 126 and 127)
indx <- which((gme$Adj_Close) > 70)
indx

gmed <- read.csv("GME.csv", header=T,
                 colClasses = c("character", "numeric", "numeric"))
model.pop <- lm(adjc[1:127,] ~ day[1:127, ])
theta.hat <- model.pop$coef

N = nrow(adjc)
delta = matrix(0, nrow = N, ncol = 2)

for (i in 1:N) { 
  temp.model = lm(Adj_Close ~ Day_Num,data=gme[-i, ])
  delta[i, ] = abs(theta.hat - temp.model$coef)
}

#par(mfrow = c(1,3))
plot(delta[, 1],
     ylab = bquote(Delta[alpha]), main = bquote("Influence on" ~ alpha), 
     pch = 19, col = adjustcolor("grey", 0.6))

obs = c(125,126,127)
text(obs, delta[obs, 1] + 0.001, obs, pos=2)

plot(delta[, 2], 
     ylab = bquote(Delta[beta]), main = bquote("Influence on" ~ beta), 
     pch = 19, col = adjustcolor("grey", 0.6))
text(obs + 3, delta[obs, 2], obs, pos=2)

delta2 = apply(X = delta, MARGIN = 1, FUN = function(z) {
  sqrt(sum(z^2)) })

plot(delta2, ylab = bquote(Delta), main = bquote("Influence on" ~ theta), 
     pch = 19, col = adjustcolor("grey", 0.6))

text(obs + 3, delta2[obs], obs, pos=2)

# e)
adjc <- gme["Adj_Close"]
day <- gme["Day_Num"]
adjc2 <- adjc[1:124, ]
day2 <- day[1:124, ]

plot(day[1:127,], adjc[1:127,],
     xlab = "Day",
     ylab = "Closing Price in USD",
     main = "Closing price vs Day", pch=16)
abline(lm(adjc[1:127,] ~ day[1:127, ]), col = "red")
abline(lm(adjc2 ~ day2), col = "green")
legend("topleft",legend=c("LS line","LS line without 3 most influencial observations"),
       col=c("red","green"), cex=0.60, bty = "n", lty = 1)


# Problem 2 part f

# ii)

# functions from guide

tukey.fn <- function(r, k) {
  val = (r^2)/2 - (r^4)/(2 * k^2) + (r^6)/(6 * k^4) 
  subr = abs(r) > k
  val[subr] = (k^2)/6
  return(val)
}

createRobustTukeyRho <- function(x, y, kval) { ## local variable
  ## Return this function
  function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    sum(tukey.fn(y - alpha - beta * x, k = kval))
  } 
}

tukey.fn.prime <- function(r, k) {
  val = r - (2 * r^3)/(k^2) + (r^5)/(k^4) 
  subr = abs(r) > k
  val[subr] = 0
  return(val)
}

createRobustTukeyGradient <- function(x, y, kval) { ## local variables
  function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    ru = y - alpha - beta * x
    rhok = tukey.fn.prime(ru, k = kval)
    -1 * c(sum(rhok * 1), sum(rhok * x))
  } 
}


# iii)

library(MASS)
adjc <- gme["Adj_Close"]
day <- gme["Day_Num"]

nlminb(c(0,1), objective=createRobustTukeyRho(day,adjc,4.685), 
       gradient=createRobustTukeyGradient(day,adjc,4.685))

# guide to confirm our answer is similar
rlm(Adj_Close ~ Day_Num, data = gme, psi = "psi.bisquare")


# iv)

adjc <- gme["Adj_Close"]
day <- gme["Day_Num"]
adjc2 <- adjc[1:124, ]
day2 <- day[1:124, ]

# create Tukey Bisquare Regression Line 
tukeyd <- 0:127
tukeycp <- 3.3106361 * tukeyd + 0.1307125

plot(day[1:127,], adjc[1:127,],
     xlab = "Day",
     ylab = "Closing Price in USD",
     main = "Closing price vs Day", pch=16)
abline(lm(adjc[1:127,] ~ day[1:127, ]), col = "red")
abline(lm(adjc2 ~ day2), col = "green")
abline(tukeyd,tukeycp, col = "blue")
legend("topleft",
       legend=c("LS line",
                "LS line without 3 most influencial observations",
                "Tukey Bisquare RegressionLine"),
       col=c("red","green", "Blue"), cex=0.60, bty = "n", lty = 1)


