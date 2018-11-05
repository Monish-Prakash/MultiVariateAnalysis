#Load needed packages
library("DMwR")
library("chemometrics")
library("VIM")
library("mice")
library('robustbase')

#The given dataset is loaded
dataset = read.csv("Russet_ineqdata.txt", sep = '\t')
#The total number of missing values is calculated
sum(is.na(dataset))
#The number of missing values per attribute is calculated
sapply(dataset, function(x) sum(is.na(x)))
#The average number of total values against the total dataset is calculated
mean(is.na(dataset))

#Another way to retrive the number of missing values per attribute
summary(dataset)

#Function that calculates the percentage of missing data
calculateMissingDataPerc <- function(x){
  sum(is.na(x))/length(x)*100
}
#The total percentage of missing data is calculated for each attribute of the dataset
apply(dataset,2,calculateMissingDataPerc)
#The total percentage of missing data is calculated for each instance of the dataset
apply(dataset,1,calculateMissingDataPerc)

#The pattern of missing data is retrieved (non-visual way)
md.pattern(dataset)

#The pattern of missing data is retrieved (plot / visual way)
aggr_plot <- aggr(dataset, col=c('navyblue','red'), numbers=TRUE, sortVars=FALSE, labels=names(dataset), 
                  cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

par(mfrow=c(2,3))
#Listwise deletion
noNAvaluesDataset <- na.omit(dataset)
position <- 1:nrow(noNAvaluesDataset)
plot(position, noNAvaluesDataset$Rent,xlab="",ylab="Rent",pch=3, col="black", main="Listwise deletion")

#Unconditional mean imputation
noNAvaluesDataset <- dataset
averageRent <- mean(na.omit(noNAvaluesDataset$Rent))
averageEcks <- mean(na.omit(noNAvaluesDataset$ecks))
noNAvaluesDataset$Rent[is.na(noNAvaluesDataset$Rent)] <- averageRent
noNAvaluesDataset$ecks[is.na(noNAvaluesDataset$ecks)] <- averageEcks
position <- 1:nrow(noNAvaluesDataset)
plot(position, noNAvaluesDataset$Rent,xlab="",ylab="Rent",pch=3, col=ifelse(position==2|position==30|position==35, "red", "black"), main="Unconditional mean imputation")


#Regression imputation
noNAvaluesDataset <- mice(dataset, method = "norm.predict", m = 1)
noNAvaluesDataset <- complete(noNAvaluesDataset)
position <- 1:nrow(noNAvaluesDataset)
plot(position, noNAvaluesDataset$Rent,xlab="",ylab="Rent",pch=3, col=ifelse(position==2|position==30|position==35, "red", "black"), main="Regression imputation")

#Stochastic imputation
noNAvaluesDataset <- mice(dataset,method="norm.nob",m=1,maxit=1,seed=1)
noNAvaluesDataset <- complete(noNAvaluesDataset)
position <- 1:nrow(noNAvaluesDataset)
plot(position, noNAvaluesDataset$Rent,xlab="",ylab="Rent",pch=3, col=ifelse(position==2|position==30|position==35, "red", "black"), main="Stochastic imputation")

#Knn imputation
noNAvaluesDataset = knnImputation(dataset, k = 1, scale = T)
position <- 1:nrow(noNAvaluesDataset)
plot(position, noNAvaluesDataset$Rent,xlab="",ylab="Rent",pch=3, col=ifelse(position==2|position==30|position==35, "red", "black"), main="Knn imputation")

#Chained equations
noNAvaluesDataset <- mice(dataset, m = 1) 
noNAvaluesDataset <- complete(noNAvaluesDataset)
position <- 1:nrow(noNAvaluesDataset)
plot(position, noNAvaluesDataset$Rent,xlab="",ylab="Rent",pch=3, col=ifelse(position==2|position==30|position==35, "red", "black"), main="Chained equations")

#Random Forests (not implemented)

#Knn imputation (find best K value)
par(mfrow=c(2,3))
valuesK = c(1,2,4,5,7,10)
position <- 1:nrow(noNAvaluesDataset)
for (kv in valuesK){
noNAvaluesDataset = knnImputation(dataset, k = kv, scale = T)
plot(position, noNAvaluesDataset$Rent,xlab="",ylab="Rent",pch=3, col=ifelse(position==2|position==30|position==35, "red", "black"), main=paste("Knn imputation ( K=",kv,")"))
}

#Function developed to compute mahalanobis Distance
mahalanobisDistance <- function(x, center, cov) {
  x <- sweep(x, 2L, center)
  cov <- solve(cov)
  setNames(rowSums(x %*% cov * x), rownames(x))
}

#Function to deal with outliers
multivariateOutliers <- function (X, quantile = 0.975, plot = TRUE) 
{
  X.mcd <- covMcd(X)
  md = sqrt(mahalanobisDistance(as.matrix(X), apply(X, 2, mean), cov(as.matrix(X))))
  rd = sqrt(mahalanobisDistance(as.matrix(X), X.mcd$center, X.mcd$cov))
  cutoff <- sqrt(qchisq(quantile, ncol(X)))
  if (plot) {
    par(mfrow = c(1, 2))
    plot(1:length(md), md, xlab = "Index of object", ylim = c(0, 
                                                              max(md, rd)), ylab = "Classical Mahalanobis distance")
    abline(h = sqrt(qchisq(quantile, ncol(X))), lty = 2)
    plot(1:length(rd), rd, xlab = "Index of object", ylim = c(0, 
                                                              max(md, rd)), ylab = "Robust Mahalanobis distance")
    abline(h = cutoff, lty = 2)
  }
  list(md = md, rd = rd, cutoff = cutoff)
}

adjustedQuantileOutliers <- function (x, delta = qchisq(0.975, df = ncol(x)), quan = 1/2, alpha = 0.05) 
{
  covr <- covMcd(x, alpha = quan)
  dist <- mahalanobisDistance(x, center = covr$center, cov = covr$cov)
  s <- sort(dist, index = TRUE)
  z <- x
  par(mfrow = c(1, 2))
  
  # Plot of cumulative probability of individuals regarding mahalanobis distance
  plot(s$x, (1:length(dist))/length(dist), col = "green", xlab = "Ordered squared robust distance", ylab = "Cumulative probability", type = "n")
  # Introduce index as point
  text(s$x, (1:length(dist))/length(dist), as.character(s$ix), col = "green", cex = 0.8)
  
  # Introduce pink line according chi-square cumualtive probability (0.01 steps)
  t <- seq(0, max(dist), by = 0.01)
  lines(t, pchisq(t, df = ncol(x)), col = "pink")
  
  # Showing vertical line regarding the mahalanobis distance according the introduced quantile (e.g 97.5%)
  abline(v = delta, col = "blue")
  text(x = delta, y = 0.4, paste(100 * (pchisq(delta, df = ncol(x))), "% Quantile", sep = ""), col = "blue", pos = 2, srt = 90, cex = 0.8)
  
  # Showing vertical line regarding the mahalanobis distance according the adjusted quantile (e.g 97.5%, calculated by the 
  # Adaptative rightweighted estimator) for the outlier threshold
  xarw <- arw(x, covr$center, covr$cov, alpha = alpha)
  if (xarw$cn < Inf) {
    abline(v = xarw$cn, col = "cadetblue1")
    text(x = xarw$cn, y = 0.4, "Adjusted threshold", col = "cadetblue1", pos = 4, srt = 90, cex = 0.8)
  }
  
  if (ncol(x) > 2) {
    p <- princomp(x, covmat = covr)
    # Take the value of the score (distance) for the first 2 principal components.
    z <- p$scores[, 1:2]
  }
  
  plot(z, col = "green", type = "n", main = "Outliers (adjusted threshold)")
  if (xarw$cn < Inf) {
    # Showing individuals treated as outliers in the two first principal components. (Disance > adaptive outlier threshold)
    text(z[dist > xarw$cn, 1], z[dist > xarw$cn, 2], dimnames(as.data.frame(x)[dist > xarw$cn, ])[[1]], col = "red", cex = 0.8)
  }
  # Showing individuals in the two first principal components. (Disance <= adaptive outlier threshold)
  text(z[dist <= xarw$cn, 1], z[dist <= xarw$cn, 2], dimnames(as.data.frame(x)[dist <= xarw$cn, ])[[1]], col = "green", cex = 0.8)
  # List of all individuals. The outliers are computed by comparing the mahalanobis distance with the maximum between the adjusted outlier threshold 
  # and the chi square quantile function with 9 degrees of freedom (9 dimensions of the dataset).
  o <- (dist > max(xarw$cn, qchisq(0.975, dim(x)[2])))
  subset(o, o == TRUE)
}

