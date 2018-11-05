library("DMwR")
library("chemometrics")
library("VIM")
library("robustbase")
library("mice")
library("mvoutlier")
library("fBasics")
library("FactoMineR")
# 1
#Load the given dataset
dataset = read.csv("Russet_ineqdata.txt", sep = '\t')
#Fill the missing values of the dataset with KNN method
filedDataset = knnImputation(dataset, k = 5, scale = T)
#Define X matrix with continuous variables of the dataset
X = filedDataset[1:8]

# 2
#Own implementation of PCA function
#X matrix, Weights vector and Mode string are given as parameters
PCAFunction <- function(X, weights, distance = 'centered') {
  #Number of individuals
  rows = nrow(X)
  #Number of columns
  cols = ncol(X)
  
  # a
  #Define the matrix N of weights of individuals (with uniform weights)
  N = diag(1, rows, rows) * weights

  # b
  #Compute the centroid G of individuals
  centroidX = as.numeric(colMeans(X))
  
  # c
  #Compute the centered or standardized X matrix
  centeredX = (X - t(replicate(rows, centroidX)))
  colsd = colStdevs(X)
  standarizedX = (X - t(replicate(rows, centroidX))) / t(replicate(rows, colsd))
  
  X = standarizedX
  if(distance == 'centered') {
    X = centeredX
    print("centered")
  }
  
  # d
  #Compute the covariance or correlation matrix of X and diagonalize it
  covX = t(as.matrix(X)) %*% as.matrix(N) %*% as.matrix(X)
  diagCovX = eigen(covX)
  
  # e
  #Do the screeplot of the eigenvalues and define the number of significant dimension
  plot(diagCovX$values, main="Eigenvalues", type='b', cex.main=1, cex.axis=1, ylab = NA, xlab=NA, cex.axis=0.8)
  abline(h = 1, lty = 2)
  # Percentage of accumulated inertia
  percetangeOfVariance = (diagCovX$values / sum(diagCovX$values)) * 100
  accumulatedPercentageOfVariance = cumsum(percetangeOfVariance)
  text(diagCovX$values, labels = round(accumulatedPercentageOfVariance, digits = 2), col=2, cex=0.75)
  
  # f
  #Compute the projections of individuals in the significant dimensions
  individualsProjection = as.matrix(X) %*% diagCovX$vectors[,1:5]
  
  # g
  #Compute the projection of variables in the significant dimensions
  variablesProjection = cor(X, individualsProjection)
  if(distance == 'centered') {
    variablesProjection = (colsd * diag(1, cols, cols) %*% cor(X, individualsProjection))
    #variablesProjection = t(t(variablesProjection) * colsd)
  }
  View(variablesProjection)
  
  # h
  #Plot the individuals in the first factorial plane of Rp. 
  #Color the individuals according thedemovariable
  demo <- as.factor(filedDataset[,9])
  
  plot(individualsProjection[,1], individualsProjection[,2], pch=18, col=as.vector(demo), cex.main=1, cex.axis=1, ylab = NA, xlab=NA, cex.axis=0.8, main="Individuals proj.")
  text(individualsProjection[,1], individualsProjection[,2], labels=row.names(individualsProjection), col=as.vector(demo), cex=0.75)
  
  # i
  #Plot the variables (as arrows) in the first factorial plane of Rn
  if(distance == 'centered') {
    plot(variablesProjection[,1], variablesProjection[,2], type = "n", main="Variables proj.", xlab="Dim 1", ylab="Dim 2"
         , cex.main=1, cex.axis=1, cex.axis=0.8)
    theta <- seq(0, 2 * pi, length = 200)
    arrows(0, 0, variablesProjection[,1], variablesProjection[,2])
    text(variablesProjection[,1], variablesProjection[,2], labels=names(X), cex=0.75)
    
  } else{
    radius <- max(abs(variablesProjection[,1:2]), 1)
    plot(c(-radius, radius), c(-radius, radius), type = "n", main="Var. proj. (normalized)", xlab="Dim 1", ylab="Dim 2"
         , cex.main=1, cex.axis=1, cex.axis=0.8)
    theta <- seq(0, 2 * pi, length = 200)
    arrows(0, 0, variablesProjection[,1], variablesProjection[,2])
    text(variablesProjection[,1], variablesProjection[,2], labels=names(X), cex=0.75)
    lines(x = radius * cos(theta), y = radius * sin(theta))
    abline(h= 0, v=0, lty=2)
  }
}

#Execution of 2
weights = rep_len(1/nrow(X), nrow(X))
res.pca1 = PCAFunction(X, weights, 'normalized')
PCAFunction(X, weights, 'centered')

# 3 / 4
weights = rep_len(1/(nrow(X)-1), nrow(X))
weights[11] = 0
res.pca2 = PCAFunction(X, weights, 'normalized')
PCAFunction(X, weights, 'centered')
cor(res.pca2, res.pca1)

# 5
#PCA considering cuba as an Outlier
res.pca = PCA(filedDataset, quali.sup=c(9), ind.sup=c(11), scale.unit=TRUE, ncp=5, graph=T)
#PCA of all data
#res.pca1 = PCA(filedDataset, quali.sup=c(9), scale.unit=TRUE, ncp=5, graph=T)

# 6
#Best country in first factorial plane
c(which.max(res.pca$ind$cos2[,1] +res.pca$ind$cos2[,2]), res.pca$ind$cos2[which.max(res.pca$ind$cos2[,1] + res.pca$ind$cos2[,2]),])
#Worst country in first factorial plane
c(which.min(res.pca$ind$cos2[,1] +res.pca$ind$cos2[,2]), res.pca$ind$cos2[which.min(res.pca$ind$cos2[,1] + res.pca$ind$cos2[,2]),])

# 7
#Three countries that most influenciate the formation of the first principal component
sortedCountries = sort(res.pca$ind$contrib[,1], decreasing = TRUE)
sortedCountries[1:3]
#Three countries that most influenciate the formation of the second principal component
sortedCountries = sort(res.pca$ind$contrib[,2], decreasing = TRUE)
sortedCountries[1:3]

# 8
#Best variable in first factorial plane
c(which.max(res.pca$var$cos2[,1] + res.pca$var$cos2[,2]), res.pca$var$cos2[which.max(res.pca$var$cos2[,1] + res.pca$var$cos2[,2]), 1:2])
#Worst variable in first factorial plane
c(which.min(res.pca$var$cos2[,1] + res.pca$var$cos2[,2]), res.pca$var$cos2[which.min(res.pca$var$cos2[,1] + res.pca$var$cos2[,2]), 1:2])

# 9
#Three variables that most influenciate the formation of the first principal component
sortedVariables = sort(res.pca$var$contrib[,1], decreasing = TRUE)
sortedVariables[1:3]
#Three variables that most influenciate the formation of the second principal component
sortedVariables = sort(res.pca$var$contrib[,2], decreasing = TRUE)
sortedVariables[1:3]

# 10
#Modalities of significance of variable demo in the first two principal components
res.pca$quali.sup$v.test[,1:2]
