library("DMwR")
library("chemometrics")
library("VIM")
library("robustbase")
library("mice")
library("mvoutlier")
library("FactoMineR")
library("timeSeries")
library("bigpca")
library("calibrate")

# 1
#Load the given dataset
dataset <- read.delim("Russet_ineqdata.txt")
#Fill the missing values of the dataset with KNN method
filedDataset = knnImputation(dataset, k = 5, scale = T)
#Define X matrix with continuous variables of the dataset
X = filedDataset[1:8]
#Number of rows
rows = nrow(X)
#Number of columns
cols = ncol(X)
#Compute the centroid G of individuals
centroidX = as.numeric(colMeans(X))
colsd = colStdevs(X)
standarizedX = (X - t(replicate(rows, centroidX))) / t(replicate(rows, colsd))

#2
significantDimensions = 5
standarizedX <- standarizedX[-c(11),]
m0 <- nipals(standarizedX, significantDimensions, it=70)

#3
biplot(m0$T, m0$P, col=c('grey', 'red'), xlab='PCA1', ylab='PCA2', cex=c(0.7,0.7), cex.axis=0.7)

#4
pc.rot <- varimax(m0$P)

iden = row.names(standarizedX); etiq = names(standarizedX)
Phi.rot = pc.rot$loadings[1:8,]
lmb.rot = diag(t(pc.rot$loadings) %*% pc.rot$loadings)

ze = rep(0,8)
plot(Phi.rot,type="n",xlim=c(-1,1),ylim=c(-1,1))
text(Phi.rot,labels=etiq, col="blue")
arrows(ze, ze, Phi.rot[,1], Phi.rot[,2], length = 0.07,col="blue")
abline(h=0,v=0,col="gray")

#5
Psi_stan.rot = as.matrix(standarizedX) %*% solve(cor(as.matrix(standarizedX))) %*% Phi.rot
Psi.rot = Psi_stan.rot %*% diag(sqrt(lmb.rot))

plot(Psi.rot,type="n")
text(Psi.rot,labels=iden)
abline(h=0,v=0,col="gray")

pca = PCA(filedDataset[,1:8], ind.sup = 11, graph=FALSE)
pca$ind$coord[,1:5] = Psi.rot
dimdesc(pca, axes=1:5)

#6
dataset2 <- read.delim("PCA_quetaltecaen.txt")

#7
symmetricMatrix = dataset2[,2:9]
rownames(symmetricMatrix) = dataset2[, 1]
colnames(symmetricMatrix) = dataset2[, 1]
for(i in 1:nrow(symmetricMatrix)) {
  for(j in 1:i) {
    jointFeeling = (symmetricMatrix[i, j] + symmetricMatrix[j, i]) / 2
    symmetricMatrix[i, j] = jointFeeling;
    symmetricMatrix[j, i] = jointFeeling;
  }
}
View(symmetricMatrix)

#8
dissimilarityMatrix = matrix(0,nrow(symmetricMatrix), ncol(symmetricMatrix))
rownames(dissimilarityMatrix) = dataset2[, 1]
colnames(dissimilarityMatrix) = dataset2[, 1]
for(i in 1:nrow(symmetricMatrix)) {
  for(j in 1:i) {
    dissimilarityMatrix[i, j] = 10 - symmetricMatrix[i, j];
    dissimilarityMatrix[j, i] = 10 - symmetricMatrix[j, i];
  }
}

#9 && 10
mds = cmdscale(dissimilarityMatrix)
plot(mds[,1], mds[,2], pch = 19, cex.axis=0.7, xlim=c(-2.5,2.5))
text(mds[,1], mds[,2], labels = dataset2[,1], pos = 4, cex=c(0.7,0.7))
abline(h= 0, v=0, lty=2)

