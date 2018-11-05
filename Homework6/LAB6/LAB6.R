library("DMwR")
library("chemometrics")
library("VIM")
library("robustbase")
library("mice")
library("mvoutlier")
library("fBasics")
library("FactoMineR")
library("factoextra")

# 1
#Load the given dataset
dataset = read.csv("mca_car.csv")
#Check missing values --> No missing values found, otherwise they would be calculated
which(is.na(dataset))
summary(dataset)

# 2
#MCA is performed
mca.car <- MCA(dataset, quanti.sup=c(17), quali.sup=c(18), ncp=19, graph=FALSE)
plot(mca.car,choix="var")

# 3
#Factors equival dimensions?
fviz_mca_biplot(mca.car,labelsize = 3)


# 4
#Average eigen values

plot(mca.car$eig[,1], type='l', ylab = NA, xlab=NA)
title("Eigenvalues", line = 0.7)
# Percentage (+ accumulated)

meanEig <- mean(mca.car$eig[,1])
recompEign <- mca.car$eig[mca.car$eig[,1]>meanEig, 1]
recompEign <- recompEign - meanEig

# PAG 67
plot(recompEign/sum(recompEign), type="l")
cumsum(100*recompEign/sum(recompEign))

# 90% INERTIA -> 9 DIM
nd = 9

# 5
d <- dist(rbind(mca.car$ind$coord[,1:9]), method = "euclidean")
hc <- hclust(d, method = "ward.D2")
barplot(hc$height, main="Heights")
abline(h=9.5, col = 2)

plot(as.dendrogram(hc), labels=NULL)
abline(h=9.5, col = 2)
nd = 9

nc = 6
c1 <- cutree(hc,nc)

cdg <- aggregate(as.data.frame(mca.car$ind$coord[,1:nd]),list(c1),mean)[,2:(nd+1)]
Bss <- sum(rowSums(cdg^2)*as.numeric(table(c1)))
Tss <- sum(rowSums(mca.car$ind$coord[,1:9]^2))
Ib6 <- 100*Bss/Tss
Ib6

kmeans = kmeans(mca.car$ind$coord[,1:9], centers=cdg)

plot(mca.car$ind$coord[,1],mca.car$ind$coord[,2],main=paste("Clustering in", toString(nc),"classes (Kmeans)"),col=kmeans$cluster,pch=20,cex=1.0)
abline(h=0,v=0,col="gray")

# 6
catdesPlot = catdes(cbind(as.factor(kmeans$cluster),dataset[, -c(17,18)]),1)
#plot.catdes(catdesPlot)

#catdesPlot
catdesPlot$category$`1`
catdesPlot$category$`2`
catdesPlot$category$`3`
catdesPlot$category$`4`
catdesPlot$category$`5`
catdesPlot$category$`6`