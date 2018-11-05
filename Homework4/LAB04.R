library("DMwR")
library("chemometrics")
library("VIM")
library("robustbase")
library("mice")
library("mvoutlier")
library("fBasics")
library("FactoMineR")
library("factoextra")
library(cluster)
library(HSAUR)
# 1
#Load the given dataset
dataset = read.csv("Russet_ineqdata.txt", sep = '\t')
#Fill the missing values of the dataset with KNN method
filedDataset = knnImputation(dataset, k = 5, scale = T)

# 2
res.pca = PCA(filedDataset, quali.sup=c(9), ind.sup=c(11), scale.unit=TRUE, ncp=5, graph=T)

# 3
plot(res.pca$eig[1:8]+0.2, main="Eigenvalues (normalized)", type='b', cex.main=1, cex.axis=1, ylab = NA, xlab=NA, cex.axis=0.8)
abline(h = mean(res.pca$eig[1:8]), lty = 2)
# Percentage (+ accumulated)
text(res.pca$eig[1:8], labels = round(res.pca$eig[17:25], digits = 2), col=2, cex=0.75)
# We pick 5 components due to 90% of inertia

# 4
d <- dist(rbind(res.pca$ind$coord), method = "euclidean")
hc <- hclust(d, method = "ward.D2")

plot(hc)


barplot(hc$height)
text(hc$height, labels=round(hc$height, digits=2), col=2, cex = 0.6)

nc = 5

fviz_dend(hc, k = nc, cex = 0.6,  k_colors = "jco", rect = TRUE, rect_border = "jco", rect_fill = TRUE)


c1 <- cutree(hc,nc)

# VISUAL PARTITION HCLUST
plot(res.pca$ind$coord[,1],res.pca$ind$coord[,2],type="n",main=paste("Clustering of expenses in", toString(nc),"classes (HC)"))
text(res.pca$ind$coord[,1],res.pca$ind$coord[,2],col=c1,labels=rownames(res.pca$ind$coord),cex = 0.6)
abline(h=0,v=0,col="gray")


cdg<- aggregate(as.data.frame(res.pca$ind$coord),list(c1),mean)[,2:(5+1)]

kmeans = kmeans(res.pca$ind$coord, centers=cdg)
Bss <- sum(rowSums(kmeans$centers^2)*kmeans$size)
Wss <- sum(kmeans$withinss)
Ib5 <- 100*Bss/(Bss+Wss)
Ib5

#fviz_cluster(object = kmeans, data = filedDataset[-c(11),], geom = "point", stand = FALSE, ellipse.type = "norm") + theme_bw()
#VISUAL PARTITION USING KMEANS
plot(res.pca$ind$coord[,1],res.pca$ind$coord[,2],type="n",main=paste("Clustering of expenses in", toString(nc),"classes (kmeans)"))
text(res.pca$ind$coord[,1],res.pca$ind$coord[,2],col=kmeans$cluster,labels=rownames(res.pca$ind$coord),cex = 0.6)
abline(h=0,v=0,col="gray")
# 5
catdes = catdes(cbind(as.factor(kmeans$cluster),filedDataset[-11,]),1)

plot.catdes(catdes)
print(catdes)

# 6
euc_dist <- function(x1, c1){ 
  total <- 0;
  for(i in 1:ncol(x1)) {
    diff = sum((x1[,i] - c1[i])^ 2);
    total = total + diff;
  }
  return(sqrt(total)) 
}

#result <- euc_dist(res.pca$ind.sup$coord, kmeans$centers[1,])
clusters <- c()
for(j in 1:nrow(kmeans$centers)){
  result <- euc_dist(res.pca$ind.sup$coord, kmeans$centers[j,])
  clusters <- cbind(clusters, result)
}
print(clusters[1,])

