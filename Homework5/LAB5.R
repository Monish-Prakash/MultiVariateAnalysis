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
dataset = read.csv("PCA_quetaltecaen.txt", sep = '\t')

# 2
K <- dataset[,-1]
rownames(K) <- dataset[, 1]

calculateF <- function(matrixIni){
  matrixF <- matrixIni/sum(matrixIni);
  return (matrixF)
}

F <- calculateF(K)
fi = rowSums(F)
fj = colSums(F)

Z = array(0, dim=c(nrow(F),ncol(F)))
Zstar = array(0, dim=c(nrow(F),ncol(F)))
X = array(0, dim=c(nrow(F),ncol(F)))
for (i in 1:nrow(F)){ 
  for (j in 1:ncol(F)) {
    Z[i, j] =  (F[i, j] - fi[i]*fj[j]) / (sqrt(fi[i] * fj[j]))
    Zstar[i, j] = F[i, j] / sqrt(fi[i] * fj[j])
    X[i, j] = (F[i, j] - fi[i]*fj[j]) / (fi[i] * sqrt(fj[j]))
  }
}

eig_Z = eigen(t(Z) %*% Z)
eig_Zstar = eigen(t(Zstar) %*% Zstar)

Psi_ia = array(0, dim=c(nrow(F),ncol(F)))
for (i in 1:nrow(F)) {
  for (a in 1:ncol(F)) {
    for (j in 1:ncol(F)) {
      Psi_ia[i, a] =  Psi_ia[i, a] + (F[i,j]/(fi[i]*sqrt(fj[j])))*eig_Z$vectors[j,a]
    }
  }
}

contrib_ia = array(0, dim=c(nrow(Psi_ia),ncol(Psi_ia)))
for (i in 1:nrow(F)) {
  for (a in 1:ncol(F)) {
    contrib_ia[i, a] = ((fi[i]*Psi_ia[i,a]^2) / eig_Z$values[a])*100
  }
}

Phi_ja = array(0, dim=c(ncol(F),nrow(F)))
for (j in 1:ncol(F)) {
  for (a in 1:ncol(F)) {
    for (i in 1:nrow(F)) {
      Phi_ja[j, a] =  Phi_ja[j, a] + (F[i,j]/fj[j])*Psi_ia[i,a]
    }
    Phi_ja[j, a] = Phi_ja[j, a] * 1/sqrt(eig_Z$values[a])
  }
}

contrib_ja = array(0, dim=c(ncol(Phi_ja),ncol(Phi_ja)))
for (j in 1:ncol(F)) {
  for (a in 1:ncol(F)) {
    contrib_ja[j, a] = ((fj[j]*Phi_ja[j,a]^2) / eig_Z$values[a])*100
  }
}
plot(cbind(Phi_ja[,1], Psi_ia[,1]), cbind(-Phi_ja[,2], -Psi_ia[,2]), main="CA factor map (developed)", xlab="Dim 1", ylab="Dim 2",  xlim=c(-0.3,0.4),  ylim=c(-0.15,0.15), pch=17)
abline(h = 0, v=0)

# 2b
ca.quetaltecaen <- CA(K)

plot(ca.quetaltecaen$eig[1:7], type='b', ylab = NA, xlab=NA)
title("Eigenvalues", line = 0.7)
# Percentage (+ accumulated)
text(ca.quetaltecaen$eig[1:7]+ 0.0003, labels = round(ca.quetaltecaen$eig[15:23], digits = 2), col=2)
# We pick 5 components due to 90% of inertia

# 3
T_I = 0
I = array(0, dim=c(nrow(F),ncol(F)))

for (i in 1:nrow(F)){ 
    for (j in 1:ncol(F)) {
        I[i, j] = ((F[i,j] - fi[i]*fj[j])/sqrt(fi[i]*fj[j]))^2
        T_I = T_I + P_I[i, j]
    }
}

T_I

for (i in 1:nrow(F)) {
    for (j in 1:ncol(F)) {
        P_I[i, j] = (I[i, j] / T_I)* 100
    }
}
P_I

P_I_Diag = sum(diag(P_I))
P_I_Diag


# 4
Z <- as.matrix(K);
for (j in 1:20){
  for (x in 1:nrow(Z)){
    n = sum(Z)
    f = calculateF(Z)
    fi = rowSums(f)
    fj = colSums(f)
    Z[x,x] <- n*fi[x]*fj[x]
    #View(Z)
  }
}

#We repeat exercise 3
T_I = 0
I = array(0, dim=c(nrow(f),ncol(f)))

for (i in 1:nrow(f)){ 
  for (j in 1:ncol(f)) {
    I[i, j] = ((f[i,j] - fi[i]*fj[j])^2/(fi[i]*fj[j]))
    T_I = T_I + I[i, j]
  }
}
#T_I

P_I = array(0, dim=c(nrow(I),ncol(I)))
for (i in 1:nrow(f)) {
  for (j in 1:ncol(f)) {
    P_I[i, j] = (I[i,j] / T_I)* 100
  }
}
P_I_Diag = sum(diag(P_I))
P_I_Diag

# 5
ca.quetaltecaen1 <- CA(Z)

plot(ca.quetaltecaen1$eig[1:7], type='b', ylab = NA, xlab=NA)
title("Eigenvalues", line = 0.7)
text(ca.quetaltecaen1$eig[1:7], labels = round(ca.quetaltecaen1$eig[15:23], digits = 2), col=2)