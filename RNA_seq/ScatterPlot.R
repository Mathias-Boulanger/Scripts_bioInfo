#!/usr/bin/R

# Charger la matrice des RPKM dans la variable RPKM_mat
#RPKM_mat <- read.table("/home/chipseq/Bureau/Formation/Rnaseq/matrice_RPKM_hg19.txt", dec='.', header=TRUE,row.names=1)
RPKM_mat <- read.table("/home/chipseq/Bureau/Formation/Rnaseq/matrice_RPKM_ERCC92.txt", dec='.', header=TRUE,row.names=1)
# Afficher les premières lignes de RPKM_mat
head(RPKM_mat)
head(RPKM_mat[,7:12])
head(RPKM_mat[,1:6])

# Calculer la moyenne des mutants et des wildtype pour chaque  gène
RPKM_mat$MeanA <- rowMeans(RPKM_mat[,1:6])
RPKM_mat$MeanB <- rowMeans(RPKM_mat[,7:12])

# Ré-afficher les premières lignes de RPKM_mat
head(RPKM_mat)


# Applique un log2 sur la matrice+1
Log2_RPKM_mat = log2(RPKM_mat+1)


############################################## IMPORTANT: Variable à changer
#####  Variables à changer

### Variable 1 
Xlab="MeanA"
X=Log2_RPKM_mat$MeanA

### Variable 2 
Ylab="MeanB"
Y=Log2_RPKM_mat$MeanB

############################################## IMPORTANT: Variable à changer
### plot
postscript(paste("Plot_log2RPKM_",Xlab, "_",Ylab,".ps", sep= ""));

plot(X,Y, xlab=paste("log2 RPKM of  ",Xlab), ylab=paste("log2 RPKM of ",Ylab), pch=20,col = ifelse((((X/Y)>=1.5 | (Y/X)>=1.5)), 'RED', 'black'))
myline.fit <- lm(X ~ Y) 
abline(myline.fit, col="red" )
cor(X,Y)

dev.off()



#cor(Log2_RPKM_mat)

