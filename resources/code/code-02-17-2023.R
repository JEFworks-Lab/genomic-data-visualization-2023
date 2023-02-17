## download: https://bit.ly/GDV23_mystery
data <- read.csv('~/Downloads/mystery_data.csv.gz', row.names=1)

?head
head(data)
class(data)
dim(data) ## 2000 cells, 1000 genes
data[1:5,1:5] ## probably mouse tissue

## QC
pos <- data[,1:2]
gexp <- data[,3:ncol(data)]
dim(pos)
dim(gexp)
hist(as.matrix(gexp))
hist(log10(colSums(gexp)+1))
hist(log10(rowSums(gexp)+1))

## spatial position
library(ggplot2)
colnames(data)[grepl('Cox', colnames(data))]
ggplot(data, aes(x=x, y=y, col=rowSums(gexp))) + geom_point(size=0.5) +
  scale_color_gradient(low = 'lightgrey', high='red') + 
  theme_classic()
colnames(data)

hist(gexp[, 'Cox20'])
hist(log10(gexp[, 'Cox20']+1))

## normalize
mat <- gexp/rowSums(gexp) * 5000
range(mat)
mat <- log10(mat+1)
## visualize variance on genes
hist(apply(mat, 2, var))

plot(apply(mat, 2, mean), apply(mat, 2, var), pch=".")

mat.scale <- scale(mat) ## scale to make variance 1
apply(mat.scale, 2, var)

## pca
pcs <- prcomp(mat)
df <- data.frame(pcs$x[,1:10], totgexp = rowSums(gexp))
p1 <- ggplot(df, aes(x=PC1, y=PC2, col=totgexp)) + geom_point()
p2 <- ggplot(df, aes(x=PC3, y=PC4), col=totgexp) + geom_point()
gridExtra::grid.arrange(p1, p2)
par(mfrow=c(1,1))
plot(pcs$sdev[1:20])
head(sort(pcs$rotation[,1]))
tail(sort(pcs$rotation[,1])) 
hist(pcs$rotation[,1])
## genes with high absolute loadings
## related to the brain

## clustering
results <- do.call(rbind, lapply(seq_len(15), function(k) {
  out <- kmeans(mat, centers=k)
  c(out$tot.withinss, out$betweenss)
}))
par(mfrow=c(1,2))
plot(seq_len(15), results[,1], main='tot.withinss')
plot(seq_len(15), results[,2], main='betweenss')

com <- kmeans(mat, centers=2, iter.max = 50)
cluster <- com$cluster

df <- data.frame(cluster=as.factor(cluster), pos, pcs$x[,1:4])
head(df)
ggplot(df, aes(x=x, y=y, col=cluster)) + geom_point()
ggplot(df, aes(x=PC1, y=PC2, col=cluster)) + geom_point()
ggplot(df, aes(x=PC2, y=PC3, col=cluster)) + geom_point()

head(sort(pcs$rotation[,1]))
tail(sort(pcs$rotation[,1]))

## tsne
emb <- Rtsne::Rtsne(mat)
df <- data.frame(emb$Y, cluster=as.factor(cluster), pos, gene=gexp[, 'Kcnip4'])
head(df)
ggplot(df, aes(x=X1, y=X2, col=cluster==14)) + geom_point()
ggplot(df, aes(x=X1, y=X2, col=gene)) + geom_point() + 
  scale_color_gradient(low = 'lightgrey', high='red')
ggplot(df, aes(x=x, y=y, col=cluster==3)) + geom_point()

## kmeans on tsne space?
com2 <- kmeans(emb$Y, centers=15)
df <- data.frame(emb$Y, cluster=as.factor(com2$cluster), pos)
head(df)
ggplot(df, aes(x=x, y=y, col=cluster)) + geom_point()

## differentially expressed genes
diffgexp <- sapply(colnames(gexp), function(g) {
  ct3 <- names(which(cluster==14))
  ctother <- names(which(cluster!=14))
  wilcox.test(gexp[ct3, g], gexp[ctother, g])$p.value
})
table(diffgexp < 0.05)
diffgexp[names(which(diffgexp < 0.05))]

logfc <- sapply(colnames(gexp), function(g) {
  ct3 <- names(which(cluster==14))
  ctother <- names(which(cluster!=14))
  log2(mean(gexp[ct3, g])/mean(gexp[ctother, g]))
})
df <- data.frame(logfc, pv=-log10(diffgexp))
ggplot(df, aes(x=logfc, y=pv)) + geom_point() 

logfc[logfc > 5]

## bit.ly/GDV23_rc
## HW5 is due tonight


