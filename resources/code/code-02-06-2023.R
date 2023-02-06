## review on loops
?lapply
out <- sapply(1:10, function(x) {
  return(x)
})
## for i in 1 through 10, perform function x

?apply
data <- read.csv('data/pikachu.csv.gz', row.names=1)
head(data)
gexp <- as.matrix(data[, 4:ncol(data)])
dim(gexp)

## for every row, compute median
med1 <- apply(gexp, 1, median) 
med2 <- sapply(1:nrow(gexp), function(x) {
  median(gexp[x,])
})

## kmeans
mat <- log10(gexp+1)
?kmeans
com <- kmeans(mat, centers = 10)
clusters <- as.factor(com$cluster) ## for categorical
com$tot.withinss
com$betweenss

table(clusters)
df <- data.frame(gene1 = mat[, 'MS4A1'], 
                 gene2 = mat[, 'LUM'], 
                 clusters = clusters)
library(ggplot2)
ggplot(df, aes(x = gene1, y = gene2, col=clusters)) +
  geom_point() 

## how can we find the optimal k for your dataset
## how should we visualize?
## applying kmeans to gene expression, pcs, 2D tSNE?

mat <- log10(gexp+1)
emb <- Rtsne::Rtsne(mat, perplexity = 100, check = FALSE)
set.seed(0)
com1 <- kmeans(emb$Y, centers=10)
set.seed(10) ## run twice
com2 <- kmeans(emb$Y, centers=10)
df <- data.frame(emb$Y, com1=as.factor(com1$cluster), com2=as.factor(com2$cluster))
p1 <- ggplot(df, aes(x = X1, y = X2, col=com1)) + geom_point(size = 0.1) + 
  theme_classic() + ggtitle('kmeans on tSNE embedding')
p2 <- ggplot(df, aes(x = X1, y = X2, col=com2)) + geom_point(size = 0.1) + 
  theme_classic() + ggtitle('kmeans on tSNE embedding')
gridExtra::grid.arrange(p1, p2)

## make sure you all can run loops
ks <- 1:20
out <- do.call(rbind, lapply(ks, function(k){
  com <- kmeans(emb$Y, centers=k)
  c(within = com$tot.withinss, between = com$betweenss)
}))
## check out stackoverflow of do.call, rbind, cbind, lapply, sapply
out
plot(ks, out[,1], type="l")
plot(ks, out[,2], type="l")





