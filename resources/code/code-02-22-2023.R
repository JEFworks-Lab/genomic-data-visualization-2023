# bit.ly/GDV_codex
data <- read.csv('~/Downloads/codex_spleen_subset.csv.gz', row.names=1)
data[1:5,1:5]
dim(data)
## how many cells ## 11512
## how many proteins ## 28
## explore this data
## plot cells in the tissue
## see if you can find any weird artifacts in this data?

pexp <- data[,4:ncol(data)]
library(ggplot2)
ggplot(data, aes(x=x, y=y, col=log10(Podoplanin+1))) + 
  geom_point(size=0.01) +
  scale_color_gradient(low = 'lightgrey', high='red') + 
  theme_void() + theme(legend.position = "none") 

hist(log10(data$CD45RO+1))
hist(log10(data$Podoplanin+1))

## dimensionality reduction and kmeans clustering 
## help me figure out what is this thing we're looking at
plot(log10(data$area+1), log10(rowSums(pexp)+1), pch=".")
mat <- pexp/rowSums(pexp)*median(rowSums(pexp))
## variance of each protein
hist(apply(mat, 2, var), breaks=100)
sort(apply(mat, 2, var))
pcs <- prcomp(mat)
sort(abs(pcs$rotation[,1]))

df <- data.frame(pcs$x[,1:2], gene=data$CD15)
ggplot(df, aes(x=PC1, y=PC2, col=gene)) + 
  geom_point(size=0.01) +
  scale_color_gradient(low = 'lightgrey', high='red') + 
  theme_void() + theme(legend.position = "none") 
## what does this data visualization tell you about
## what additional normalizations or processing you may want to do?
## Xinyue: we should scale

emb <- Rtsne::Rtsne(mat)$Y
com <- kmeans(mat, centers=5)$cluster
df <- data.frame(x=data$x, y=data$y, emb, pcs$x, com=as.factor(com), area=log10(data$area+1))
#head(df)
ggplot(df, aes(x=X1, y=X2, col=com)) + 
  geom_point(size=0.01) 
ggplot(df, aes(x=PC1, y=PC2, col=com)) + 
  geom_point(size=0.01) 
ggplot(df, aes(x=x, y=y, col=com)) + 
  geom_point(size=0.01) 

