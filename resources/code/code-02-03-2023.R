data <- read.csv('data/squirtle.csv.gz', row.names = 1)

## pull out gene expression
gexp <- data[,4:ncol(data)]

## normalize
head(sort(apply(gexp, 2, var), decreasing=TRUE))
foo <- scale(gexp)
head(sort(apply(foo, 2, var), decreasing=TRUE))

totgenes <- rowSums(gexp)
good.cells <- names(totgenes)[totgenes > 0]
mat <- gexp[good.cells,]/totgenes[good.cells]
mat <- mat * 200
dim(mat)

## tSNE
library(Rtsne) ## install.packages('Rtsne')
?Rtsne
set.seed(0) ## for reproducibility
emb <- Rtsne(mat)
dim(emb$Y)

## ggplot
library(ggplot2)
df <- data.frame(emb$Y, gene1=mat[,'ERBB2'], gene2=mat[,'LUM'])
head(df)
p1 <- ggplot(df, aes(x = X1, y = X2, col=gene1)) + geom_point(size = 0.1) + 
  scale_color_gradient(low = 'lightgrey', high='red') + ggtitle('ERBB2') + theme_classic()

## try this out on your data
## patterns of genes in our 2D embedding space?
## try out different normalizations? or scaling?
## combining with linear dimensionality reduction?
p2 <- ggplot(df, aes(x = X1, y = X2, col=gene2)) + geom_point(size = 0.1) + 
  scale_color_gradient(low = 'lightgrey', high='red') + ggtitle('LUM') + theme_classic()
library(gridExtra)
grid.arrange(p1, p2, ncol=2)



