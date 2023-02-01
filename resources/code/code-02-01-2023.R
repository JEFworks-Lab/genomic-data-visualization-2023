############### dummy example
geneA <- c(1,2,3,11,12,13,7,7,7,7,7,7)
geneB <- rev(geneA)

gexp <- cbind(geneA, geneB)

?prcomp
pcs <- prcomp(gexp)
names(pcs)
pcs$sdev
plot(gexp)
plot(pcs$x)

plot(pcs$x[,1], gexp[,1])
plot(pcs$x[,2], gexp[,2])

pcs$rotation ## check loading

######################
data <- read.csv('data/bulbasaur.csv.gz', row.names = 1)
dim(data)
head(data)

## quality control
gexp <- data[, 4:ncol(data)]
head(gexp)
hist(log10(colSums(gexp)+1))
hist(log10(rowSums(gexp)+1))

## PCA
mat <- log10(gexp+1)
pcs <- prcomp(mat) ## raw counts
##### for those more advanced
## try out scale? 
## normalize to get proportion 
## normalize by area
## normalize by median
###### explore results
names(pcs)
dim(gexp)
length(pcs$sdev)
## scree plot
plot(1:100, pcs$sdev[1:100], type="l")
## look at loadings
head(sort(pcs$rotation[,1], decreasing=TRUE)) ## loading (coefficients on PCs)
head(sort(pcs$rotation[,1], decreasing=FALSE))

head(sort(pcs$rotation[,2], decreasing=TRUE))
head(sort(pcs$rotation[,2], decreasing=FALSE))
## if we plot ERBB2 vs LUM
ggplot(data, aes(x = ERBB2, y = KRT7)) + geom_point()

###### make a data visualization to explore our first two PCs
library(ggplot2)
df <- data.frame(pcs$x[,1:2], gene = log10(gexp[,'LUM']+1))
head(df)
ggplot(data = df, aes(x = PC1, y = PC2, col=gene)) + 
  geom_point() + 
  theme_classic() +
  scale_colour_gradient2(
    low = "red",
    mid = "white",
    high ="blue",
    midpoint = 1)

