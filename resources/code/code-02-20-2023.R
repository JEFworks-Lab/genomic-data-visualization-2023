# bit.ly/GDV23_visiumBC
data <- read.csv('~/Downloads/visium_breast_cancer.csv.gz', row.names=1)
## how many spots are there?
dim(data) ## 2518
## how many genes are there? 
pos <- data[,1:2]
gexp <- data[,3:ncol(data)]
dim(gexp)
## how many mitochondrial genes are there? 'MT-' 
hist(log10(colSums(gexp)+1))
table(grepl('^MT-', colnames(gexp)))
## is this human or mouse? ## human

## QC
hist(log10(rowSums(gexp)+1))

p1 <- ggplot(data, aes(x=V6, y=V5, col=HES4)) + geom_point(size=0.5) +
  scale_color_gradient(low = 'lightgrey', high='red') + theme_classic()
p2 <- ggplot(data, aes(x=V6, y=V5, col=CD68)) + geom_point(size=0.5) +
  scale_color_gradient(low = 'lightgrey', high='red') + theme_classic()
gridExtra::grid.arrange(p1, p2)

ggplot(data, aes(x=V6, y=V5, col=rowSums(gexp))) + geom_point(size=0.5) +
  scale_color_gradient(low = 'lightgrey', high='red') + theme_classic()

## normalize
mat <- gexp/rowSums(gexp)
mat <- mat*mean(rowSums(gexp))

## after normalization
p1 <- ggplot(data, aes(x=V6, y=V5, col=mat[,'HES4'])) + geom_point(size=0.5) +
  scale_color_gradient(low = 'lightgrey', high='red') + theme_classic()
p2 <- ggplot(data, aes(x=V6, y=V5, col=mat[,'CD68'])) + geom_point(size=0.5) +
  scale_color_gradient(low = 'lightgrey', high='red') + theme_classic()
gridExtra::grid.arrange(p1, p2)

## if I don't normalize, what do you think would happen?
## probably find a cluster that just corresponds 
## to spots with a lot of cells?

## PCA
pcs <- prcomp(mat) ## note no log transformation for now
## irlba, and other faster SVD methods, faster than prcomp
## plot variance explained
length(pcs$sdev)
plot(pcs$sdev[1:30], type="l")

ggplot(data.frame(pcs$x[,1:30]), aes(x=PC20, y=PC21, col=mat[,'CD68'])) + 
  geom_point(size=0.5) +
  scale_color_gradient(low = 'lightgrey', high='red') + 
  theme_classic()

## kmeans 
## should we cluster on the genes?
## on the PCs?
## how many PCs?
set.seed(0)
com1 <- kmeans(pcs$x[,1:30], centers=5)
com2 <- kmeans(mat, centers=5)
com3 <- kmeans(mat[,1:30], centers=5)

df <- data.frame(pos, 
                 com1=as.factor(com1$cluster), 
                 com2=as.factor(com2$cluster), 
                 com3=as.factor(com3$cluster))
p1 <- ggplot(df, aes(x=V6, y=V5, col=com1)) + geom_point(size=0.5)
p2 <- ggplot(df, aes(x=V6, y=V5, col=com2)) + geom_point(size=0.5)
p3 <- ggplot(df, aes(x=V6, y=V5, col=com3)) + geom_point(size=0.5)
gridExtra::grid.arrange(p1, p2)
gridExtra::grid.arrange(p1, p3)
gridExtra::grid.arrange(p2, p3)

heatmap(table(com1$cluster, com3$cluster))
