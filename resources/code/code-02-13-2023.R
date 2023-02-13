data <- read.csv('data/pikachu.csv.gz', row.names=1)
pos <- data[,1:2]
gexp <- data[, 4:ncol(data)]

good.cells <- rownames(gexp)[rowSums(gexp) > 10]
pos <- pos[good.cells,]
gexp <- gexp[good.cells,]

## food for thought: what do you think would happen
## if my gene panel was chosen to be mostly immune genes
## and my tissue has both immune cells and non-immune cells
totgexp <- rowSums(gexp)
mat <- gexp/totgexp
mat <- mat*median(totgexp)
mat <- log10(mat + 1)

set.seed(1)
emb <- Rtsne::Rtsne(mat)
com <- kmeans(mat, centers=5) ## you should have a better way

## plot my clusterings in embedding and tissue space
library(ggplot2)
df <- data.frame(pos, emb$Y, celltype=as.factor(com$cluster))
head(df)
p1 <- ggplot(df, aes(x = x_centroid, y = y_centroid, col=celltype)) + 
  geom_point(size = 0.1) + theme_classic()
p2 <- ggplot(df, aes(x = X1, y = X2, col=celltype)) + geom_point(size = 0.1) +
  theme_classic()
library(gridExtra)
grid.arrange(p1,p2)

## pick a cluster
cluster.of.interest <- names(which(com$cluster == 4))
cluster.other <- names(which(com$cluster != 4))
## loop through my genes and test each one
genes <- colnames(mat)
pvs <- sapply(genes, function(g) {
  a <- mat[cluster.of.interest, g]
  b <- mat[cluster.other, g]
  wilcox.test(a,b,alternative="two.sided")$p.val
})
names(which(pvs < 1e-8))
head(sort(pvs), n=20)

## visualize this gene
g <- 'KRT7'
df <- data.frame(pos, emb$Y, gene=mat[,g])
head(df)
p1 <- ggplot(df, aes(x = x_centroid, y = y_centroid, col=gene)) + 
  geom_point(size = 0.1) + theme_classic() + scale_color_continuous(low='lightgrey', high='red')
p2 <- ggplot(df, aes(x = X1, y = X2, col=gene)) + geom_point(size = 0.1) +
  theme_classic() + scale_color_continuous(low='lightgrey', high='red')
library(gridExtra)
grid.arrange(p1,p2)

### calculate a fold change
## for students who are more expert, see if you can do this in one loop with pvs
log2fc <- sapply(genes, function(g) {
  a <- mat[cluster.of.interest, g]
  b <- mat[cluster.other, g]
  log2(mean(a)/mean(b))
})
pvs['KRT7']
log2fc['KRT7']
pvs['CD3E']
log2fc['CD3E']

## volcano plot
df <- data.frame(pvs, log2fc)
head(df)
ggplot(df, aes(y=-log10(pvs), x=log2fc)) + geom_point()

library(ggrepel) ## install.packages("ggrepel")
ggplot(df, aes(y=-log10(pvs), x=log2fc)) + geom_point() +
  ggrepel::geom_label_repel(label=rownames(df))

## one hint for HW5: we are dealing with breast tissue
## bit.ly/GDV23_rc
