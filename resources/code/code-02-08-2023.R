data <- read.csv('data/bulbasaur.csv.gz', row.names=1)
data[1:5,1:5]
gexp <- data[, 4:ncol(data)]
mat <- log10(gexp + 1)
com <- kmeans(mat, centers=20)
emb <- Rtsne::Rtsne(mat)
df <- data.frame(emb$Y, 
                 clusters=as.factor(com$cluster),
                 gene=mat[, 'MS4A1'])
p1 <- ggplot(df, aes(x = X1, y = X2, col=clusters)) +
  geom_point(size=0.1)
p2 <- ggplot(df, aes(x = X1, y = X2, col=gene)) +
  geom_point(size=0.1) + 
  scale_color_continuous(low = 'lightgrey', high='red')
grid.arrange(p1, p2)

## i have clusters
## i want to characterize differentially upregulated genes
ggplot(df, aes(x = X1, y = X2, col=clusters == 2)) +
  geom_point(size=0.1)

## grab cells in cluster 2
cluster.cells <- which(com$cluster == 2)
other.cells <- which(com$cluster != 2)
par(mfrow=c(2,1))
hist(mat[cluster.cells, 'MS4A1'])
hist(mat[other.cells, 'MS4A1'])
## test for differential expression
?wilcox.test
wilcox.test(mat[cluster.cells, 'MS4A1'],
            mat[other.cells, 'MS4A1'],
            alternative = 'greater')

## write a loop to test all genes
out <- sapply(colnames(mat), function(g) {
  wilcox.test(mat[cluster.cells, g],
              mat[other.cells, g],
              alternative = 'greater')$p.value
})
diff.up.genes <- names(out)[out < 0.05]
sort(out[diff.up.genes]) ## look at list

## check out some of these significantly upregulated genes
test.gene <- 'KIT'
df <- data.frame(data[,1:2], emb$Y, 
                 clusters=as.factor(com$cluster),
                 gene=mat[, test.gene])
p1 <- ggplot(df, aes(x = X1, y = X2, col=clusters)) +
  geom_point(size=0.1)
p2 <- ggplot(df, aes(x = X1, y = X2, col=gene)) +
  geom_point(size=0.1) + 
  scale_color_continuous(low = 'lightgrey', high='red') + 
  ggtitle(test.gene)
grid.arrange(p1, p2)

## pick a cluster
## find genes upregulated in that cluster
## advanced: are all your clusters "real"?
## what happens if you pick a really really large k?

## we can also compare clusters to each other
cluster1.cells <- which(com$cluster == 1)
cluster2.cells <- which(com$cluster == 2)

wilcox.test(mat[cluster1.cells, 'MS4A1'],
            mat[cluster2.cells, 'MS4A1'],
            alternative = 'greater')
?t.test
t.test(mat[cluster1.cells, 'MS4A1'],
            mat[cluster2.cells, 'MS4A1'],
            alternative = 'greater')$p.value
