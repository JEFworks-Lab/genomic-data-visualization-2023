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
pcs <- prcomp(mat)

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

## what if I want to use animation
library(gganimate) ## install.packages("gganimate")
## make a data frame that corresponds to cells in space
a <- data.frame(scale(pcs$x[,1:2]), celltype=as.factor(com$cluster), type="pcs 1 2")
b <- data.frame(scale(pcs$x[,3:4]), celltype=as.factor(com$cluster), type="pcs 3 4")
colnames(a) <- colnames(b) <- c('x', 'y', 'celltype', 'type')
df <- rbind(a,b)
p <- ggplot(df, aes(x=x, y=y, col=celltype)) + 
  geom_point(size=0.1)
anim <- p + gganimate::transition_states(type)
anim + labs(title = '{closest_state}')

## if you are getting a folder of images and not a gif:
## install.packages('magick') ## https://cran.r-project.org/web/packages/magick
## other options: install.packages('gifski')
## may need to restart

## try it out for your data
## animate between different normalizations?
## different PC vs tSNE?
## maybe you want visualize different genes

## bit.ly/GDV23_rc

