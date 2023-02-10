#### Lesson 1 (1-25-23): general preliminary investigation of data ####

data <- read.csv('data/charmander.csv.gz')
dim(data) # row and columns
head(data, n=6) # return n rows of data with column names 
data[1:5,1:5] #matrix indexing

summary(data[1:5, 1:5]) #summary statistics depending on data class
class(data$X)

library(ggplot2)
ggplot(data = data) +
  geom_point(aes(x=x_centroid, y = y_centroid), size = 0.1) + theme_classic()

ggplot(data = data)+
  geom_histogram(aes(x=area)) + theme_classic()

#### Lesson 2 ( 1-27-23): exploring spatial transcriptomics data ####

## store gene expression in dataframe for convenience
gexp <- data[,5:ncol(data)]

gexp[1:5,1:5] 

## each row is a cell and each column is gene
dim(gexp)

## total number of genes expressed per cell
length(rowSums(gexp)) ## should return number of cells

## total of copies of gene for each gene
length(colSums(gexp)) ## should return number of genes

##set graphical parameters so subsequent figures will be drawn in 1 x 2
par(mfrow=c(1,2))
hist(rowSums(gexp), main = 'total genes per cell')
hist(colSums(gexp), main = 'total copies of each gene across all cells')

##log transform gene expression
mat <- log10(gexp+1)

##set graphical parameters so subsequent figures will be drawn in 1 x 2
par(mfrow=c(1,3))
hist(colSums(mat), main = 'log total copies of each gene across all cells')
hist(colSums(gexp), main = 'total copies of each gene across all cells')
hist(log10(colSums(gexp)), main = 'total copies of each gene across all cells') #log scale on x-axis

## determining gene prevalence by totaling how many cells have non-zero expression for each gene
gprev <- colSums(gexp != 0)

par(mfrow=c(1,1))
hist(gprev, main = 'prevalence of genes across cells')

## ordering genes in decreasing order of prevalence
prevalent.genes <- sort(gprev, decreasing = TRUE)
prevalent.genes[1:5]
prevalent.genes[length(prevalent.genes)]

## plotting cell positions colored by which cells have non-zero expression of ERBB2 and CRHBP

p1 <- ggplot(data=data)+
  geom_point(aes(x=x_centroid, y=y_centroid, col = ERBB2 > 0), size=0.1) + theme_classic()

p2 <- ggplot(data=data)+
  geom_point(aes(x=x_centroid, y=y_centroid, col = CRHBP > 0), size=0.1) + theme_classic()

grid.arrange(p1,p2)

## determining variance of gene expression by variance function on each column of gexp
gvar <- apply(gexp,2, var)
length(gvar)

par(mfrow=c(1,1))
hist(gvar, main = 'variance of genes')

## ordering gene in decreasing order of variance
variable.genes <- sort(gvar, decreasing = TRUE)
variable.genes[1:5]
variable.genes[(length(variable.genes)-4):length(variable.genes)]

par(mfrow=c(1,4))
hist(gexp$ERBB2)
hist(gexp$LUM)
hist(gexp$MPO)
hist(gexp$CRHBP)

## scale gene expression
mat <-scale(gexp) #center and divide the (centered) columns of x by their standard deviations

class(mat)
class(gexp)

mat <- data.frame(mat)

## determine variance of scale gene expression
gvar2 <- apply(mat, 2, var)

par(mfrow=c(1,2))
hist(gvar2, main = 'variance of scale genes')
hist(gvar, main = 'variance of genes')

## testing scaling in wrong dimension
mat2 <-t(scale(t(gexp)))
dim(mat2)

mat2[1:5,1:5]

gvar.rev <- apply(mat2, 2, var, na.rm=TRUE)
length(gvar.rev)
gvar.rev[1:5]
hist(gvar.rev)

#### Lesson 4 (02-01-23): PCA ####

mat <-scale(log10(gexp+1), center = FALSE) #scale and log10 transform
mat[1:5,1:5]

## PCA
pcs <- prcomp(mat)

#should be as many PCs as genes
dim(mat)
dim(pcs$x)

## scree plot: look at how much variance is explained by first n PCs, n=10
par(mfrow=c(1,1))
plot(1:10, pcs$sdev[1:10], type = 'l')

## look at loadings (coefficients of linear combination of genes for each PC)
head(sort(abs(pcs$rotation[,1]), decreasing = TRUE)) #highest absolute loading, genes that contribute the most to PC1
head(sort(abs(pcs$rotation[,2]), decreasing = TRUE)) #genes that contribute most PC2

## sum of loadings^2 should total to 1
sum((pcs$rotation[,1])^2)

## data visualization to explore the first two PCs
df <- data.frame(pcs$x[,1:2], gene=mat[,'KRT8'])

ggplot(data = df, aes(x=PC1, y= PC2, col=gene)) +
  geom_point()+ theme_classic()+
  scale_color_gradient(low='lightgrey',  high = 'red')

#### Lesson 5 (02-03-23): tSNE ####
library(Rtsne)
set.seed(0) #reproducibility
emb <- Rtsne(mat, check_duplicates = FALSE)
dim(mat)
dim(emb$Y)

df2 <- data.frame(emb$Y, gene=mat[,'KRT8'])

ggplot(data = df2, aes(x=X1, y= X2, col=gene)) +
  geom_point()+ theme_classic()+
  scale_color_gradient(low='lightgrey',  high = 'red')
