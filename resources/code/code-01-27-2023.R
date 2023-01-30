?read.csv
data <- read.csv('data/charmander.csv.gz')
dim(data)
data[1:5,1:5]

## not necessary 
pos <- data[,2:3]
area <- data[,4]
gexp <- data[,5:ncol(data)]
rownames(pos) <- names(area) <- rownames(gexp) <- data[,1]
head(pos)
head(area)

## your homework will be due Monday midnight 
## is to make a data visualization
## please use ggplot
library(ggplot2)
ggplot(data = data) +
  geom_point(aes(x = x_centroid, y = y_centroid,
                 col = MS4A1 > 0), size=0.1) +
  theme_classic() 

## relationship between MS4A1 vs. GZMB
ggplot(data = data) +
  geom_point(aes(x = MS4A1, y = GZMB), size=0.1) +
  theme_classic() 

## Ryan's question: 
## using multiple matrices
## using variables
g <- 'MS4A1'
ggplot(data = pos) +
  geom_bin_2d(aes(x =  x_centroid, y = gexp[,g])) +
  theme_classic() 







