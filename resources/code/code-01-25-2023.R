?iris
?dim
dim(iris)
class(iris)
?head
head(iris, n=10)
iris[1,]
iris[,1]
iris[1:3,1:3]

## what are the data types?
# 1. categorical data type: species
# 2. quantitative data type: sepal length, sepal width, petal length, petal width
## advanced students: pick a different dataset

## let's make an exploratory data visualization
#install.packages("ggplot2")
library(ggplot2)
?ggplot
ggplot(data = iris) +
  geom_point(aes(x = Sepal.Length, y = Sepal.Width, 
                 col = Petal.Width))

## TODO: make a data visualization
## What do you think of this data visualization?
ggplot(data = iris) +
  geom_bar(stat = "identity", position = "dodge2",
           aes(x = Species, y = Petal.Width)) +
  theme_classic()
## I'm trying to make salient the relationship 
## between petal width and species

## What feature of the data were you trying to explore?
#### association between petal.length versus petal width
#### and whether it differs between species
## What do I do to make your data visualization?

## tutorial used: https://community.rstudio.com/t/insert-regression-model-into-ggplot2/2439
## I worked with Wendy
ggplot(data = iris) +
  geom_point(aes(x = Petal.Length, y = Petal.Width,
                 col = Species)) +
  theme_classic() +
  geom_smooth(method = "lm",
              aes(x = Petal.Length, y = Petal.Width,
              col = Species))



ggplot(data = iris) +
  geom_bin_2d(aes(x = Petal.Length, y = Petal.Width,
                  col = Species)) +
  theme_classic() 


