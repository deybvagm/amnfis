library(ggplot2)
library(dplyr)
x <- c(0,1/2,sqrt(2)/2,sqrt(3)/2,1,sqrt(3)/2,sqrt(2)/2,1/2,0)
y <- c(1,sqrt(3)/2, sqrt(2)/2, 1/2, 0, -1/2, -sqrt(2)/2, -sqrt(3)/2, -1)

x <- c(x, x/2)
y <- c(y, y/2)

df <- data.frame(x = x, y = y)
df

df %>% ggplot(aes(x, y )) + geom_point()
