library(ggplot2)
library(dplyr)
x <- c(0,1/2,sqrt(2)/2,sqrt(3)/2,1,sqrt(3)/2,sqrt(2)/2,1/2,0)
y <- c(1,sqrt(3)/2, sqrt(2)/2, 1/2, 0, -1/2, -sqrt(2)/2, -sqrt(3)/2, -1)

x <- c(x, x/2)
y <- c(y, y/2)

df <- data.frame(x = x, y = y)
df

df %>% ggplot(aes(x, y )) + geom_point()

df$x[1:5]

x1c1 <- mean(df$x[1:5])
x2c1 <- mean(df$y[1:5])
x1c2 <- mean(df$x[6:10])
x2c2 <- mean(df$y[6:10])

nl_cl <- matrix(c(x1c1, x1c2, x2c1, x2c2), nrow = 2)
nl_X <- as.matrix(df)
nl_d <- c(rep(0, times = 9), rep(1, times = 9))

nl_m <- amnfis(nl_X, nl_d, nl_cl)
nl_f <- amnfis.simulate(nl_m, nl_X, nl_cl)

transform_output(nl_f)
