


X <- mnist[,-785]
labels <- mnist[,785]


x1 <- cmdscale(dist(X))
plot(x1,col=labels)

install.packages("Rtsne")
library(Rtsne)

fit <- Rtsne::Rtsne(dist(X),perplexity=150)
plot(fit$Y,col=labels)


library(MASS)

x3 <- MASS::sammon(dist(X))
plot(x3$points,col=labels)



smap <- sequential_map(X,d=2,niter=500,print_iter = TRUE)
head(X$X)
