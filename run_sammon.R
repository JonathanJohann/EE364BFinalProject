


library(Rtsne)
library(stats)
library(MASS)

source("l2_map.R")

dataset = "clusters"
iterations = 10

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / ifelse(sd(x, na.rm=TRUE)==0,1,sd(x, na.rm=TRUE))
}

generate_clusters <- function(mean_vals,n_per_cluster){
  df <- c()
  labels <- c()
  p <- dim(mean_vals)[2]
  for(i in 1:dim(mean_vals)[1]){
    tmp <- MASS::mvrnorm(n_per_cluster,mu=mean_vals[i,],Sigma=diag(p))
    tmp <- cbind(tmp,rep(i,n_per_cluster))
    df <- rbind(df,tmp)
  }
  return(df)
}

clusters_data <- function(seed=123123){
  set.seed(seed)
  mean_vals <- expand.grid(
    v1 <- c(5,0),
    v2 <- c(5,0),
    v3 <- c(5,0),
    v4 <- c(5,0)
  )
  dat <- generate_clusters(as.matrix(mean_vals),50)
  labels <- dat[,5]
  dat <- dat[,-5]
  return(dat)
}

if(dataset=="clusters"){
  X = as.matrix(clusters_data())
} else {
  df2 <- read.csv("mnist_data.csv") 
  labels <- df2[,"label"]
  mnist <- df2[,-dim(df2)[2]]
  X = as.matrix(mnist)
}

for(i in 1:dim(X)[2]){
  X[,i] <- scale_this(X[,i])
}

iterative_map <- function(X,d=2,niter=1,niter2=1000,...){
  x_temp <- cmdscale(dist(X),k=2)
  for(i in 1:niter){
    eps = 1/(dist(X) - dist(x_temp))^2
    eps = ifelse(eps>1e+6,1e+6,eps)
    x_temp <- l2_map(X,d=2,W=eps,niter=niter2,...)$X
  }
  return(x_temp)
}
# 260 is max for clusters
# 160 is max for MNIST sub sample

for(i in 1:iterations){
  seed <- sample(1000000,size=1)
  
  range_kepts <- c()
  
  
  set.seed(seed)
  x3 <- sammon(dist(X),k = 2)$points
  rk <- range_kept(x3,X,(dim(X)[1]-1))
  
  
  
  result <- data.frame(method="sammon",
                       eval=rk)
  result[,"seed"] <- seed
  result[,"dataset"] <- dataset
  
  write.csv(result,file=paste(dataset,seed,".csv",sep=""))
  
}
