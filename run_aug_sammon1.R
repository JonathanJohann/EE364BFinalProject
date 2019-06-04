


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


# 260 is max for clusters
# 160 is max for MNIST sub sample


config <- expand.grid(
  niter=1000,
  squared=c(TRUE),#,FALSE
  augment=c(TRUE),
  stepsize=c(1e-5)#1e-1,1e-3,1e-4,
)

for(i in 1:iterations){
  seed <- sample(1000000,size=1)
  
  range_kepts <- c()
  
  for(q in 1:dim(config)[1]){
    cfg <- config[q,]
    set.seed(seed)
    x1 <- rank_constrained_projection(X,squared=cfg$squared,augment=cfg$augment,stepsize=cfg$stepsize)
    rk <- range_kept(x1,X,(dim(X)[1]-1))
    range_kepts <- c(range_kepts,rk)
    intermediate <- data.frame(dataset=dataset,
                               seed=seed,
                               eval=rk)
    intermediate <- cbind(intermediate,cfg)
    write.csv(intermediate,file=paste("rank_constrained_intermediate",dataset,seed,".csv",sep=""))
  }

  result <- data.frame(method="aug_sammon",
                       eval=range_kepts)
  result[,"seed"] <- seed
  result[,"dataset"] <- dataset
  result <- cbind(result,config)
  write.csv(result,file=paste("rank_constrained",dataset,seed,".csv",sep=""))
  
}
