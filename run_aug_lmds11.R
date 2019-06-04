


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
  squared=c(TRUE),#,FALSE),
  t=c(10^-2),#10^seq(-3,0,0.5),
  k=c(20),#c(3,5,10,20,40,80,100,200),
  stepsize=c(1e-3)
)

for(i in 1:iterations){
  seed <- sample(1000000,size=1)
  
  range_kepts <- c()
  
  for(q in 1:dim(config)[1]){
    cfg <- config[q,]
    set.seed(seed)
    x1 <- rank_constrained_projection(X,lmds=TRUE,augment=TRUE,
                                      squared=cfg$squared,stepsize=cfg$stepsize,
                                      k=cfg$k,tau=cfg$t)
    rk <- range_kept(x1,X,(dim(X)[1]-1))
    range_kepts <- c(range_kepts,rk)
    intermediate <- data.frame(method="aug_lmds",
                               dataset=dataset,
                               seed=seed,
                               eval=rk)
    intermediate <- cbind(intermediate,cfg)
    write.csv(intermediate,file=paste("rank_lmds_constrained_intermediate",dataset,seed,".csv",sep=""))
  }
  
  result <- data.frame(method="aug_lmds",
                       eval=range_kepts)
  result[,"seed"] <- seed
  result[,"dataset"] <- dataset
  result <- cbind(result,config)
  write.csv(result,file=paste("rank_lmds_constrained",dataset,seed,".csv",sep=""))
  
}
