


library(Rtsne)

source("l2_map.R")

dataset = "clusters"

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

seed <- sample(1000000,size=1)

range_kepts <- c()
if(dataset=="clusters"){
  perplexitys <- seq(26,260,26)
}else{
  perplexitys <- seq(16,160,16)
}
for(i in 1:10){
  set.seed(seed)
  fit <- Rtsne::Rtsne(X,perplexity=perplexitys[i],theta=0.0)
  rk <- range_kept(fit$Y,X,(dim(X)[1]-1))
  range_kepts <- c(range_kepts,rk)
  print(i)
}

result <- data.frame(perp=perplexitys,
                     eval=range_kepts)
result[,"seed"] <- seed
result[,"dataset"] <- dataset

write.csv(result,file=paste(dataset,seed,".csv",sep=""))
