


library(ggplot2)
library(Matrix)
library(MASS)

source("l2_map.R")

dataset = "clusters"


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


mat_data <- function(d1=5,d2=5,d3=5){
  X = as.matrix(expand.grid(x1=1:d1,
                            x2=1:d2))
  P = matrix(rnorm(2*d3),2,d3)
  Z = X %*% P
  
  out = list(X=X,
             Z=Z)
  return(out)
}

cube_data <- function(){
  X = as.matrix(expand.grid(x1=c(-3,3),
                            x2=c(-3,3),
                            x3=c(-3,3)))
  out = list(X=X)
  return(out)
}



set.seed(123123)
if(dataset=="clusters"){
  X <- clusters_data()
}else{
  X <- read.csv("mnist_data.csv")
  X <- as.data.frame(X)
  X <- X[,-785]
}
seeds <- sample(999999,size=50,replace=TRUE)

best_x2 <- -Inf
best_x3 <- -Inf
best_x4 <- -Inf
best_x4_1 <- -Inf
best_x4_2 <- -Inf
best_x4_3 <- -Inf
best_x4_4 <- -Inf
best_x4_5 <- -Inf
best_x5 <- -Inf
best_x5_1 <- -Inf
best_x5_2 <- -Inf
best_x5_3 <- -Inf
best_x5_4 <- -Inf
best_x5_5 <- -Inf
best_x6 <- -Inf
best_x6_1 <- -Inf
best_x6_2 <- -Inf
best_x6_3 <- -Inf
best_x6_4 <- -Inf
best_x6_5 <- -Inf
best_x7 <- -Inf
best_x7_1 <- -Inf
best_x7_2 <- -Inf
best_x7_3 <- -Inf
best_x7_4 <- -Inf
best_x7_5 <- -Inf

x2_star <- 0
x3_star <- 0
x4_star <- 0
x4_1_star <- 0
x4_2_star <- 0
x4_3_star <- 0
x4_4_star <- 0
x4_5_star <- 0
x5_star <- 0
x5_1_star <- 0
x5_2_star <- 0
x5_3_star <- 0
x5_4_star <- 0
x5_5_star <- 0
x6_star <- 0
x6_1_star <- 0
x6_2_star <- 0
x6_3_star <- 0
x6_4_star <- 0
x6_5_star <- 0
x7_star <- 0
x7_1_star <- 0
x7_2_star <- 0
x7_3_star <- 0
x7_4_star <- 0
x7_5_star <- 0

iterative_map <- function(X,d=2,niter=1,niter2=1000,...){
  x_temp <- cmdscale(dist(X),k=2)
  for(i in 1:niter){
    eps = 1/(dist(X) - dist(x_temp))^2
    eps = ifelse(eps>1e+6,1e+6,eps)
    x_temp <- l2_map(X,d=2,W=eps,niter=niter2,...)$X
  }
  return(x_temp)
}
#inverse_p <- function(Z){
#  P = ecdf(Z[Z!=0])
#  eps = 1e-06
#  out <- ifelse(P(Z)!=0,1/P(Z),1/eps)
#  return(out)
#}
#ip_weights <- inverse_p(dist(X))


n = dim(X)[1]-1
X = as.matrix(X)
for(i in 1:1){
  set.seed(seeds[i])
  x2 <- cmdscale(dist(X),k=2)
  x3 <- sammon(dist(X),k = 2)$points
  eps = 1/(dist(X) - dist(x2))^2
  eps = ifelse(eps>1e+6,1e+6,eps)
  #single iteration
  x4 <- l2_map(X,d=2,W=eps,niter=1000,method=0)$X
  x4_1<- l2_map(X,d=2,W=eps,niter=1000,method=1)$X
  x4_2<- l2_map(X,d=2,W=eps,niter=1000,method=2)$X
  x4_3<- l2_map(X,d=2,W=eps,niter=1000,method=3,beta=0.3)$X
  x4_4<- l2_map(X,d=2,W=eps,niter=1000,method=3,beta=0.7)$X
  x4_5 <- l2_map(X,d=2,W=eps,niter=1000,method=4)$X
  
  #normal
  x5 <- l2_map(X,d=2,niter=1000,method=0)$X
  x5_1<- l2_map(X,d=2,niter=1000,method=1)$X
  x5_2<- l2_map(X,d=2,niter=1000,method=2)$X
  x5_3<- l2_map(X,d=2,niter=1000,method=3,beta=0.3)$X
  x5_4<- l2_map(X,d=2,niter=1000,method=3,beta=0.7)$X
  x5_5 <- l2_map(X,d=2,niter=1000,method=4)$X
  
  #iterative
  x6 <- iterative_map(X=X,niter=10,method=0)
  x6_1 <- iterative_map(X=X,niter=10,method=1)
  x6_2<- iterative_map(X=X,niter=10,method=2)
  x6_3<- iterative_map(X=X,niter=10,method=3,beta=0.3)
  x6_4 <- iterative_map(X=X,niter=10,method=3,beta=0.7)
  x6_5<- iterative_map(X=X,niter=10,method=4)
  
  #inverse probability
  #x7 <- l2_map(X,d=2,W=ip_weights,niter=1000,method=0)$X
  #x7_1 <- l2_map(X,d=2,W=ip_weights,niter=1000,method=1)$X
  #x7_2 <- l2_map(X,d=2,W=ip_weights,niter=1000,method=2)$X
  #x7_3 <- l2_map(X,d=2,W=ip_weights,niter=1000,method=3,beta=0.3)$X
  #x7_4 <- l2_map(X,d=2,W=ip_weights,niter=1000,method=3,beta=0.7)$X
  #x7_5 <- l2_map(X,d=2,W=ip_weights,niter=1000,method=4)$X
  
  rk2 <- range_kept(X1=x2,X2=X,k=n)
  rk3 <- range_kept(X1=x3,X2=X,k=n)
  
  rk4 <- range_kept(X1=x4,X2=X,k=n)
  rk4_1 <- range_kept(X1=x4_1,X2=X,k=n)
  rk4_2 <- range_kept(X1=x4_2,X2=X,k=n)
  rk4_3 <- range_kept(X1=x4_3,X2=X,k=n)
  rk4_4 <- range_kept(X1=x4_4,X2=X,k=n)
  rk4_5 <- range_kept(X1=x4_5,X2=X,k=n)
  
  rk5 <- range_kept(X1=x5,X2=X,k=n)
  rk5_1 <- range_kept(X1=x5_1,X2=X,k=n)
  rk5_2 <- range_kept(X1=x5_2,X2=X,k=n)
  rk5_3 <- range_kept(X1=x5_3,X2=X,k=n)
  rk5_4 <- range_kept(X1=x5_4,X2=X,k=n)
  rk5_5 <- range_kept(X1=x5_5,X2=X,k=n)
  
  rk6 <- range_kept(X1=x6,X2=X,k=n)
  rk6_1 <- range_kept(X1=x6_1,X2=X,k=n)
  rk6_2 <- range_kept(X1=x6_2,X2=X,k=n)
  rk6_3 <- range_kept(X1=x6_3,X2=X,k=n)
  rk6_4 <- range_kept(X1=x6_4,X2=X,k=n)
  rk6_5 <- range_kept(X1=x6_5,X2=X,k=n)
  
  #rk7 <- range_kept(X1=x7,X2=X,k=n)
  #rk7_1 <- range_kept(X1=x7_1,X2=X,k=n)
  #rk7_2 <- range_kept(X1=x7_2,X2=X,k=n)
  #rk7_3 <- range_kept(X1=x7_3,X2=X,k=n)
  #rk7_4 <- range_kept(X1=x7_4,X2=X,k=n)
  #rk7_5 <- range_kept(X1=x7_5,X2=X,k=n)
  
  if(rk2>best_x2){
    best_x2 = rk2
    x2_star = x2
  }
  if(rk3>best_x3){
    best_x3 = rk3
    x3_star = x3
  }
  if(rk4>best_x4){
    best_x4 = rk4
    x4_star = x4
  }
  if(rk4_1>best_x4_1){
    best_x4_1 = rk4_1
    x4_1_star = x4_1
  }
  if(rk4_2>best_x4_2){
    best_x4_2 = rk4_2
    x4_2_star = x4_2
  }
  if(rk4_3>best_x4_3){
    best_x4_3 = rk4_3
    x4_3_star = x4_3
  }
  if(rk4_4>best_x4_4){
    best_4_x4 = rk4_4
    x4_4_star = x4_4
  }
  if(rk4_4>best_x4_4){
    best_x4_4 = rk4_4
    x4_4_star = x4_4
  }
  if(rk4_5>best_x4_5){
    best_x4_5 = rk4_5
    x4_5_star = x4_5
  }
  
  
  
  if(rk5>best_x5){
    best_x5 = rk5
    x5_star = x5
  }
  if(rk5_1>best_x5_1){
    best_x5_1 = rk5_1
    x5_1_star = x5_1
  }
  if(rk5_2>best_x5_2){
    best_x5_2 = rk5_2
    x5_2_star = x5_2
  }
  if(rk5_3>best_x5_3){
    best_x5_3 = rk5_3
    x5_3_star = x5_3
  }
  if(rk5_4>best_x5_4){
    best_5_x4 = rk5_4
    x5_4_star = x5_4
  }
  if(rk5_5>best_x5_5){
    best_x5_5 = rk5_5
    x5_5_star = x5_5
  }
  
  
  if(rk6>best_x6){
    best_x6 = rk6
    x6_star = x6
  }
  if(rk6_1>best_x6_1){
    best_x6_1 = rk6_1
    x6_1_star = x6_1
  }
  if(rk6_2>best_x6_2){
    best_x6_2 = rk6_2
    x6_2_star = x6_2
  }
  if(rk6_3>best_x6_3){
    best_x6_3 = rk6_3
    x6_3_star = x6_3
  }
  if(rk6_4>best_x6_4){
    best_6_x4 = rk6_4
    x6_4_star = x6_4
  }
  if(rk6_5>best_x6_5){
    best_x6_5 = rk6_5
    x6_5_star = x6_5
  }
  
  
  #if(rk7>best_x7){
  #  best_x7 = rk7
  #  x7_star = x7
  #}
  #if(rk7_1>best_x7_1){
  #  best_x7_1 = rk7_1
  #  x7_1_star = x7_1
  #}
  #if(rk7_2>best_x7_2){
  #  best_x7_2 = rk7_2
  #  x7_2_star = x7_2
  #}
  #if(rk7_3>best_x7_3){
  #  best_x7_3 = rk7_3
  #  x7_3_star = x7_3
  #}
  #if(rk7_4>best_x7_4){
  #  best_7_x4 = rk7_4
  #  x7_4_star = x7_4
  #}
  #if(rk7_5>best_x7_5){
  #  best_x7_5 = rk7_5
  #  x7_5_star = x7_5
  #}
  
  print("Done")
}

method <- c("MDS","Sammon",rep("Inverse Square",6),rep("Iterative Inverse Square",6),rep("L2 Map",6))
stepsize <- c("-","-",rep(c("0.1","1/k","1/sqrt(k)","0.1","0.1","0.1"),3))
descent_method <- c("-","-",rep(c("-","-","-","beta=0.3","beta=0.7","adagrad"),3))
eval <- c(rk2,rk3,rk4,rk4_1,rk4_2,rk4_3,rk4_4,rk4_5,rk5,rk5_1,rk5_2,rk5_3,rk5_4,rk5_5,rk6,rk6_1,rk6_2,rk6_3,rk6_4,rk6_5)
df <- data.frame(method=method,
                 stepsize=stepsize,
                 descent_method=descent_method,
                 eval=eval)
df[,"seed"] <- seeds[1]

write.csv(df,filename=paste(dataset,seeds[1],".csv",sep=""))

