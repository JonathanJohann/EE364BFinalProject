range_kept <- function(X1,X2,k){
  
  #' Since we consistently want to see a spectrum of
  #' how many nearest neighbors are kept, the focus
  #' of this function is to capture 1 to k of the
  #' percent of nearest neighbors kept.
  
  n <- dim(X1)[1]
  
  distances1 <- dist(X1)
  distances2 <- dist(X2)
  
  rank1 <- apply(distances1,2,sort)
  rank2 <- apply(distances2,2,sort)
  
  kept <- matrix(0,nrow=(k-1),ncol=1)
  
  for(i in 1:(k-1)){
    neighborhood1 <- ifelse(distances1>rank1[(i+1),],0,1)
    neighborhood2 <- ifelse(distances2>rank2[(i+1),],0,1)
    percent_neighbors_same <- (sum(neighborhood1*neighborhood2)-n)/(i*n)
    kept[i] <- percent_neighbors_same
  }
  return(sum(kept))
}
#### Calculate the distance matrix of X: n by p ####
dist <- function(X)
{
  X2 <- apply (X, 1, function(z)sum(z^2))
  dist2 <- outer(X2, X2, "+")- 2*X%*%t(X)
  dist2 <- ifelse(dist2<0,0,dist2)
  return (sqrt(dist2))
}


rank_constrained_projection <- function(X,d=2,niter=500,print_iter=TRUE,squared=FALSE,
                              augment=FALSE,beta=0.3,method=0,lambda=1,
                              mu=1,nu=0,tau=1,stepsize=1e-4,lmds=FALSE,k=6){
  X <- as.matrix(X)
  N = dim(X)[1]
  P = matrix(rnorm(n=dim(X)[2]*d),nrow=dim(X)[2],ncol=d)
  X1 <- X %*% P
  D1 <- dist(X1)
  D0 <- dist(X)
  Ri <- apply(D0,1,rank)
  Ro <- apply(D1,1,rank)
  it <- 0
  jt <- 1
  last_update <- 0
  cache = matrix(0,nrow=N,ncol=d)
  Do <- dist(X)
  n <- nrow(Do)
  Dnu <- Do^(0)
  Dnulam <- Do^(1)
  diag(Dnu) <- 0
  diag(Dnulam) <- 0
  if(lmds==TRUE){
    Do <- dist(X)
    n <- nrow(Do)
    Daux <- apply(Do,2,sort)[k+1,]# choose the kth smallest distance expect 0.
    Inb <- ifelse(Do>Daux, 0,1)
    k.v <-  apply(Inb, 1, sum)
    k <- (sum(k.v)-n)/n
    Inb.sum <- matrix(rep(k.v, n),ncol=n)
    Mka <- 0
    Inb1 <- pmax(Inb,t(Inb)) # expanding the neighbors for symmetry.
    Dnu <- ifelse(Inb1==1, Do^nu, 0)
    Dnulam <- ifelse(Inb1==1, Do^(nu+1/lambda), 0)
    diag(Dnu) <- 0
    diag(Dnulam) <- 0
    cc <- (sum(Inb1)-n)/n/n*median(Dnulam[Dnulam!=0])
    t <- tau*cc
  }
  while(it<niter){
    updated <- FALSE
    i <- (as.integer(jt) %% N) + 1
    j <- ((as.integer(jt/N))%%N)+1
    j2 <- ((as.integer(jt/(N^2)))%%N)+1
    close_indices <- which(Ri[,i]==j)
    far_indices <- which(Ri[,i]==j2)
    if(j>j2){
      if(length(close_indices)>0){
        for(q in 1:length(close_indices)){
          if(length(far_indices)>0){
            for(qq in 1:length(far_indices)){
              kc <- close_indices[q]
              kf <- far_indices[qq]
              ref <- X1[i,]
              clo <- X1[kf,]
              far <- X1[kc,]
              l <- 1
              grad <- matrix(0,nrow=N,ncol=d)
              if(squared){
                grad[i,] <- 2*(ref-clo) - 2*(ref-far)
                grad[kf,] <- -2*(ref-clo)
                grad[kc,] <- 2*(ref-far)
              }
              else{
                grad[i,] <- 2*(ref-clo)/sqrt(sum((ref-clo)^2)) - 2*(ref-far)/sqrt(sum((ref-far)^2))
                grad[kf,] <- -2*(ref-clo)/sqrt(sum((ref-clo)^2))
                grad[kc,] <- 2*(ref-far)/sqrt(sum((ref-far)^2))
              }
              if(j2<j){
                print(i)
                print(kc)
                print(kf)
                print(ref)
                print(clo)
                print(far)
                while((sum((ref-clo)^2)>sum((ref-far)^2))&(l<10)){
                  alpha <- 0.00001/l
                  l <- l + 1
                  X1 <- X1 - alpha * grad
                  updated <- TRUE
                  print("Updating...")
                }
              }
              if(augment){
                X0 <- X1
                D1mu2 <- D1^(mu-2)
                diag(D1mu2) <- 0
                D1mulam2 <- D1^(mu+1/lambda-2)
                diag(D1mulam2) <- 0
                if(lmds==TRUE){
                  M <- Dnu*D1mulam2-D1mu2*(Dnulam+t*(!Inb1))
                } else{
                  M <- (Dnu*D1mulam2-D1mu2*Dnulam)/ifelse(D0==0,1,D0)
                }
                E <- matrix(rep(1,n*d),n,d)
                Grad <- X0*(M%*%E)-M%*%X0
                normgrad <- (norm(X0)/norm(Grad))*Grad
                if(method==0){
                  #constant
                  X1 <- X0 - stepsize*normgrad  
                } else if(method == 1){
                  #1/k
                  stepsize <- 1/i
                  X1 <- X0 - stepsize*normgrad  
                } else if(method == 2){
                  #1/sqrt(k)
                  stepsize <- 1/sqrt(i)
                  X1 <- X0 - stepsize*normgrad
                } else if(method == 3){
                  #heavy ball
                  X1 <- X0 - stepsize*normgrad + beta * (X0 - lagX0)
                  lagX0 <- X0
                } else if(method == 4){
                  #adagrad
                  cache = cache + normgrad^2
                  X1 <- X0 - stepsize*normgrad / (sqrt(cache) + epsilon)
                }
                D1 <- dist(X1)
                D1mulam <- D1^(mu+1/lambda)
                diag(D1mulam) <- 0
                D1mu <- D1^mu
                diag(D1mu) <- 0
              }
              if(updated){
                last_update <- jt
                it<- it + 1
                if(print_iter){
                  print(paste("Current Update -- ",it,sep=""))
                  if(it%%100==1){
                    plot(X1)
                    print(is.nan(X1))
                  }
                }
              }
              if(jt-last_update>(N*N*(N+1))){
                output <- list()
                output$X <- X1
                return(X1)
              }
            }
          }
        }
      }
    }
    jt <- jt + 1
    print(jt)
  }
  output <- list()
  output$X <- X1
  return(X1)
}




l2_map <- function(X,k=1,d=3,lambda=1,mu=1,nu=0,tau=1,niter=500,W=NULL,method=0,beta=0.5,epsilon=1e-6,rankings=TRUE){
  Do <- dist(X)
  n <- nrow(Do)
  Dnu <- Do^(0)
  Dnulam <- Do^(1)
  diag(Dnu) <- 0
  diag(Dnulam) <- 0
  P = matrix(rnorm(dim(X)[2]*d),nrow=dim(X)[2],ncol=d)
  X1 <- X %*% P
  D1 <- dist(X1)
  X1 <- X1*norm(Do)/norm(D1)
  stepsize <-0.1
  i <- 1
  lagX0 <- 0
  cache = matrix(0,nrow=n,ncol=d)
  while (i < niter)
  {
    X0 <- X %*% P
    D1mu2 <- D1^(mu-2)
    diag(D1mu2) <- 0
    D1mulam2 <- D1^(mu+1/lambda-2)
    diag(D1mulam2) <- 0
    M <- Dnu*D1mulam2-D1mu2*Dnulam
    E <- matrix(rep(1,n*d),n,d)
    if(rankings){
      knnX <- apply(D1,1,rank)
      knnY <- apply(Do,1,rank)
      R <- ifelse(knnX==knnY,1,0)
      M <- M*R
    }  
    if(!is.null(W)){
      Grad <- X0*((W*M)%*%E)-(W*M)%*%X0
    }
    else{
      Grad <- X0*(M%*%E)-M%*%X0
    }
    normgrad <- (norm(X0)/norm(Grad))*Grad
    if(method==0){
      #constant
      X1 <- X0 - stepsize*normgrad  
    } else if(method == 1){
      #1/k
      stepsize <- 1/i
      X1 <- X0 - stepsize*normgrad  
    } else if(method == 2){
      #1/sqrt(k)
      stepsize <- 1/sqrt(i)
      X1 <- X0 - stepsize*normgrad
    } else if(method == 3){
      #heavy ball
      X1 <- X0 - stepsize*normgrad + beta * (X0 - lagX0)
      lagX0 <- X0
    } else if(method == 4){
      #adagrad
      cache = cache + normgrad^2
      X1 <- X0 - stepsize*normgrad / (sqrt(cache) + epsilon)
    }
    D1 <- dist(X1)
    D1mulam <- D1^(mu+1/lambda)
    diag(D1mulam) <- 0
    D1mu <- D1^mu
    diag(D1mu) <- 0
    P = ginv(t(X)%*%X)%*%t(X)%*%X1
    i <- i+1
    if(i%%100 == 0){
      print(i)
    }
  }
  result <- list()
  result$X <- X1
  result$P <- P
  return(result)
}
