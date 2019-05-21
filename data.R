


library(ggplot2)
library(Matrix)
library(MASS)

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
df <- cube_data()#mat_data()
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
    x_temp <- l2_map(df$X,d=2,W=eps,niter=niter2,...)$X
  }
  return(x_temp)
}
inverse_p <- function(Z){
  P = ecdf(Z[Z!=0])
  eps = 1e-06
  out <- ifelse(P(Z)!=0,1/P(Z),1/eps)
  return(out)
}
ip_weights <- inverse_p(dist(X))


n = dim(df$X)[1]-1
X = as.matrix(df$X)
for(i in 1:20){
  set.seed(seeds[i])
  x2 <- cmdscale(dist(X),k=2)
  x3 <- sammon(dist(X),k = 2)$points
  eps = 1/(dist(X) - dist(x2))^2
  eps = ifelse(eps>1e+6,1e+6,eps)
  #single iteration
  x4 <- l2_map(df$X,d=2,W=eps,niter=1000,method=0)$X
  x4_1<- l2_map(df$X,d=2,W=eps,niter=1000,method=1)$X
  x4_2<- l2_map(df$X,d=2,W=eps,niter=1000,method=2)$X
  x4_3<- l2_map(df$X,d=2,W=eps,niter=1000,method=3,beta=0.3)$X
  x4_4<- l2_map(df$X,d=2,W=eps,niter=1000,method=3,beta=0.7)$X
  x4_5 <- l2_map(df$X,d=2,W=eps,niter=1000,method=4)$X
  
  #normal
  x5 <- l2_map(df$X,d=2,niter=1000,method=0)$X
  x5_1<- l2_map(df$X,d=2,niter=1000,method=1)$X
  x5_2<- l2_map(df$X,d=2,niter=1000,method=2)$X
  x5_3<- l2_map(df$X,d=2,niter=1000,method=3,beta=0.3)$X
  x5_4<- l2_map(df$X,d=2,niter=1000,method=3,beta=0.7)$X
  x5_5 <- l2_map(df$X,d=2,niter=1000,method=4)$X
  
  #iterative
  x6 <- iterative_map(X=df$X,niter=10,method=0)
  x6_1 <- iterative_map(X=df$X,niter=10,method=1)
  x6_2<- iterative_map(X=df$X,niter=10,method=2)
  x6_3<- iterative_map(X=df$X,niter=10,method=3,beta=0.3)
  x6_4 <- iterative_map(X=df$X,niter=10,method=3,beta=0.7)
  x6_5<- iterative_map(X=df$X,niter=10,method=4)
  
  #inverse probability
  x7 <- l2_map(df$X,d=2,W=ip_weights,niter=1000,method=0)$X
  x7_1 <- l2_map(df$X,d=2,W=ip_weights,niter=1000,method=1)$X
  x7_2 <- l2_map(df$X,d=2,W=ip_weights,niter=1000,method=2)$X
  x7_3 <- l2_map(df$X,d=2,W=ip_weights,niter=1000,method=3,beta=0.3)$X
  x7_4 <- l2_map(df$X,d=2,W=ip_weights,niter=1000,method=3,beta=0.7)$X
  x7_5 <- l2_map(df$X,d=2,W=ip_weights,niter=1000,method=4)$X
  
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
  
  rk7 <- range_kept(X1=x7,X2=X,k=n)
  rk7_1 <- range_kept(X1=x7_1,X2=X,k=n)
  rk7_2 <- range_kept(X1=x7_2,X2=X,k=n)
  rk7_3 <- range_kept(X1=x7_3,X2=X,k=n)
  rk7_4 <- range_kept(X1=x7_4,X2=X,k=n)
  rk7_5 <- range_kept(X1=x7_5,X2=X,k=n)
  
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
  
  
  if(rk7>best_x7){
    best_x7 = rk7
    x7_star = x7
  }
  if(rk7_1>best_x7_1){
    best_x7_1 = rk7_1
    x7_1_star = x7_1
  }
  if(rk7_2>best_x7_2){
    best_x7_2 = rk7_2
    x7_2_star = x7_2
  }
  if(rk7_3>best_x7_3){
    best_x7_3 = rk7_3
    x7_3_star = x7_3
  }
  if(rk7_4>best_x7_4){
    best_7_x4 = rk7_4
    x7_4_star = x7_4
  }
  if(rk7_5>best_x7_5){
    best_x7_5 = rk7_5
    x7_5_star = x7_5
  }
  
  print("Done")
}





plot_cube <- function(xi,ref_cube,title){
  ref_cube <- as.matrix(ref_cube)
  line_segments <- ifelse(dist(ref_cube)==min(dist(ref_cube)[dist(ref_cube)!=0]),1,0)
  print(line_segments)
  xi <- data.frame(xi)
  colnames(xi) <- c("x","y")
  p <- ggplot(data=xi,aes(x=x,y=y)) + geom_point()
  for(i in 1:dim(line_segments)[1]){
    for(j in i:dim(line_segments)[1]){
      if(line_segments[i,j]==1){
        df <- data.frame(x1 = xi[i,1],y1=xi[i,2],x2=xi[j,1],y2=xi[j,2])
        p <- p + geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2,colour="segment"),data=df)
      }
      
    }
  }
  p <- p + ggtitle(title)
  return(p)
}

p1 <- plot_cube(x2_star,X,paste("MDS : ",rk2,sep=""))
p2 <- plot_cube(x3_star,X,paste("Sammon : ",rk3,sep=""))

p3 <- plot_cube(x4_star,X,paste("Inverse Square : ",rk4,sep=""))
p3_1 <- plot_cube(x4_1_star,X,paste("1/e^2 w/ 1/k : ",rk4_1,sep=""))
p3_2 <- plot_cube(x4_2_star,X,paste("1/e^2 w/ 1/sqrt(k) : ",rk4_2,sep=""))
p3_3 <- plot_cube(x4_3_star,X,paste("1/e^2 w/ heavy ball, beta = 0.3 : ",rk4_3,sep=""))
p3_4 <- plot_cube(x4_4_star,X,paste("1/e^2 w/ heavy ball, beta = 0.7 : ",rk4_4,sep=""))
p3_5 <- plot_cube(x4_5_star,X,paste("1/e^2 w/ adagrad : ",rk4_5,sep=""))


p4 <- plot_cube(x5_star,X,paste("L2 Map : ",rk5,sep=""))
p4_1 <- plot_cube(x5_1_star,X,paste("L2 Map w/ 1/k : ",rk5_1,sep=""))
p4_2 <- plot_cube(x5_2_star,X,paste("L2 Map w/ 1/sqrt(k) : ",rk5_2,sep=""))
p4_3 <- plot_cube(x5_3_star,X,paste("L2 Map w/ heavy ball, beta = 0.3 : ",rk5_3,sep=""))
p4_4 <- plot_cube(x5_4_star,X,paste("L2 Map w/ heavy ball, beta = 0.7 : ",rk5_4,sep=""))
p4_5 <- plot_cube(x5_5_star,X,paste("L2 Map w/ adagrad : ",rk5_5,sep=""))

p5 <- plot_cube(x6_star,X,paste("Iterative Map : ",rk6,sep=""))
p5_1 <- plot_cube(x6_1_star,X,paste("Iterative Map w/ 1/k : ",rk6_1,sep=""))
p5_2 <- plot_cube(x6_2_star,X,paste("Iterative Map w/ 1/sqrt(k) : ",rk6_2,sep=""))
p5_3 <- plot_cube(x6_3_star,X,paste("Iterative Map w/ heavy ball, beta = 0.3 : ",rk6_3,sep=""))
p5_4 <- plot_cube(x6_4_star,X,paste("Iterative Map w/ heavy ball, beta = 0.7 : ",rk6_4,sep=""))
p5_5 <- plot_cube(x6_5_star,X,paste("Iterative Map w/ adagrad : ",rk6_5,sep=""))




p6 <- plot_cube(x7_star,X,paste("Inverse Prob : ",rk7,sep=""))
p6_1 <- plot_cube(x7_1_star,X,paste("Inverse Prob w/ 1/k : ",rk7_1,sep=""))
p6_2 <- plot_cube(x7_2_star,X,paste("Inverse Prob w/ 1/sqrt(k) : ",rk7_2,sep=""))
p6_3 <- plot_cube(x7_3_star,X,paste("Inverse Prob w/ heavy ball, beta = 0.3 : ",rk7_3,sep=""))
p6_4 <- plot_cube(x7_4_star,X,paste("Inverse Prob w/ heavy ball, beta = 0.7 : ",rk7_4,sep=""))
p6_5 <- plot_cube(x7_5_star,X,paste("Inverse Prob w/ adagrad : ",rk7_5,sep=""))






grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] +
                    theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x +
                 theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend, ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}
grid_arrange_shared_legend(p6,p6_1,p6_2,p6_3,p6_4,p6_5,ncol=3,nrow=2,position="right")

