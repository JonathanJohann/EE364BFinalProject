

install.packages("dslabs")
library(dslabs)
df <- dslabs::read_mnist()

set.seed(123123)
train <- df$train
indices <- c()
for(i in 1:10){
  indices <- c(indices,sample(which(train$labels==(i-1)),size=50,replace=FALSE))
}

mnist <- train$images[indices,]
mnist <- cbind(mnist,train$labels[indices])
colnames(mnist) <- c(paste("pixel",1:784,sep=""),"label")
write.csv(mnist,file="mnist_data.csv")
