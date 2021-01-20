
levels(new_house$TimeToSubway)<-1:3
levels(new_house$TimeToBusStop)<-1:5
new_house[, 6] <- as.numeric(as.character( new_house[, 6] ))
new_house[, 7] <- as.numeric(as.character( new_house[, 7] ))
levels(new_house$Year_Trans)<-c('below4years','above4years','above4years')
norm <-function(x) { (x -min(x))/(max(x)-min(x))   }
house_norm <- as.data.frame(lapply(new_house[,c(1:19)], norm))
round(prop.table(table(new_house$Year_Trans)) * 100, digits = 1)
train <- house_norm[1:2094,] 
test<- house_norm[2095:5878,]

train_target<-new_house[1:2094,20]
test_target<-new_house[2095:5878,20]

library(class)

knn.45 <- knn(train, test, cl=train_target, k=45)
knn.46 <- knn(train, test, cl=train_target, k=46)

ACC.45 <- 100 * sum(test_target == knn.45)/NROW(test_target)
ACC.46 <- 100 * sum(test_target == knn.46)/NROW(test_target)

library(caret)

confusionMatrix(table(knn.45,test_target))
confusionMatrix(table(knn.46,test_target))


i=1
k.optm=1
for (i in 1:47){
knn.mod <- knn(train=train, test=test, cl=train_target, k=i)
k.optm[i] <- 100 * sum(test_target == knn.mod)/NROW(test_target)
k=i
cat(k,'=',k.optm[i],'')
   }
plot(k.optm, type="b", xlab="K- Value",ylab="Accuracy level")

