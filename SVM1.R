levels(new_house$TimeToSubway)<-1:3
levels(new_house$TimeToBusStop)<-1:5
new_house[, 6] <- as.numeric(as.character( new_house[, 6] ))
new_house[, 7] <- as.numeric(as.character( new_house[, 7] ))
levels(new_house$Year_Trans)<-c('below4years','above4years','above4years')
levels(new_house$Year_Trans)<-0:1
new_house[,20]<- as.numeric(as.character(new_house[,20]))

intrain <- createDataPartition(y = new_house$Year_Trans, p= 0.35, list = FALSE)
train<-new_house[intrain,]
test<-new_house[-intrain,]
train[["Year_Trans"]] = factor(train[["Year_Trans"]])
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
svm_Linear <- train(Year_Trans ~., data = train, method = "svmLinear",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)
test_pred <- predict(svm_Linear, newdata = test)
confusionMatrix(table(test_pred, test$Year_Trans))

