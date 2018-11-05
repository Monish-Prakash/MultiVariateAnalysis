library(readr)
library(rpart)
library(rpart.plot)
library(ROCR)
library(DMwR)
library(mice)
library(VIM)
library(caret)
library(randomForest)

wd = getwd()
if(grepl("charlyo", wd)) {
    setwd("~/Desktop/MVA/Practiques")
} else {
    setwd("~/Documents/UPC-MASTER/2nd Semester/MVA/Homework7")
}
Sys.setlocale("LC_TIME", "en_US")  # OS X, in UTF-8

# 1
#audit <- read_excel("audit.xlsx", na = "NA")
audit <- read_delim("audit.csv", ";", escape_double = FALSE, trim_ws = TRUE)

# 2
audit <- audit[,-c(1, 12)]

#Check num. of missing values
#which(is.na(audit))
aggr_plot <- aggr(audit, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"), cex.numbers=0.2)

filedDataset <- audit
#Convert categorical variables to factor, otherwise mice not working well
cols <- c("Employment", "Education", "Marital", "Occupation", "Gender", "Accounts", "Adjusted")
filedDataset[cols] <- lapply(filedDataset[cols], factor)
#Fill the missing values of the dataset 
predfiledDataset <- mice(filedDataset, m = 5, maxit = 10, method = "cart", seed = 15000)
#predfiledDataset$loggedEvents
filtered_audit <- complete(predfiledDataset, 1)
#Check num. of missing values
which(is.na(filtered_audit))

# 3
#Select test and training data
rows <- nrow(filtered_audit)
training_data <- filtered_audit[1:((2/3)*rows),]
test_data <- filtered_audit[((2/3)*rows):rows,]

# 4
p1 = rpart(Adjusted ~ ., data=training_data,control=rpart.control(cp=0.001, xval=10))
printcp(p1)
plotcp(p1)

p1$cptable = as.data.frame(p1$cptable)
min_xerror = min(p1$cptable$xerror)
xstd = p1$cptable[p1$cptable$xerror == min_xerror, ]$xstd

filtered_cptable = p1$cptable[(min_xerror+xstd) > p1$cptable$xerror, ]
optimal_splits =  filtered_cptable[1, ]$nsplit
optimal_CP =  filtered_cptable[1, ]$CP

data.frame(optimal_CP = optimal_CP, optimal_splits = optimal_splits)

p2 <- prune(p1,cp=optimal_CP)
rpart.plot(p2)

# 5
barplot(p2$variable.importance, cex.names = 0.5)

# 6
#id <- which(!(test_data$Accounts %in% levels(as.factor(training_data$Accounts))))
#test_data$Accounts[id] <- NA
#str(test_data)
#levels(as.factor(predictions))
predictions <- predict(p2, newdata = test_data[,-c(11)], type = "class")

confusionMatrix <- confusionMatrix(predictions, test_data$Adjusted)

TP <- confusionMatrix$table[1,1]
FN <- confusionMatrix$table[1,2]
FP <- confusionMatrix$table[2,1]
TN <- confusionMatrix$table[2,2]

#Compute the accuracy, precision, recall and AUC on the test individuals.
errorRate = (FN + FP)/length(predictions)
acc = ((1-errorRate) *100)

precisionP = (TP / (TP + FP))
precisionN = (TN / (FN + TN))
precision = (((precisionP + precisionN)/2) *100)

recall = (TP / (TP + FN) *100)

predictions1  <-  as.data.frame(predict(p2, newdata = test_data[,-c(11)], type="prob"))
pred <- prediction(predictions1$`1`, test_data$Adjusted)
#roc <- performance(pred,measure="tpr",x.measure="fpr")
#plot(roc, main="ROC curve")
#abline(0,1,col="blue")
auc = performance(pred,"auc")
auc = as.numeric(auc@y.values)

results = data.frame(accuracy = acc, precision = precision, recall = recall, AUC = auc)
results

# 7
audit.rf = randomForest(formula = Adjusted ~ ., data = training_data, 
                        xtest=test_data[,-11], ytest=test_data[,11], 
                        importance=TRUE, na.action=na.fail)

confusion_matrix = as.data.frame(audit.rf$test$confusion)

precision_TRUE =  (confusion_matrix['1', '1']* 100)/sum(confusion_matrix$'1')

precision_FALSE = (confusion_matrix['0','0']* 100)/sum(confusion_matrix$'0')

avg_precision = (precision_TRUE + precision_FALSE)/2

avg_precision

audit.rf[order(audit.rf$importance)]
barplot(as.data.frame(audit.rf$importance)$MeanDecreaseGini, names.arg=row.names(as.data.frame(audit.rf$importance)), main="MeanDecreaseGini", cex.names=0.5)

plot(audit.rf)
error.rate = mean(as.data.frame(audit.rf$err.rate)$OOB)

100* (1-error.rate)

