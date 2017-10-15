# Ideation challenge code for InBIx-2017
# Classfication of Boruta selected features using random forest
# DBT-BIF team
# Author: Deepak Sharma
# Date: 25 September, 2017


# Install packages if not installed ---------------------------------------

# install.packages("randomForest")
# install.packages("foreign")
# install.packages("caret")

# Import required libraries -----------------------------------------------

library(randomForest)
library(foreign)
library(caret)

# Parameter configuration -------------------------------------------------

k <- 5                    # K-cross validation
nt <- 3000                # Number of trees
mt <- 656                 # Number of mtry
mn <- 450                 # Maximum no. of nodes
ss <- c(300,300)          # Sample sizes 
cwt <- c(0.385,0.614)     # Class weights
oobt <- 1000              # Out of bag times
ns <- 2                   # Node size

# Import and data processing ----------------------------------------------

features <- read.arff("data/boruta_selected_features.arff") # Read data file
features <- features[sample(nrow(features)),]   # Shuffling of data records
indexes = sample(1:nrow(features), size=1/k*nrow(features)) # Indices
test = features[indexes,]   # Test data
train = features[-indexes,] # Training data

# Data summary ------------------------------------------------------------

cat("No of missing values : ",sum(is.na(features)),"\n")
cat("No. of variables : ",vars <- ncol(features),"\n")
cat("Training cases : ",nrow(train),"\n")
cat("Test cases : ",nrow(test),"\n") 
cat("Positive cases :",sum(features[,657]==1),"\n")
cat("Negative cases :",sum(features[,657]==0),"\n")

# Random Forest -----------------------------------------------------------

attach(train)

# Training case accuracy --------------------------------------------------

tr_forest <- randomForest(output ~., data = train,
          ntree=nt, mtry=mt,importance=TRUE, proximity=TRUE,
          maxnodes=mn,sampsize=ss,classwt=cwt,
          keep.forest=TRUE,oob.prox=TRUE,oob.times= oobt,
          replace=TRUE,nodesize=ns, do.trace=1
          )
tra <- unname(confusionMatrix(as.table(tr_forest$confusion[,-3]))$overall[1])
cat("Training accuracy is :",tra,"\n")

# Test case accuracy ------------------------------------------------------

ts_forest <- randomForest(output ~.,
          data = train, xtest=test[,-vars], ytest=test[,vars],
          ntree=nt, mtry=mt,importance=TRUE, proximity=TRUE,
          maxnodes=mn,sampsize=ss,classwt=cwt,
          keep.forest=TRUE,oob.prox=TRUE,oob.times= oobt,
          replace=TRUE,nodesize=ns, do.trace=1
          )
tsa <- unname(confusionMatrix(as.table(ts_forest$confusion[,-3]))$overall[1])
cat("Testing case accuracy is :",tsa,"\n")

detach(train)

# 5-cross validation ------------------------------------------------------

data <- features
data <- data[sample(nrow(data)),]
folds <- cut(seq(1,nrow(data)),breaks = k,labels = FALSE)
accur_log <- matrix(0,k)
for(i in 1:k){
  indices <- which(folds==i,arr.ind = TRUE)
  test_d <- data[indices,]
  train_d <- data[-indices,]
  attach(train_d)
  cat("Running ",i,"/",k,"fold\n")
  cv_forest <- randomForest(output ~.,
              data = train_d, xtest=test[,-vars], ytest=test[,vars],
              ntree=nt, mtry=mt,importance=TRUE, proximity=TRUE,
              maxnodes=mn,sampsize=ss,classwt=cwt,
              keep.forest=TRUE,oob.prox=TRUE,oob.times= oobt,
              replace=TRUE,nodesize=ns, do.trace=1
              )
  accur_log[i] <- unname(confusionMatrix(as.table(cv_forest$confusion[,-3]))$overall[1])
  detach(train_d)
  Sys.sleep(1)
}
cva <- mean(accur_log)
cat("Average accuracy after ",k,"cross validation is :",cva,"\n")

# Brief report ------------------------------------------------------------

out <- matrix(0,3,2)
cat(out[1,1] <- "Training accuracy : ",out[1,2] <- tra,"\n")
cat(out[2,1] <- "Testing accuracy : ",out[2,2] <- tsa,"\n")
cat(out[3,1] <- "Cross validation accuracy : ",out[3,2] <- cva,"\n")
write.csv(out,"Result.csv")
detach(train)
