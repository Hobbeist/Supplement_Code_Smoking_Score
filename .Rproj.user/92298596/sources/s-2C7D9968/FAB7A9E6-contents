#=====================================================================================
# Machine learning-based DNA methylation score for fetal exposure to maternal smoking: 
# development and validation in samples collected from adolescents and adults
# Training:
#     1. Random Forest
#     2. Gradient Boosting Machine
#     3. Support Vector Machine
#=====================================================================================

# Loading all the necessary libraries
# load the libraries
library(tidyverse)
library(caTools)
library(doParallel)
library(mlbench)
library(caret)
library(ggcorrplot)
library(pROC)
library(ROCR)
library(MLmetrics)
library(Hmisc)

# Getting the full set with the maximum n
# Load the Raine study data.
# This data has the following dimensions: 995 (Caucasian) study participants,
# 1511 CpGs, RaineID and maternal smoking variable (preg_SMK_YN: 0: no exposure, 1: exposed)
raine.data <- readRDS("/path/to/raine/data.rds")

#======================================================
# 1. Setting up training and test data sets
#======================================================

# This creates stratified random splits of the data.
set.seed(3456)
trainIndex <- createDataPartition(raine.data$preg_SMK_Y_N, p = .7,
                                  list = FALSE)

# Split the initial data set into training and test data for the first time
# testing data is used in the end to test the model created with cross validation in the training data
training <- raine.data[ trainIndex[,1],]
testing  <- raine.data[-trainIndex[,1],]


# Split into features and class: training set
X_train <- training[,grep("cg",names(training))]
Y_train <- training$preg_SMK_Y_N

#Split into features and class: testing
X_test <- testing[,grep("cg",names(testing))]
Y_test <- testing$preg_SMK_Y_N

#====================================================
# 2. Specify the cross validation settings for all model
#====================================================

#Set the training process and metric
fitControl <- trainControl(method="repeatedcv",
                           number=10,
                           repeats=3,
                           classProbs = TRUE,
                           #Synthetic minority Over-sampling technique: https://jair.org/index.php/jair/article/view/10302
                           sampling="smote")
metric = "Kappa"

#==============================================================
# 3. Modelling
#==============================================================
## 3.1. Random Forest
#==============================================================
# Parallel computing
cl <- makeCluster(detectCores())

# Register the number of parallel workers (here all CPUs)
registerDoParallel(cl)

set.seed(7)
fit.rf <- train(X_train,Y_train,
                method="rf",
                metric=metric,
                trControl=fitControl,
                tuneLength=10)

# return number of parallel workers
getDoParWorkers()
# stop the cluster
stopCluster(cl)

# Save the final model if desired
saveRDS(fit.rf, "fit_rf.rds")


#==============================================================
## 3.2. Gradient Boosting Machine
#==============================================================
# Parallel computing
cl <- makeCluster(detectCores())

# Register the number of parallel workers (here all CPUs)
registerDoParallel(cl)
set.seed(7)
fit.gbm <- train(X_train,Y_train,
                 method="gbm",
                 metric=metric,
                 trControl=fitControl,
                 verbose=FALSE,
                 tuneLength=10)

# return number of parallel workers
getDoParWorkers()
# stop the cluster
stopCluster(cl)

# Save the final model if desired
saveRDS(fit.gbm, "fit_gbm.rds")

#==============================================================
## 3.3. Support Vector Machine
#==============================================================
# Parallel computing
cl <- makeCluster(detectCores())

# Register the number of parallel workers (here all CPUs)
registerDoParallel(cl)

fit.svmLinear <- train(X_train,Y_train,
                       method="svmLinear",
                       metric=metric,
                       preProc=c("center", "scale"),
                       trControl=fitControl,
                       tuneLength=10)



# return number of parallel workers
getDoParWorkers()
# stop the cluster
stopCluster(cl)

# Save the final model if desired
saveRDS(fit.svmLinear, "fit_svmLinear.rds")

#========================================================
# END
#========================================================
