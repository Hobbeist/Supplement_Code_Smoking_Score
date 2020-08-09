#=====================================================================================
# Predict with Random Forest, Gradient Boosting Machine and Support vector machine
#
#=====================================================================================

#==============================================================
# 1. Load libraries
library(tidyverse)
library(caret)

#==============================================================
# 2.Loading your data
# This requires to have exactly the 1511 CpGs available that are
# listed in Excel Table S1! The CpGs need to be columns, the 
# study participants need to be rows
#==============================================================

X_test <- read_csv("Your/data/file.csv")

#==============================================================
# 3. Load the models saved in the training step
#==============================================================
# 3.1 Random Forest
rf.fit <- readRDS("models/randomForest_model.rds")

# These two lines create the smoking score based on random forest and a class wise prediction (smoke exposed versus not exposed)
glmnet.score        <- predict(rf.fit, X_test, type = "prob")
glmnet.score.class  <- predict(rf.fit, X_test)

#==============================================================
# 3.2 Gradient Boosting Machine
gbm.fit <- readRDS("models/gradientBoosting_model.rds")

# These two lines create the smoking score based on random forest and a class wise prediction (smoke exposed versus not exposed)
glmnet.score        <- predict(gbm.fit, X_test, type = "prob")
glmnet.score.class  <- predict(gbm.fit, X_test)

#==============================================================
# 3.1 Support Vector Machine
svm.fit <- readRDS("models/supportVector_model.rds")

# These two lines create the smoking score based on random forest and a class wise prediction (smoke exposed versus not exposed)
glmnet.score        <- predict(svm.fit, X_test, type = "prob")
glmnet.score.class  <- predict(svm.fit, X_test)


#========================================================
# END
#========================================================