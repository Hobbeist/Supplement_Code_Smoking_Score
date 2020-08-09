#=====================================================================================
# Machine learning-based DNA methylation score for fetal exposure to maternal smoking: 
# development and validation in samples collected from adolescents and adults
#
# Training elastic net regression: The final model
#=====================================================================================

# load libraries
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

# Load the CpG names for the modelling
FinalFeatures <- readRDS("data/FinalFeatures.rds")


# Excluding the highly correlated CpGs
# Find the highly correlated CpGs
# calculate correlation matrix
correlationMatrix <- cor(raine.data[,(grep("cg", names(raine.data)))])

# Set the seed for starting point to ensure the results are repeatable
set.seed(7)

# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)

# Save the names of the highly correlated CpGs
corr.CpGs <- colnames(correlationMatrix[,highlyCorrelated])

set.seed(7)
(CpGnamesfromcorr <- findCorrelation(correlationMatrix, cutoff=0.75, names=T))


# Substract the variables that are suggested to be correlated >.75 
DATA2 <- raine.data[,-which(names(raine.data) %in% CpGnamesfromcorr)]

# Remove data tables to free space in memory
rm(list=c("raine.data"))



# Set up training and test data

# This creates stratified random splits of the data 
set.seed(3456)
trainIndex <- createDataPartition(DATA2$preg_SMK_Y_N, p = .8, 
                                  list = FALSE) 


# Split the initial data set into training and test data for the first time
# testing data is used in the end to test the model created with cross validation in the training data
training <- DATA2[ trainIndex[,1],]
testing  <- DATA2[-trainIndex[,1],]


# create the stratified, minority class sensitive data split
cvIndex <- createMultiFolds(y = training$preg_SMK_Y_N, 
                            k=10, 
                            times=5 )   

# Split into features and class: training set
X_train <- training[,grep("cg",names(training))]
Y_train <- training$preg_SMK_Y_N

# Split into features and class: testing
X_test <- testing[,grep("cg",names(testing))]
Y_test <- testing$preg_SMK_Y_N


# Create the training specifications
my_trainControl <- trainControl(method="repeatedcv",
                                number=10,
                                repeats=5,
                                classProbs = TRUE,
                                #Synthetic minority Over-sampling technique: https://jair.org/index.php/jair/article/view/10302
                                sampling="smote")


# Train glmnet model: using tuneLength

# parallel computing
cl <- makeCluster(detectCores())
# register the number of parallel workers (here all CPUs)
registerDoParallel(cl)

set.seed(7)
# train glmnet model
model_glmnet <- train(x = X_train, y = Y_train, 
                      method = "glmnet",
                      # Specifying the family necessary?
                      # I think it is, as binary outcome
                      family = "binomial",
                      # Kappa is recommended for 
                      # imbalanced data
                      metric = "Kappa",
                      # Make this 100 to explore a large set to 
                      # get pre-selection for tuneGrid
                      tuneLength = 20,
                      preProc=c("center", "scale"),
                      trControl = my_trainControl)

# return number of parallel workers
getDoParWorkers() 
# insert parallel calculations here
# stop the cluster 
stopCluster(cl)

# This shows the model results in a tidy data format
model_glmnet$results %>% 
  as.tibble %>%  
  arrange(desc(Kappa))

# plot ROCs
plot(model_glmnet)

#=====================================================================================
# Apply the model
#=====================================================================================

## 1. Apply to Raine test set
model <- model_glmnet

# Plot the Cross validation steps, apply the score to training set from Raine 
# and print confusion matrix
ggplot(model)
pred.lasso <- predict(model, X_test)
confusionMatrix(pred.lasso, Y_test)

glmnet.score        <- predict(model, X_test, type = "prob")
glmnet.score.class  <- predict(model, X_test)#, type = "prob")


plotting            <- data.frame(glmnet.score=glmnet.score[,2], smoke=Y_test)
plotting.class      <- data.frame(glmnet.score=glmnet.score.class, smoke=Y_test)

# Plot the score stratified by smoke exposure group
ggplot(plotting, aes(x=smoke, y=glmnet.score, fill=smoke))+geom_boxplot()+theme_minimal()

# Barplot for the numbers in the groups based on predicted classification and original
ggplot(plotting.class,aes(smoke, fill=glmnet.score))+geom_bar(position = "dodge")+theme_minimal()

model$bestTune



#=====================================================================================
## 2. Apply to NFBC86
#=====================================================================================
# Load the NFBC86 CpGs
validation.86 <- readRDS("path/to/NFBC86_data.rds")
lasso.validation <- na.omit(validation.86[,c(FinalFeatures,"mat_smk")])

final_predictions <- predict(model, newdata=lasso.validation[,names(X_test)])
confusionMatrix(final_predictions, lasso.validation$mat_smk)



# Relevant plots 
glmnet.score.nfbc <- predict(model, lasso.validation[,names(X_test)], type = "prob")
glmnet.score.class.nfbc <- predict(model, lasso.validation[,names(X_test)])#, type = "prob")

plotting.nfbc <- data.frame(glmnet.score=glmnet.score.nfbc[,2], smoke=lasso.validation$mat_smk)
plotting.class.nfbc <- data.frame(glmnet.score=glmnet.score.class.nfbc, smoke=lasso.validation$mat_smk)

# Boxplot and Barplot
ggplot(plotting.nfbc, aes(x=smoke, y=glmnet.score, fill=smoke))+geom_boxplot()+theme_minimal()
ggplot(plotting.class.nfbc,aes(smoke, fill=glmnet.score))+geom_bar(position = "dodge")+theme_minimal()

#AUC
result.roc <- roc(lasso.validation$mat_smk, predict(model, lasso.validation[,grep("cg",names(lasso.validation))], type="prob")$smoke_exp) # Draw ROC curve.
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")


#=====================================================================================
## 3. Apply to NFBC66
#=====================================================================================

#Load the NFBC66 CpGs
validation.66 <- readRDS("/path/to/NFBC66_data.rds")
lasso.validation <- na.omit(validation.66[,c(names(X_test),"mat_smk")])
final_predictions <- predict(model, newdata=lasso.validation[,grep("cg",names(lasso.validation))])
confusionMatrix(final_predictions, lasso.validation$mat_smk)



# Relevant plots 
glmnet.score.nfbc       <- predict(model, lasso.validation[,grep("cg",names(lasso.validation))], type = "prob")
glmnet.score.class.nfbc <- predict(model, lasso.validation[,grep("cg",names(lasso.validation))])#, type = "prob")

plotting.nfbc       <- data.frame(glmnet.score=glmnet.score.nfbc[,2], smoke=lasso.validation$mat_smk)
plotting.class.nfbc <- data.frame(glmnet.score=glmnet.score.class.nfbc, smoke=lasso.validation$mat_smk)

# Boxplot and Barplot
ggplot(plotting.nfbc, aes(x=smoke, y=glmnet.score, fill=smoke)) + 
  geom_boxplot() + 
  theme_minimal()

ggplot(plotting.class.nfbc,aes(smoke, fill=glmnet.score)) + 
  geom_bar(position = "dodge") + 
  theme_minimal()

#AUC
result.roc <- roc(lasso.validation$mat_smk, predict(model, lasso.validation[,grep("cg",names(lasso.validation))], type="prob")$smoke_exp) # Draw ROC curve.
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")

#========================================================
# END
#========================================================