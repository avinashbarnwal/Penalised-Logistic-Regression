library(RCurl)
library(prettyR)
library(caret)
library(e1071)
library(ROCR)
library(nnet)
library(pROC)
setwd("/Users/avinashbarnwal/Desktop/Spring 2017/Categorical Data Analysis/Project")

url <- "https://raw.githubusercontent.com/avinashbarnwal/Penalised-Logistic-Regression/master/yeast_dataset.csv"
yeast_data <- getURL(url)
yeast_data <- read.csv(textConnection(yeast_data))
yeast_data<-yeast_data[-1]
classes <- yeast_data[, "Localization"]
no_var=ncol(yeast_data)-1

################################# TEST AND TRAIN #####################################
train_set <- createDataPartition(classes, p = 0.8, list = FALSE)
yeast_data_train <- yeast_data[train_set, ]
yeast_data_test  <- yeast_data[-train_set, ]

##################################FIRST ORDER#########################################

              ###############ORDINARY LOGISTIC REGRESSION################
first_order_predictors=colnames(yeast_data)[1:no_var]
first_order_com=paste(first_order_predictors,collapse="+")
null_formula=paste("Localization~","1",sep="")
full_formula=paste("Localization~",first_order_com,sep="")

fullmod <- multinom(full_formula, data=yeast_data_train,MaxNWts=10000)
nothing <- multinom(null_formula,data=yeast_data_train,MaxNWts=10000)
forwards = step(nothing,scope=list(lower=formula(nothing),upper=formula(fullmod)),data=yeast_data_train, direction="forward",
                k=log(nrow(yeast_data_train)))
summary(forwards)
coefficients_name     = colnames(coefficients(forwards))[-1]
first_order_var       = coefficients_name
coefficients_name_com = paste(coefficients_name,collapse="+")
first_formula         = paste("Localization~",coefficients_name_com,sep="")
final_model           = multinom(first_formula, data=yeast_data_train,MaxNWts=10000)

              ##################PENALISED REEGRESSION############
set.seed(124)
glmnet_grid <- expand.grid(alpha = c(0,  .1,  .2, .4, .6, .8, 1),lambda = seq(.01, 1, length = 20))
glmnet_ctrl <- trainControl(method = "cv", number = 10)
glmnet_fit <- train(Localization ~ ., data = yeast_data_train,method = "glmnet",preProcess = c("center", "scale"),
                    tuneGrid = glmnet_grid,trControl = glmnet_ctrl)

trellis.par.set(caretTheme())
plot(glmnet_fit, scales = list(x = list(log = 2)))

#The final values used for the model were alpha = 1 and lambda = 0.01. 
               ##################COMPARISON#####################

###PENALISED REGRESSION###      

pred_penalised <- predict(glmnet_fit, newdata = yeast_data_test)
pred_penalised_table<-table(yeast_data_test$Localization,pred_penalised)
total=sum(apply(pred_penalised_table,1,function(x){sum(x)}))
accuracy =sum(diag(pred_penalised_table))
accuracy_perc=accuracy/total
# 0.5290102

###ORDINARY REGRESSION###

pred_ordinary = predict(final_model, yeast_data_test, type = "class")
pred_ordinary_table = table(yeast_data_test$Localization,pred_ordinary)
total=sum(apply(pred_ordinary_table,1,function(x){sum(x)}))
accuracy =sum(diag(pred_ordinary_table))
accuracy_perc=accuracy/total
# 0.5767918


##################################SECOND ORDER#########################################

##################################VARIABLE CREATION####################################

set.seed(124)
yeast_data2=yeast_data

var_name2=vector()
flag=0
for(i in 1:no_var){
  for(j in i:no_var){
    flag=flag+1
    var_name2[flag]=paste(colnames(yeast_data2)[i],"_",colnames(yeast_data2)[j],sep="")
    yeast_data2[,var_name2[flag]]=yeast_data2[,i]*yeast_data2[,j]
  }
}


train_set <- createDataPartition(classes, p = 0.8, list = FALSE)
yeast_data2_train <- yeast_data2[train_set, ]
yeast_data2_test  <- yeast_data2[-train_set, ]


###############ORDINARY LOGISTIC REGRESSION################

second_order=append(first_order_var,var_name2)
second_order_com=paste(second_order,collapse="+")
null_formula=first_formula
full_formula=paste("Localization~",second_order_com,sep="")

fullmod <- multinom(full_formula, data=yeast_data2_train,MaxNWts=10000)
nothing <- multinom(null_formula,data=yeast_data2_train,MaxNWts=10000)

forwards = step(nothing,scope=list(lower=formula(nothing),upper=formula(fullmod)),data=yeast_data2_train, direction="forward",
                k=log(nrow(yeast_data2_train)))
summary(forwards)
coefficients_name     = colnames(coefficients(forwards))[-1]
second_order_var      = coefficients_name
coefficients_name_com = paste(coefficients_name,collapse="+")
second_formula         = paste("Localization~",coefficients_name_com,sep="")
final_model           = multinom(second_formula, data=yeast_data2_train,MaxNWts=10000)


##################PENALISED REEGRESSION############
set.seed(124)
glmnet_grid <- expand.grid(alpha = c(0,  .1,  .2, .4, .6, .8, 1),lambda = seq(.01, 2, length = 40))
glmnet_ctrl <- trainControl(method = "cv", number = 10)
glmnet_fit <- train(Localization ~ ., data = yeast_data2_train,method = "glmnet",preProcess = c("center", "scale"),
                    tuneGrid = glmnet_grid,trControl = glmnet_ctrl)

trellis.par.set(caretTheme())
plot(glmnet_fit, scales = list(x = list(log = 2)))
#The final values used for the model were alpha = 0.4 and
#lambda = 0.01.
##################COMPARISON#####################

###PENALISED REGRESSION###      

pred_penalised <- predict(glmnet_fit, newdata = yeast_data2_test)
pred_penalised_table<-table(yeast_data2_test$Localization,pred_penalised)
total=sum(apply(pred_penalised_table,1,function(x){sum(x)}))
accuracy =sum(diag(pred_penalised_table))
accuracy_perc=accuracy/total
# 0.6109215

###ORDINARY REGRESSION###

pred_ordinary = predict(final_model, yeast_data2_test, type = "class")
pred_ordinary_table = table(yeast_data2_test$Localization,pred_ordinary)
total=sum(apply(pred_ordinary_table,1,function(x){sum(x)}))
accuracy =sum(diag(pred_ordinary_table))
accuracy_perc=accuracy/total
# 0.6040956



##################################THIRD ORDER#########################################

##################################VARIABLE CREATION####################################

var_name3=vector()
flag=0
yeast_data3=yeast_data2
for(i in 1:no_var){
  for(j in i:no_var){
    for(k in j:no_var){
      flag=flag+1
      var_name3[flag]=paste(colnames(yeast_data3)[i],"_",colnames(yeast_data3)[j],"_",colnames(yeast_data3)[k],sep="")
      yeast_data3[,var_name3[flag]]=yeast_data3[,i]*yeast_data3[,j]*yeast_data3[,k]
    }
  }
}

train_set         <- createDataPartition(classes, p = 0.8, list = FALSE)
yeast_data3_train <- yeast_data3[train_set, ]
yeast_data3_test  <- yeast_data3[-train_set, ]



###############ORDINARY LOGISTIC REGRESSION################

third_order=append(second_order_var,var_name3)
third_order_com=paste(third_order,collapse="+")
null_formula=second_formula
full_formula=paste("Localization~",third_order_com,sep="")

fullmod <- multinom(full_formula, data=yeast_data3_train,MaxNWts=10000)
nothing <- multinom(null_formula, data=yeast_data3_train,MaxNWts=10000)

forwards = step(nothing,scope=list(lower=formula(nothing),upper=formula(fullmod)),data=yeast_data3_train, direction="forward",
                k=log(nrow(yeast_data3_train)))
summary(forwards)
coefficients_name     = colnames(coefficients(forwards))[-1]
third_order_var       = coefficients_name
coefficients_name_com = paste(coefficients_name,collapse="+")
third_formula         = paste("Localization~",coefficients_name_com,sep="")
final_model           = multinom(third_formula, data=yeast_data3_train,MaxNWts=10000)


##################PENALISED REEGRESSION############
set.seed(124)
glmnet_grid <- expand.grid(alpha = c(0,  .1,  .2, .4, .6, .8, 1),lambda = seq(.01, 2, length = 40))
glmnet_ctrl <- trainControl(method = "cv", number = 10)
glmnet_fit <- train(Localization ~ ., data = yeast_data3_train,method = "glmnet",preProcess = c("center", "scale"),
                    tuneGrid = glmnet_grid,trControl = glmnet_ctrl)

trellis.par.set(caretTheme())
plot(glmnet_fit, scales = list(x = list(log = 2)))
#The final values used for the model were alpha = 0.1 and
#lambda = 0.01. 
##################COMPARISON#####################

###PENALISED REGRESSION###      

pred_penalised <- predict(glmnet_fit, newdata = yeast_data3_test)
pred_penalised_table<-table(yeast_data3_test$Localization,pred_penalised)
total=sum(apply(pred_penalised_table,1,function(x){sum(x)}))
accuracy =sum(diag(pred_penalised_table))
accuracy_perc=accuracy/total
# 0.5631399

###ORDINARY REGRESSION###

pred_ordinary = predict(final_model, yeast_data3_test, type = "class")
pred_ordinary_table = table(yeast_data3_test$Localization,pred_ordinary)
total=sum(apply(pred_ordinary_table,1,function(x){sum(x)}))
accuracy =sum(diag(pred_ordinary_table))
accuracy_perc=accuracy/total
# 0.5460751

##################################FOURTH ORDER#########################################

##################################VARIABLE CREATION####################################

var_name4=vector()
flag=0
yeast_data4=yeast_data3
for(i in 1:no_var){
  for(j in i:no_var){
    for(k in j:no_var){
      for(l in k:no_var){
        flag=flag+1
        var_name4[flag]=paste(colnames(yeast_data4)[i],"_",colnames(yeast_data4)[j],"_",colnames(yeast_data4)[k],"_",colnames(yeast_data4)[l],sep="")
        yeast_data4[,var_name4[flag]]=yeast_data4[,i]*yeast_data4[,j]*yeast_data4[,k]*yeast_data4[,l]
      }
    }
  }
}

set.seed(124)
train_set         <- createDataPartition(classes, p = 0.8, list = FALSE)
yeast_data4_train <- yeast_data4[train_set, ]
yeast_data4_test  <- yeast_data4[-train_set, ]

###############ORDINARY LOGISTIC REGRESSION################

fourth_order=append(third_order_var,var_name4)
fourth_order_com=paste(fourth_order,collapse="+")
null_formula=third_formula
full_formula=paste("Localization~",fourth_order_com,sep="")

fullmod <- multinom(full_formula, data=yeast_data4_train,MaxNWts=10000)
nothing <- multinom(null_formula, data=yeast_data4_train,MaxNWts=10000)

forwards = step(nothing,scope=list(lower=formula(nothing),upper=formula(fullmod)),data=yeast_data4_train, direction="forward",
                k=log(nrow(yeast_data4_train)))


summary(forwards)
coefficients_name     = colnames(coefficients(forwards))[-1]
fourth_order_var      = coefficients_name
coefficients_name_com = paste(coefficients_name,collapse="+")
fourth_formula        = paste("Localization~",coefficients_name_com,sep="")
final_model           = multinom(third_formula, data=yeast_data4_train,MaxNWts=10000)


##################PENALISED REEGRESSION############
set.seed(124)
glmnet_grid <- expand.grid(alpha = c(0,  .1,  .2, .4, .6, .8, 1),lambda = seq(.01, 1, length = 20))
glmnet_ctrl <- trainControl(method = "cv", number = 10)
glmnet_fit <- train(Localization ~ ., data = yeast_data4_train,method = "glmnet",preProcess = c("center", "scale"),
                    tuneGrid = glmnet_grid,trControl = glmnet_ctrl)

trellis.par.set(caretTheme())
plot(glmnet_fit, scales = list(x = list(log = 2)))
#The final values used for the model were alpha = 0.1 and
#lambda = 0.01. 
##################COMPARISON#####################

###PENALISED REGRESSION###      

pred_penalised <- predict(glmnet_fit, newdata = yeast_data4_test)
pred_penalised_table<-table(yeast_data4_test$Localization,pred_penalised)
total=sum(apply(pred_penalised_table,1,function(x){sum(x)}))
accuracy =sum(diag(pred_penalised_table))
accuracy_perc=accuracy/total
# 0.6006826
#alpha = 0.4 and
#lambda = 0.01. 

###ORDINARY REGRESSION###

pred_ordinary = predict(final_model, yeast_data4_test, type = "class")
pred_ordinary_table = table(yeast_data4_test$Localization,pred_ordinary)
total=sum(apply(pred_ordinary_table,1,function(x){sum(x)}))
accuracy =sum(diag(pred_ordinary_table))
accuracy_perc=accuracy/total
#0.5836




