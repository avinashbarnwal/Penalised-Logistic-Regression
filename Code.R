require(RCurl)
require(prettyR)
require(caret)
require(e1071)
library(ROCR)
setwd("/Users/avinashbarnwal/Desktop/Spring 2017/Categorical Data Analysis/Project")


url <- "https://raw.githubusercontent.com/gastonstat/CreditScoring/master/CleanCreditScoring.csv"

cs_data <- getURL(url)
cs_data <- read.csv(textConnection(cs_data))
str(cs_data)
describe(cs_data)
write.table(cs_data,file="credit_scoring.csv",sep=",")

cs_data<-cs_data[, -match(c( "seniorityR", "timeR", "ageR","expensesR", "incomeR", "assetsR", "debtR", "amountR", "priceR", "finratR","savingsR"), colnames(cs_data))]
										

classes <- cs_data[, "Status"]

train_set <- createDataPartition(classes, p = 0.8, list = FALSE)
cs_data_train <- cs_data[train_set, ]
cs_data_test <- cs_data[-train_set, ]

set.seed(124)
cv_splits <- createFolds(classes, k = 10, returnTrain = TRUE)
glmnet_grid <- expand.grid(alpha = c(0,  .1,  .2, .4, .6, .8, 1),
                           lambda = seq(.01, .2, length = 20))
glmnet_ctrl <- trainControl(method = "cv", number = 10)
glmnet_fit <- train(Status ~ ., data = cs_data_train,
                    method = "glmnet",preProcess = c("center", "scale"),
                    tuneGrid = glmnet_grid,trControl = glmnet_ctrl)
glmnet_fit

trellis.par.set(caretTheme())
plot(glmnet_fit, scales = list(x = list(log = 2)))
pred_penalised <- predict(glmnet_fit, newdata = cs_data_test)

table(cs_data_test$Status,pred_penalised)

pr <- prediction(pred_penalised, cs_data_test$Status)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)


##Ordinary Logistic Regression
model <- glm(Status ~.,family=binomial(link='logit'),data=cs_data_train)
summary(model)
anova(model, test="Chisq")


pred_ordinary <- predict(model, newdata = cs_data_test,type='response')
fitted.results <- ifelse(pred_ordinary > 0.5,1,0)
table(cs_data_test$Status,fitted.results)

pr <- prediction(pred_ordinary, cs_data_test$Status)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)




prostate <- read.delim("http://web.as.uky.edu/statistics/users/pbreheny/603/prostate.txt")
