
## Load mlbench package
#install.packages("klaR")
library(mlbench)
library(dplyr)
library(psych)
library(bestglm)
library(glmnet)
library(nclSLR)
library(MASS)
library(ggplot2)
library(klaR)
## Load the data
data(BreastCancer)
?BreastCancer
## Print first few rows
head(BreastCancer)
BreastCancer=BreastCancer[-1]

#Cleaning the data.
#Cl.thickness,Cell.size,Epith.c.size,Cell.shape (scores),should be numeric with reason for choosing it as a continouus data
#resopnsevariables changed to numeric
BreastCancer$Cl.thickness=as.numeric(BreastCancer$Cl.thickness)
BreastCancer$Cell.size=as.numeric(BreastCancer$Cell.size)
BreastCancer$Cell.shape=as.numeric(BreastCancer$Cell.shape)
BreastCancer$Marg.adhesion=as.numeric(BreastCancer$Marg.adhesion)
BreastCancer$Epith.c.size=as.numeric(BreastCancer$Epith.c.size)
BreastCancer$Bare.nuclei=as.numeric(BreastCancer$Bare.nuclei)
BreastCancer$Bl.cromatin=as.numeric(BreastCancer$Bl.cromatin)
BreastCancer$Normal.nucleoli=as.numeric(BreastCancer$Normal.nucleoli)
BreastCancer$Mitoses=as.numeric(BreastCancer$Mitoses)


#changing class values to either 1 or 0 (predictable variable)

BreastCancer$Class <- as.character(BreastCancer$Class)
BreastCancer$Class[BreastCancer$Class == "benign"] = 0
BreastCancer$Class[BreastCancer$Class == "malignant"] = 1
BreastCancer$Class=as.numeric(BreastCancer$Class)



BreastCancer <- BreastCancer %>% na.omit()
summary(BreastCancer)
count(BreastCancer)


# the relationships between the response variable and predictor variables and about the relationships between predictor variables?
par(mfrow=c(1,1))
hist(BreastCancer$Cl.thickness, main="Frequency plot of Uniformity of Clump Thickness", xlab="Clump Thickness", ylab="Frequency") 
hist(BreastCancer$Cell.size, main="Cell Size", xlab="x", ylab="Frequency") 
hist(BreastCancer$Cell.shape, main=" Frequency plot of Uniformity of Cell Shape", xlab="Uniformity of Cell Shape", ylab="Frequency") 
hist(BreastCancer$Marg.adhesion, main="", xlab="x", ylab="Frequency") 
hist(BreastCancer$Epith.c.size, main="", xlab="x", ylab="Frequency") 
hist(BreastCancer$Bare.nuclei, main=" Frequency plot of Bare Nuclei", xlab="Bare Nuclei", ylab="Frequency") 
hist(BreastCancer$ Normal.nucleoli, main="", xlab="x", ylab="Frequency") 
hist(BreastCancer$Bl.cromatin, main="", xlab="x", ylab="Frequency") 
hist(BreastCancer$Mitoses, main=" Frequency plot of Mitoses", xlab="Mitoses", ylab="Frequency") 
pairs(BreastCancer[],col=BreastCancer[,10]+1)
summary(BreastCancer)
describe(BreastCancer)

BreastCancer_scaled=scale(BreastCancer[1:9])
BreastCancer_scaled_1=scale(BreastCancer[1:10])
BreastCancer_scaled_co = cor(BreastCancer_scaled_1)


# selecting the best variables for 
par(mar=c(10,10,4,14))
barplot(sort(BreastCancer_scaled_co[c(1:9),c(10)]),
        main = "",
        xlab = "Dependency",
        ylab = "",
        las=1,
        names.arg = c("Mitoses","Epith.c.size","Marg.adhesion","Cl.thickness","Normal.nucleoli","Bl.cromatin","Cell.size","Cell.shape","Bare.nuclei"),
        col = "darkred",
        horiz = 1)
# the graph show that Mitoses dependency is relatively low.

############################################################################################################################################################

#pca Analysis
pca2 = prcomp(x=BreastCancer_scaled)
summary(pca2)
pca2$pc1
class(pca2)
#variable reduction table using pca analysis and plotting them 
#plot
plot(pca2$x[,1], pca2$x[,2], xlab="First PC", ylab="Second PC")
# Add labels representing the cancer types 
text(pca2$x[,1], pca2$x[,2], labels=round(BreastCancer[,10]), cex=0.8, pos=3)

legend(x="topleft", pch=1, legend = c("0-benign", "1-malignant"))
#from the plot -1 means the benign and 1 means malignant

##########################################################################################################################################################






#############################data preparation ###############################
X1=scale(BreastCancer[1:9])
y=BreastCancer[,10]
BreastCancer_data=data.frame(X1,y)
n=nrow(BreastCancer_data)
p=ncol(BreastCancer_data)-1
set.seed(123)
cv=trainTestPartition(BreastCancer_data,trainFrac = 0.7)
y=cv$yTr
BreastCancer_data_tr=data.frame(cv$XyTr)
BreastCancer_data_te=data.frame(cv$XyTe)
y_te=cv$yTe
X1=data.frame(cv$XTr)
#############################logistic regression ###############################

logreg_fit_1 = glm(y ~ ., data=data.frame(BreastCancer_data_tr), family="binomial")
summary(logreg_fit_1)

##############subset selection ##################################33

best_fit_aic=bestglm(BreastCancer_data_tr,family=binomial,IC='AIC')
best_fit_bic=bestglm(BreastCancer_data_tr,family=binomial,IC='BIC')
best_fit_aic$Subsets
best_fit_bic$Subsets


## Identify best-fitting models
(best_AIC = best_fit_aic$ModelReport$Bestk)
(best_BIC = best_fit_bic$ModelReport$Bestk)
par(mfrow=c(1,2))
plot(0:p, best_fit_aic$Subsets$AIC, xlab="Number of predictors", ylab="AIC", type="b")
points(best_AIC, best_fit_aic$Subsets$AIC[best_AIC+1], col="red", pch=16)
plot(0:p, best_fit_bic$Subsets$BIC, xlab="Number of predictors", ylab="BIC", type="b")
points(best_BIC, best_fit_bic$Subsets$BIC[best_BIC+1], col="red", pch=16)

##########new model with updated variables##############
pstar = 1
## Check which predictors are in the 1-predictor model
best_fit_aic$Subsets[pstar+5,]
(indices = as.logical(best_fit_aic$Subsets[pstar+5, 2:(p+1)]))
BreastCancer_data_red = data.frame(X1[,indices], y)
## Obtain regression coefficients for this model
logreg1_fit = glm(y ~ ., data=BreastCancer_data_red, family="binomial")
summary(logreg1_fit)

#################Regularisation lasso##########################
set.seed(123)
grid = 10^seq(-1,0, length.out=1000)
## Fit a model with LASSO penalty for each value of the tuning parameter
lass_fit = glmnet(X1, y, family="binomial", alpha=1, standardize=FALSE, lambda=grid)
## Examine the effect of the tuning parameter on the parameter estimates
plot(lass_fit, xvar="lambda", col=rainbow(p), label=TRUE)

lass_cv_fit = cv.glmnet(as.matrix(X1), y, family="binomial", alpha=1, standardize=FALSE, lambda=grid,
                         type.measure="class")
plot(lass_cv_fit)
## Identify the optimal value for the tuning parameter
(lambda_lass_min = lass_cv_fit$lambda.min)

which_lambda_lass = which(lass_cv_fit$lambda == lambda_lass_min)
## Find the parameter estimates associated with optimal value of the tuning parameter
coef(lass_fit, s=lambda_lass_min)


############################################################################

####################Regularization ridge#################################
ridge_fit = glmnet(X1, y, family="binomial", alpha=0, standardize=FALSE, lambda=grid)
## Examine the effect of the tuning parameter on the parameter estimates
plot(ridge_fit, xvar="lambda", col=rainbow(p), label=TRUE)

#par(mfrow=c(1,2))
ridge_cv_fit = cv.glmnet(as.matrix(X1), y, family="binomial", alpha=0, standardize=FALSE, lambda=grid,
                        type.measure="class")
plot(ridge_cv_fit)
## Identify the optimal value for the tuning parameter
(lambda_ridge_min = ridge_cv_fit$lambda.min)

which_lambda_ridge = which(ridge_cv_fit$lambda == lambda_ridge_min)
## Find the parameter estimates associated with optimal value of the tuning parameter
coef(ridge_fit, s=lambda_ridge_min)


######################LDA ANALYSIS #############################################3
#k=linDA(variables=data.frame(X1),group=y)
#k
model <- lda(BreastCancer_data_tr$y~., data = data.frame(BreastCancer_data_tr),CV=FALSE)
model
#######################Qda Analysis##############################################
#k=linDA(variables=data.frame(X1),group=y)
#k

model_q <- qda(BreastCancer_data_tr$y~., data = data.frame(BreastCancer_data_tr),CV=FALSE)
model_q
summary(model_q)

############k fold cv based on test error##########################################
#model selection algorith
m=matrix(data=NA,nrow=10,ncol=5)






for (x in 1:10) 
{
 
  
  cv_2=trainTestPartition(BreastCancer_data_te,trainFrac = 0.7)
  grid_1 = 10^seq(-1,0, length.out=500)
  
  ###############lda#####################
  model_lda=lda(cv_2$yTr~., data = data.frame(cv_2$XTr),CV=FALSE)
  model_qda=qda(cv_2$yTr~., data = data.frame(cv_2$XTr))
  model_lasso=glmnet(data.frame(cv_2$XTr),cv_2$yTr, family="binomial", alpha=1, standardize=FALSE, lambda=grid_1)
  lass_cv_fit_1 = cv.glmnet(as.matrix(cv_2$XTr), cv_2$yTr, family="binomial", alpha=1, standardize=FALSE, lambda=grid_1,
                          type.measure="class")
  lambda_lass_min_1 = lass_cv_fit_1$lambda.min
  
  ridge_fit_1 = glmnet(as.matrix(cv_2$XTr), cv_2$yTr, family="binomial", alpha=0, standardize=FALSE, lambda=grid_1)
  ridge_cv_fit_1 = cv.glmnet(as.matrix(cv_2$XTr), cv_2$yTr, family="binomial", alpha=0, standardize=FALSE, lambda=grid_1,
                           type.measure="class")
  lambda_ridge_min_1 = ridge_cv_fit_1$lambda.min
  
  reg_red_mod = glm(cv_2$yTr ~ Cl.thickness + Epith.c.size + Bare.nuclei + Bl.cromatin + Normal.nucleoli, data=data.frame(cv_2$XyTr), family="binomial")
  
  
  
  lda_test=predict(model_lda,cv_2$XyTe)
  qda_test=predict(model_qda,cv_2$XyTe)
  lasso_test = predict(model_lasso, cv_2$XTe, s=lambda_lass_min_1, type="response")
  ridge_test = predict(ridge_fit_1, cv_2$XTe, s=lambda_ridge_min_1, type="response")
  reg_red = predict(reg_red_mod, cv_2$XyTe, type="response")
  
  yhat_test_1=lda_test$class
  yhat_test_2=qda_test$class
  yhat_test_3 = ifelse(lasso_test > 0.5, 1, 0)
  yhat_test_4 = ifelse(ridge_test > 0.5, 1, 0)
  yhat_test_5 = ifelse(reg_red > 0.5, 1, 0)
  
  confusion_lda=table(observed=cv_2$yTe,predicted=yhat_test_1)
  confusion_qda=table(observed=cv_2$yTe,predicted=yhat_test_2)
  confusion_lasso=table(observed=cv_2$yTe,predicted=yhat_test_3)
  confusion_ridge=table(observed=cv_2$yTe,predicted=yhat_test_4)
  confusion_reg_red = table(Observed=cv_2$yTe, Predicted=yhat_test_5)
 
  m[x,1]= 1-mean(cv_2$yTe == yhat_test_1)
  m[x,2]= 1-mean(cv_2$yTe == yhat_test_2)
  m[x,3]= 1-mean(cv_2$yTe == yhat_test_3)
  m[x,4]= 1-mean(cv_2$yTe == yhat_test_4)
  m[x,5]= 1-mean(cv_2$yTe == yhat_test_5)
  
}
m_dt_fr=data.frame(m)
#renaming 
names(m_dt_fr)[names(m_dt_fr) == 'X1'] <- 'model_lda'
names(m_dt_fr)[names(m_dt_fr) == 'X2'] <- 'model_qda'
names(m_dt_fr)[names(m_dt_fr) == 'X3'] <- 'model_lasso'
names(m_dt_fr)[names(m_dt_fr) == 'X4'] <- 'model_ridge'
names(m_dt_fr)[names(m_dt_fr) == 'X5'] <- 'model_log_red_reg'
m_dt_fr
describe(m_dt_fr)
summary(m_dt_fr)

