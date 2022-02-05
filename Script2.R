# Author: Julian Sauvage
# Date: 2/4/2022
# Class: Econ 613
# Assignment 2: OLS and Probit

setwd('/Users/juliansauvage/Desktop/Econ613/Assignment 2/Data/')
library(tidyverse)

ind2009 <- read.csv("datind2009.csv")

df = data.frame(NULL)
df<- ind2009 %>%
  subset(select = c(age, wage)) %>%
  drop_na() 

#Calculate correlation between Y (wage) and X (age)
df_mean <- matrix(data=1, nrow=nrow(df)) %*% cbind(mean(df$age),mean(df$wage))
diff_matrix <- df - df_mean
covar_matrix <- length(df)^-1 * t(diff_matrix) %*% as.matrix(diff_matrix)
square_root <- diag(diag(covar_matrix)^(-1/2))
square_root %*% covar_matrix %*% square_root

#check against built-in correlation function
cor(df)

#Calculate regression coefficient Beta
ind2009 <- ind2009 %>%
  drop_na()
#Separate independent variable from dependent, and append a column representing the intercept
df<- ind2009 %>%
  subset(select = c(age)) %>%
  cbind(1, .)

#Solve for coefficients using ((X^T*X)^-1)*(X^T*Y)
beta <- solve(t(df)%*%as.matrix(df)) %*% (t(df)%*%as.matrix(ind2009$wage))

#Check coefficients against built-in linear regression function
lm1 <- lm(wage ~ age, ind2009)
summary(lm1)

#Calculate Standard Errors of Reg. Coefficients
#Begin by estimating wages and residuals
wagehat = as.matrix(df) %*% as.matrix(beta)
epsHat = wagehat - ind2009$wage

#Calculate variance/covariance matrix
s = t(epsHat) %*% epsHat / (nrow(ind2009)-ncol(ind2009)-1)
v = s[1,1] * solve(t(df)%*% as.matrix(df))
#Standard errors are the square root of diagonal elements
se = diag(v)^(1/2)

standard_est = cbind(beta, se)

# Calculate coefficients/standard error using bootstrap

R    = 49                      # number of bootstraps, round 1, 
R2   = 499                     # number of bootstraps, round 2
nind = nrow(ind2009)           # number of individuals
nvar = length(df)   # number of variables

outs = mat.or.vec(R,nvar)
set.seed(123)

for (i in 1:R)
{
  samp     = sample(1:nind,nind,rep=TRUE)
  dat_samp = ind2009[samp,]
  reg1     = lm(wage ~ age,data = dat_samp)
  outs[i,] = reg1$coefficients
}

mean_est = apply(outs,2,mean)
se_est   = apply(outs,2,sd)

est = cbind(mean_est,
            se_est)

outs2 = mat.or.vec(R2,nvar)
set.seed(123)

for (i in 1:R2)
{
  samp     = sample(1:nind,nind,rep=TRUE)
  dat_samp = ind2009[samp,]
  reg1     = lm(wage ~ age,data = dat_samp)
  outs2[i,] = reg1$coefficients
}

mean_est = apply(outs2,2,mean)
se_est   = apply(outs2,2,sd)

est = cbind(standard_est, est, mean_est,se_est)

colnames(est) = c("Manual: est", "Manual: sd", "Bootstrap49: est","Bootstrap49: se", "Bootstrap499: est","Bootstrap499: se")
est

##############
# Exercise 2 #
##############
setwd("/Users/juliansauvage/Desktop/Econ613/Assignment 1/Data/")

allind <- data.frame(NULL)

# years is the sequence to iterate over
years <- 2004:2019

# Read in all datasets, and append all datasets
# Creates allhh and allind to store household/individual data respectively
for (i in years) {
  # temporary name variables represent filenames to be read in
  nameind <- paste0("datind", i,".csv")

  # Read in individual data, and append to existing ind df
  tmpdat <- read.csv(nameind, header=T, stringsAsFactors=F)
  allind <- rbind(allind, tmpdat)
}



allind$ag <- as.factor(ifelse(allind$age<26, '18-25', 
                                 ifelse(allind$age<31, '26-30',
                                        ifelse(allind$age<36, '31-35', 
                                               ifelse(allind$age<41, '36-40', 
                                                      ifelse(allind$age<46, '41-45',
                                                             ifelse(allind$age<51, '46-50',
                                                                    ifelse(allind$age<56, '51-55',
                                                                           ifelse(allind$age<61, '56-60','60+')))))))))

df <- allind %>%
  mutate(year = factor(year)) %>%
  filter(!is.na(ag)) %>%
  filter(!is.na(wage)) %>%
  filter(wage != 0)


ggplot(data=df %>% group_by(year, ag) %>% summarize(mean_income = mean(wage)), aes(x=year, y=mean_income, group=ag, colour=ag)) +
  geom_line() + labs(title = 'Average Wages by Age Group (2004-2019)', x = 'Year', y = 'Wages', color = 'Age Group')

ggplot(data=df  %>% group_by(year, ag) %>% mutate(mean_income = mean(wage)),aes(x=year, y=wage, fill=ag)) +
  geom_boxplot(outlier.shape = NA) + labs(title = 'Wage Distribution by Age Group (2004-2019)', x = 'Year', y = 'Wages', fill = 'Age Group', color = 'Average Wage') +
  geom_line(aes(x=year, y=mean_income, group=ag, colour=ag)) +
  coord_cartesian(ylim = quantile(df$wage, c(0.1, 0.99)))

lm2 <- lm(wage ~ age + (factor(year)), allind)
summary(lm2)


##############
# Exercise 3 #
##############
rm(list=ls())
ind2007 <- read.csv("datind2007.csv")



ind2007 <- ind2007 %>%
  filter(empstat != 'Retired' & empstat != 'Inactive') %>%
  mutate(employed = ifelse(empstat == 'Employed', 1, 0)) %>%
  drop_na(employed)

df <- ind2007 %>%
  subset(select = c(age)) %>%
  drop_na() %>%
  cbind(1, .)

beta <- solve(t(df)%*%as.matrix(df)) %*% (t(df)%*%as.matrix(ind2007$employed))

problike = function(par,x1,yvar)
{
  xbeta           = par[1] + par[2]*x1
  pr              =  pnorm(xbeta) 
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like           = yvar*log(pr) + (1-yvar)*log(1-pr)
  return(-sum(like))
}

problike(beta, ind2007$age, ind2007$employed)

prob_model <- glm(employed ~ age, family = binomial(link = "probit"), 
                data = ind2007)
logLik(prob_model)

res  = optim(beta,fn=problike,method="BFGS",control=list(trace=6,REPORT=1,maxit=1000),x1=ind2007$age,yvar=ind2007$employed,hessian=TRUE)
fisher_info = solve(res$hessian)       
prop_sigma  = sqrt(diag(fisher_info))
prop_sigma

est = cbind(summary(prob_model)$coefficients[, 1],summary(prob_model)$coefficients[, 2],res$par,prop_sigma)
colnames(est) = c("Built-in Probit: est","Built-in Probit :se","Own Probit : est","Own Probit :se")
est

##############
# Exercise 4 #
##############
rm(list=ls())

allind <- data.frame(NULL)

# years is the sequence to iterate over
years <- 2005:2015

# Read in all datasets, and append all datasets
# Creates allhh and allind to store household/individual data respectively
for (i in years) {
  # temporary name variables represent filenames to be read in
  nameind <- paste0("datind", i,".csv")
  
  # Read in individual data, and append to existing ind df
  tmpdat <- read.csv(nameind, header=T, stringsAsFactors=F)
  allind <- rbind(allind, tmpdat)
}

allactive <- allind %>%
  filter(empstat != 'Retired' & empstat != 'Inactive') %>%
  mutate(employed = ifelse(empstat == 'Employed', 1, 0))

df <- allactive %>%
  subset(select = c(age)) %>%
  drop_na() %>%
  cbind(1, .)

beta <- solve(t(df)%*%as.matrix(df)) %*% (t(df)%*%as.matrix(allactive$employed))

#Probit Model
problike = function(par,x1,yvar)
{
  xbeta           = par[1] + par[2]*x1
  pr              =  pnorm(xbeta) 
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like           = yvar*log(pr) + (1-yvar)*log(1-pr)
  return(-sum(like))
}

#Logit Model
loglike = function(par,x1,yvar)
{
  xbeta           = par[1] + par[2]*x1
  pr              =  exp(xbeta)/(1+exp(xbeta)) 
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like           = yvar*log(pr) + (1-yvar)*log(1-pr)
  return(-sum(like))
}

#Linear Prob model
linearlike = function(par,x1,yvar)
{
  xbeta           = par[1] + par[2]*x1
  pr              =  xbeta 
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like           = yvar*log(pr) + (1-yvar)*log(1-pr)
  return(-sum(like))
}

# Optimize probit model and compare with estimates from glm()
prob_fun = problike(beta, allactive$age, allactive$employed)

res1  = optim(beta,fn=problike,method="BFGS",control=list(trace=6,REPORT=1,maxit=1000),x1=allactive$age,yvar=allactive$employed,hessian=TRUE)
fisher_info1 = solve(res1$hessian)       
prop_sigma1  = sqrt(diag(fisher_info1))
prop_sigma1

prob_model <- glm(employed ~ age, family = binomial(link = "probit"), 
                  data = allactive)

est1 = cbind(summary(prob_model)$coefficients[, 1],summary(prob_model)$coefficients[, 2],res1$par,prop_sigma1)
colnames(est1) = c("Built-in Probit: est","Built-in Probit :se","Own Probit : est","Own Probit :se")
est1

#Optimize logit model and compare with estimates from glm() 
log_fun = loglike(beta, allactive$age, allactive$employed)

res2  = optim(beta,fn=loglike,method="BFGS",control=list(trace=6,REPORT=1,maxit=1000),x1=allactive$age,yvar=allactive$employed,hessian=TRUE)
fisher_info2 = solve(res2$hessian)       
prop_sigma2  = sqrt(diag(fisher_info2))
prop_sigma2

log_model <- glm(employed ~ age, family = binomial(link = "logit"), 
                  data = allactive)

est2 = cbind(summary(log_model)$coefficients[, 1],summary(log_model)$coefficients[, 2],res2$par,prop_sigma2)
colnames(est2) = c("Built-in Logit: est","Built-in Logit :se","Own Logit : est","Own Logit :se")
est2

#Optimize linear model and compare with estimates from glm()
linearlike(beta, allactive$age, allactive$employed)

res3  = optim(beta,fn=linearlike,method="BFGS",control=list(trace=6,REPORT=1,maxit=1000),x1=allactive$age,yvar=allactive$employed,hessian=TRUE)
fisher_info3 = solve(res3$hessian)       
prop_sigma3  = sqrt(-diag(fisher_info3))
prop_sigma3

lin_model <- glm(employed ~ age, data = allactive)
est3 = cbind(summary(lin_model)$coefficients[, 1],summary(lin_model)$coefficients[, 2],res3$par,prop_sigma3)
colnames(est3) = c("Built-in Linear: est","Built-in Linear :se","Own Linear : est","Own Linear :se")
est3

##############
# Exercise 5 #
##############

# Compute marginal effect of age from probit model
probmarg = function(par,x1,yvar)
{
  xbeta           = par[1] + par[2]*x1
  pr              =  pnorm(xbeta)
  marg            = t(as.matrix(pr))*par[2]
  return(sum(marg)/length(marg))
}

p_margin = probmarg(beta, allactive$age, allactive$employed)

# Compute marginal effect of age from logit model
logmarg = function(par,x1,yvar)
{
  xbeta           = par[1] + par[2]*x1
  pr              =  exp(xbeta)/(1+exp(xbeta)) 
  marg            = t(as.matrix(pr))*par[2]
  return(sum(marg)/length(marg))
}

l_margin = logmarg(beta, allactive$age, allactive$employed)



