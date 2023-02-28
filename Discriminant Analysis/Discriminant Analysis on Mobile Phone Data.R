## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
train = read.csv("train.csv")
knitr::kable(head(train[,1:12],2))
knitr::kable(head(train[,12:21],2))


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
sum(is.na(train))


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# https://stackoverflow.com/questions/17171148/non-redundant-version-of-expand-grid
# used to get the combinations of plots to look at since pairs() is overwhelmed by the 210 plots
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
    x <- unique(x)

    y <- unique(y)

    g <- function(i)
    {
        z <- setdiff(y, x[seq_len(i-include.equals)])

        if(length(z)) cbind(x[i], z, deparse.level=0)
    }

    do.call(rbind, lapply(seq_along(x), g))
}

combinations = expand.grid.unique(1:21,1:21)
# WARNING: 210 graphs will be printed if you uncomment and run for loop below
#for(i in 1:210){
#  plot(train[,combinations[i,1]],train[,combinations[i,2]],xlab = colnames(train)[combinations[i,1]],ylab=colnames(train)[combinations[i,2]])
#}


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#for (i in 1:20){
#  plot(train[,i],train[,21],xlab=colnames(train)[i],ylab = colnames(train)[21])
#}


## ----echo=F, fig.align='center', fig.cap="Histograms of All Variables", warning=FALSE, out.width="70%"----------------------------------------------------------------------------
par(mfrow=c(2,4))
for (i in 1:8){
  hist(train[,i],main=colnames(train)[i],xlab=colnames(train)[i])
}


## ----echo=F, fig.align='center', fig.cap="Histograms of All Variables", warning=FALSE, out.width="70%"----------------------------------------------------------------------------
par(mfrow=c(2,4))
for (i in 9:16){
  hist(train[,i],main=colnames(train)[i],xlab=colnames(train)[i])
}


## ----echo=F, fig.align='center', fig.cap="Histograms of All Variables", warning=FALSE, out.width="70%"----------------------------------------------------------------------------
par(mfrow=c(2,4))
for (i in 17:21){
  hist(train[,i],main=colnames(train)[i],xlab=colnames(train)[i])
}


## ----include=F, echo=F------------------------------------------------------------------------------------------------------------------------------------------------------------
#train %>% summarise_if(is.numeric, min)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
add_one = c('blue','dual_sim','fc','four_g','pc','px_height','sc_w','three_g','touch_screen','wifi','price_range')
for(i in add_one){
  train[,i] = train[,i] + 10
}


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
train[,'fc'] = log(train[,'fc'])
train[,'m_dep'] = log(train[,'m_dep'])
train[,'clock_speed'] = log(train[,'clock_speed'])
train[,'n_cores'] = log(train[,'n_cores'])
train[,'px_height'] = sqrt(train[,'px_height'])
train[,'sc_w'] = log(train[,'sc_w'])
train[,'sc_h'] = log(train[,'sc_h'])


## ----echo=F, fig.align='center', fig.cap=" Transformation of Variables", warning=FALSE, out.width="70%"---------------------------------------------------------------------------
knitr::opts_chunk$set(fig.pos = 'H')
par(mfrow=c(2,4))
columns_transformed = c('px_height','sc_w','clock_speed','fc','m_dep','n_cores','sc_h')
for (i in columns_transformed){
  hist(train[,i],main=i,xlab=i)
}


## ----include=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------
library(MASS)
flda=function(x,class){
cat("Fisher Linear Discriminant:\n")
a = lda(x, class); d = predict(a)
t=table(class, d$class);
er=100*(sum(t)-sum(diag(t)))/nrow(x)
#cat("Error Rate =",er,"%\n")
return(list(t=t,er=er,d=d))
}


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
fl = flda(train[,-21],train[,21])
knitr::kable(fl$t,caption="Confusion Matrix (Internal)")


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Internal Classification Error Rate =",fl$er,"%\n")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Shuffling data
set.seed(1652001)
s = sample(1:2000,1500,replace=F)
new_train = train[s,]
new_test = train[-s,]


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
d = lda(new_train[,-21],new_train[,21])
p = predict(d,new_test[,-21])
t = table(p$class,new_test[,21])


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::kable(t,caption="Confusion Matrix (External)")
cat("External Misclassification Error:",100*(sum(t)-sum(diag(t)))/nrow(new_test),"%")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loo=function(x,class){
  n=length(class)
  rslt={}
  for(i in 1:n){
     a = lda(x[-i,], class[-i])
     b = predict(a,x[i,])
     rslt[i]=b$class #[i]==class[i]
  }
return(rslt)
}


## ----echo = F---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# WARNING: This takes a lot of time to run
rslt = loo(train[,-21],train[,21])


## ----echo = F---------------------------------------------------------------------------------------------------------------------------------------------------------------------
t = table(rslt,train[,21])
knitr::kable(t,caption="Confusion Matrix (External - Leave One Out)")


## ----echo = F---------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Leave One Out External Misclassification Error:",100*(sum(t)-sum(diag(t)))/nrow(train),"%")


## ----include=T--------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(nnet)
mn = multinom(price_range ~ ., data = train)
results = predict(mn)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
t = table(train$price_range, results)
knitr::kable(t,caption="Confusion Matrix (Internal)")


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Internal Misclassification Error:",100*(sum(t)-sum(diag(t)))/nrow(train),"%")


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
mn = multinom(price_range ~ ., data = new_train)
results = predict(mn,new_test[,-21])
t = table(results,new_test[,21])


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::kable(t,caption="Confusion Matrix (External)")


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("External Misclassification Error:",100*(sum(t)-sum(diag(t)))/nrow(new_test),"%")


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
loo_mn=function(x,class){
  n=length(class)
  rslt={}
  for(i in 1:n){
     mn = multinom(price_range ~ ., data = train[-i,])
     b = predict(mn,x[i,])
     rslt[i]=b[1] #[i]==class[i]
  }
return(rslt)
}



## ----include=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------
# WARNING: This takes a lot of time to run
results_mn = loo_mn(train[,-21],train[,21])
t = table(results_mn,train[,21])


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::kable(t,caption="Confusion Matrix (External - Leave One Out)")


## ----echo = F---------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("External Misclassification Error:",100*(sum(t)-sum(diag(t)))/nrow(train),"%")

