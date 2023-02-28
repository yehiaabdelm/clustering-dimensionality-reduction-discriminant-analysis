knitr::opts_chunk$set(echo = TRUE)

df = read.csv("Life Expectancy Data.csv")

df_2014 = df[df[,2]==2014,] # taking only observations from 2014
df_2014 = df_2014[ , colSums(is.na(df_2014)) == 0] # removing columns with NA
df_2014$Developing = model.matrix( ~ Status - 1, data=df_2014 )[,2] # creating binary variable
df_2014 = df_2014[,-c(2,3)]
colnames(df_2014) <- c("Country","LifeExp","AdultMort","InfDeaths","PercExp","Measles","UnderFive","Polio","Diphtheria","HIV.AIDS","Developing")
row.names(df_2014) <- NULL


pairs(df_2014[,c(2,3,5)]) # LifeExp, AdultMort, PercExp

pairs(data.frame(df_2014[,2],log(df_2014[,5])),labels = c("LifeExp","log(PercExp)"))

pairs(df_2014[,c(3,10)]) #AdultMort vs HIV.AIDS

colMeans(df_2014[,2:11])

plot(df_2014[,c(4,7)]) # InfDeaths vs Underfive

plot(df_2014[,c(8,9)])



cor = round(cor(df_2014[,2:11]),2)
cor

reg1 = lm(df_2014$LifeExp~df_2014$HIV.AIDS)
plot(df_2014$HIV.AIDS,df_2014$LifeExp,xlab="HIV.AIDS",ylab="LifeExp")
abline(reg=reg1)

reg2 = lm(df_2014$LifeExp~df_2014$Developing)
plot(df_2014$Developing,df_2014$LifeExp,xlab="Developing",ylab="LifeExp")
abline(reg=reg2)

reg3 = lm(df_2014$InfDeaths~df_2014$Measles)
plot(df_2014$InfDeaths,df_2014$Measles,xlab="InfDeaths",ylab="Measles")
abline(reg3)

# Eliptical Distances
edist = function(x){
  x = as.matrix(x)
  ev = eigen(var(x))
  d = diag(1/sqrt(abs(ev$values)))
  b = ev$vectors%*%d%*%t(ev$vectors)
  z = x%*%b
  d = dist(z);  return(d) }
ed = edist(df_2014[,2:11])

mahalanobis_dis=function(df){
  dist = c()
  sigma_inv = solve(cov(df))
  xbar = colMeans(df)
  for(i in 1:nrow(df)){
    x_xbar = as.matrix(df[i,]-xbar)
    dist = c(dist,sqrt(x_xbar%*%(sigma_inv)%*%t(x_xbar)))
  }
  return(dist)
}

MD = mahalanobis_dis(df_2014[,2:11])
df_2014$MD = MD

sorted_MD = df_2014
sorted_MD = sorted_MD[order(-MD),]
row.names(sorted_MD) <- NULL
knitr::kable(head(sorted_MD[,c(1,2,11,12)],20),caption = "Data with Mahalanobis Distances") %>%
    kable_styling(latex_options = "HOLD_position")

e2d = function(loc, cov, dis) { 
A = solve(cov) 
eA = eigen(A) 
ev = eA$values 
lambda1 = max(ev) 
lambda2 = min(ev) 
eigvect = eA$vectors[, order(ev)[2]] 
z = seq(0, 2 * pi, by = 0.01) 
z1 = dis/sqrt(lambda1) * cos(z) 
z2 = dis/sqrt(lambda2) * sin(z) 
alfa = atan(eigvect[2]/eigvect[1]) 
r <- matrix(c(cos(alfa), -sin(alfa), sin(alfa), cos(alfa)), ncol = 2) 
z = t(t(cbind(z1, z2) %*% r) + loc) 
}

ellipse=function(x,center,cov,dis,col,type=1){
  if(class(x)== "data.frame") x=as.matrix(x)
  if(class(x)!="matrix") stop("Data is not a matrix.")
  if(ncol(x)!=2) stop("Number of variables should be 2.")
  if(type==1) dis = sqrt(qchisq(dis, 2)) 
  y=e2d(center,cov,dis); z=rbind(x,y) 
  plot(x, pch=19,cex=0.8,col=col,xlim=c(min(z[,1]),max(z[,1])),ylim=c(min(z[,2]),max(z[,2])) )
  lines(y,col=4)
}

x = df_2014[,c(2,3)]
cov = cov(x)
eig = eigen(cov)
mu = colMeans(x)
ellipse(x,mu,cov,col=df_2014[,11]+2,max(mahalanobis_dis(df_2014[,c(2,3)])),type=2)
legend(90,600,legend=c("Developing","Developed"),fill=c(3,2))
abline(a=910, b=-10.5)
abline(a=150, b=-1.7)

require(robustX); library(robustbase)
x = as.matrix(df_2014[,2:10])
output = mvBACON(x)

df_2014_bd = df_2014
df_2014_bd$bd = output$dis

top_countries = rbind(df_2014[df_2014[,"Country"]=="Nigeria",],df_2014[df_2014[,c("Country")]=="India",],df_2014[df_2014[,c("Country")]=="China",],df_2014[df_2014[,c("Country")]=="Philippines",],df_2014[df_2014[,c("Country")]=="Democratic Republic of the Congo",])
top_countries

y = cbind(1:99,df_2014_bd$bd)
colnames(y) <- c("Index","Distance")
plot(y[,1:2],xlab="Index",ylab="BD", pch=19, main = "BACON Distances")
abline(h=output$limit,  col= "red",  lty=2)

t(output$center)

robustcor = diag(1/sqrt(diag(output$cov)))%*%output$cov%*%diag(1/sqrt(diag(output$cov)))
colnames(robustcor) <- c("LifeExp", "AdultMort" ,"InfDeaths" ,"PercExp" ,"Measles", "UnderFive", "Polio", "Diphtheria" ,"HIV.AIDS")
row.names(robustcor) <- c("LifeExp", "AdultMort" ,"InfDeaths" ,"PercExp" ,"Measles", "UnderFive", "Polio", "Diphtheria" ,"HIV.AIDS")
round(robustcor,2)

ht2=function(x, y) {
n=dim(x)[1];m=dim(y)[1]; p=dim(x)[2]
xcov=cov(x);
ycov=cov(y)
Sp=(n-1)*xcov+(m-1)*ycov; Sp=Sp/(n+m-2)
xcenter=colMeans(x); ycenter=colMeans(y)
d=xcenter-ycenter 
T2=t(d)%*%solve(Sp)%*%d 
T2=T2*n*m/(n+m) 
F=T2*(n+m-p-1)/(p*(n+m-2)) 
pv=1-pf(F,p,n+m-p-1) 
list(xcenter=xcenter,ycenter=ycenter,xcov=xcov,ycov=ycov, Sp=Sp,T2=T2,F=F,df=c(p,n+m-p-1),pv=pv) }

library(Hmisc)

hist.data.frame(df_2014[,2:11])

out1 = ht2(df_2014[df_2014$Developing==1,2:10],df_2014[df_2014$Developing==0,2:10])

out1$pv

df_2014_woo = df_2014[output$subset,]

out2 = ht2(df_2014_woo[df_2014_woo$Developing==1,2:10],df_2014_woo[df_2014_woo$Developing==0,2:10])

t(out2$xcenter)

t(out2$ycenter)


out2$pv

tf = df_2014
tf$AdultMort = log(df_2014$AdultMort+1)
tf$InfDeaths = log(df_2014$InfDeaths+1)
tf$PercExp = log(df_2014$PercExp+1)
tf$Measles = log(df_2014$Measles+1)
tf$UnderFive = log(df_2014$UnderFive+1)
tf$Polio = (df_2014$Polio)^8
tf$Diphtheria = (df_2014$Diphtheria)^6
#tf = as.matrix(tf[,2:10])

#out3 = ht2(tf[tf$Developing==1,2:10],tf[tf$Developing==0,2:10])
tfm = as.matrix(tf[,2:10])

eig = eigen(t(tfm)%*%tfm)

cat(sqrt(eig$values[1]/eig$values[2]),sqrt(eig$values[1]/eig$values[3]),sqrt(eig$values[1]/eig$values[4]),sqrt(eig$values[1]/eig$values[5]),sqrt(eig$values[1]/eig$values[6]),sqrt(eig$values[1]/eig$values[6]),sqrt(eig$values[1]/eig$values[7]),sqrt(eig$values[1]/abs(eig$values[7])),sqrt(eig$values[1]/abs(eig$values[8])),sqrt(eig$values[1]/abs(eig$values[9])))
