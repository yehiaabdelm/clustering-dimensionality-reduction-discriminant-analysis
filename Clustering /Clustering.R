## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----include=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Hungarian Method Code
minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
  require(clue)
  idsA <- unique(clusteringA)  # distinct cluster ids in a
  idsB <- unique(clusteringB)  # distinct cluster ids in b
  nA <- length(clusteringA)  # number of instances in a
  nB <- length(clusteringB)  # number of instances in b
  if (length(idsA) != length(idsB) || nA != nB) {
    stop("number of clusters or number of instances do not match")
  }
 
  nC <- length(idsA)
  tupel <- c(1:nA)
 
  # computeing the distance matrix
  assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
  for (i in 1:nC) {
    tupelClusterI <- tupel[clusteringA == i]
    solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
      nA_I <- length(tupelA_I)  # number of elements in cluster I
      tupelB_I <- tupel[clusterIDsB == i]
      nB_I <- length(tupelB_I)
      nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
      return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
    }, clusteringB, tupelClusterI)
    assignmentMatrix[i, ] <- solRowI
  }
 
  # optimization
  result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
  attr(result, "assignmentMatrix") <- assignmentMatrix
  return(result)
}
# R2
R2=function(x,clusters,k=3){
  n=nrow(x); tss=var(x); tss=(n-1)*sum(diag(tss));
  wss=0
  for(j in 1:k){
    cj=x[clusters==j,]; nj=nrow(cj); vj=var(cj); wssj=0
    if(is.matrix(cj)) wssj=(nj-1)*sum(diag(vj)); wss=wss+wssj
  }
  r2=1-wss/tss; cat("R2 = ",r2,"\n")
  return(r2)
}


## ----echo=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------
df = read.csv("https://vincentarelbundock.github.io/Rdatasets/csv/DAAG/ais.csv",stringsAsFactors = T)
df = df[,-1]
head(df)


## ----echo=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------
levels(df$sport)


## ----echo=F,out.width="60%",fig.align='center',fig.cap="Pairs Plot"---------------------------------------------------------------------------------------------------------------
colors = c('black', 'red','blue','purple','pink','green','cornflowerblue','darkgoldenrod1','peachpuff3','turquoise3')[unclass(df$sport)]
pairs(df[,1:4],col=colors)


## ----echo=F,out.width="60%",fig.align='center',fig.cap="Pairs Plot"---------------------------------------------------------------------------------------------------------------
pairs(df[,5:8],col=colors)


## ----echo=F,out.width="60%",fig.align='center',fig.cap="Pairs Plot"---------------------------------------------------------------------------------------------------------------
colors_sex = c('black', 'red')[unclass(df$sex)]
pairs(df[,1:4],col=colors_sex)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
pow_end <- function(i) {
   if(i %in% c("B_Ball","Row","T_400m","Tennis","W_Polo")) {
        # endurance
        return(1)
   } else if(i %in% c("Gym","Netball","Swim","T_Sprnt","Field")) {
        # power
        return(2)
    }
}
p = sapply(df$sport,pow_end)
df$powend = p



## ----echo=F,out.width="60%",fig.align='center',fig.cap="Pairs Plot"---------------------------------------------------------------------------------------------------------------
pairs(df[,1:4],col=df$powend)


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# we are going to check the accuracy of the clustering algorithm with these 
sports = as.numeric(df$sport)
s = as.numeric(df$sex)
powend = df$powend


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
table(df$sport)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_sport = data.frame(df$sport)
sports_ = cbind(sports,df_sport)
count = c(25,37,23,22,19,29,15,11,4,17)
cbind(unique(sports_),count)


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# removing the last three columns that represent the classes
cluster_df = df[,-seq(12,14,1)]


## ----echo=F,out.width="60%",fig.align='center',fig.cap="Hierarchical Clustering (k=2)"--------------------------------------------------------------------------------------------
x = as.matrix(cluster_df)
x=scale(x, center=TRUE, scale=TRUE)
dissimilarity = dist(x,method='euclidian')
hc = hclust(dissimilarity,method='ward.D')
plot(hc)
clusters=cutree(hc,k=2)
rect.hclust(hc, k=2, border="red")


## ----echo=F,out.width="60%",fig.align='center',fig.cap="Hierarchical Clustering (k=2)"--------------------------------------------------------------------------------------------
# euclidian ward 2.9
# euclidian wardd2 3.9
# wucldian complete 6.9
# manhattan ward 5.9
# manhattan ward 4.9
dissimilaritye = dist(x,method='euclidian')
hce = hclust(dissimilaritye,method='ward.D2')
plot(hce)
clusterse=cutree(hce,k=2)
rect.hclust(hce, k=2, border="red")


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm=table(clusters,s); #  plot(x,pch=19,col=hd.clust)
print(cm)
#knitr::kable(cm,caption = "Confusion Matrix (Clusters Represent Gender)")
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")
R2(x,clusters,2)


## ----echo=F,out.width="60%",fig.align='center',fig.cap="Hierarchical Clustering (k=10)"-------------------------------------------------------------------------------------------
# manhattan ward d2 60
dissimilarity10 = dist(x,method='manhattan')
hc10 = hclust(dissimilarity10,method='ward.D2')
plot(hc10)
clusters10=cutree(hc10,k=10)
rect.hclust(hc, k=10, border="red")


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm=table(clusters10,sports);
match.label=minWeightBipartiteMatching(clusters10,sports); id=c(match.label); cm=cm[,id]
#knitr::kable(cm,caption = 'Confusion Matrix (Clusters Represent Sports)')
#cm=table(clusters10,sports);
print(cm)
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")
R2(x,clusters10,10)


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Splitting data into male and female
df_mf = df[,-c(14)]

male = df_mf[df_mf[,12]=='m',]
male_sports = as.numeric(male[,13])
male = male[,-c(12,13)]
male = as.matrix(male)
male = scale(male, center=TRUE, scale=TRUE)

female = df_mf[df_mf[,12]=='f',]
female_sports = as.numeric(female[,13])
female = female[,-c(12,13)]
female = as.matrix(female)
female = scale(female, center=TRUE, scale=TRUE)


## ----echo=F,out.width="60%",fig.align='center',fig.cap="Hierarchical Clustering (k=8)"--------------------------------------------------------------------------------------------
dissimilarity8male = dist(male,method='euclidian')
hc8male = hclust(dissimilarity8male,method='ward.D2')
plot(hc8male)
clusters8male=cutree(hc8male,k=8)
rect.hclust(hc8male, k=8, border="red")


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm=table(clusters8male,male_sports);
match.label = minWeightBipartiteMatching(clusters8male,male_sports); id = c(match.label); cm = cm[,id]
#knitr::kable(cm,caption = 'Confusion Matrix (Clusters Represent Male Sports)')
print(cm)
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")
R2(male,clusters8male,8)


## ----echo=F,out.width="60%",fig.align='center',fig.cap="Hierarchical Clustering (k=8)"--------------------------------------------------------------------------------------------
dissimilarity9female = dist(female,method='manhattan')
hc9female = hclust(dissimilarity9female,method='ward.D')
plot(hc9female)
clusters9female=cutree(hc9female,k=9)
rect.hclust(hc9female, k=9, border="red")


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm = table(clusters9female,female_sports);
match.label = minWeightBipartiteMatching(clusters9female,female_sports); id = c(match.label); cm = cm[,id]
#knitr::kable(cm,caption = 'Confusion Matrix (Clusters Represent Female Sports)')
print(cm)
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")
R2(female,clusters9female,9)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm = table(clusters,powend);
match.label = minWeightBipartiteMatching(clusters,powend); id = c(match.label); cm = cm[,id]
print(cm)
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")
R2(x,clusters,2)


## ----echo=F,out.width="60%",fig.align='center',fig.cap="Determining the Number of Clusters"---------------------------------------------------------------------------------------
# Number of Clusters
wss=c()
for (i in 2:10) {
  wss[i] = kmeans(x,centers=i)$tot.withinss}
plot(wss, type="b", pch=19, xlab="k",ylab="WSS", main="The L-Curve")


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
k=2; kmc = kmeans(x, k);
clusters=kmc$cluster
table(clusters)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm=table(clusters,s); #  plot(x,pch=19,col=hd.clust)
match.label = minWeightBipartiteMatching(clusters,s); id = c(match.label); cm = cm[,id]
print(cm)
#knitr::kable(cm,caption = "Confusion Matrix (Clusters Represent Gender - k=2)")
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")
R2(x,clusters,2)


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
k=10; kmc = kmeans(x, k);
clusters=kmc$cluster
table(clusters)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm=table(clusters,sports); #  plot(x,pch=19,col=hd.clust)
match.label = minWeightBipartiteMatching(clusters,sports); id = c(match.label); cm = cm[,id]
print(cm)
#knitr::kable(cm,caption = "Confusion Matrix (Clusters Represent Sports)")
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")
R2(x,clusters,2)


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
k=8; kmc = kmeans(male, k);
clusters=kmc$cluster
table(clusters)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm=table(clusters,male_sports);
match.label = minWeightBipartiteMatching(clusters,male_sports); id = c(match.label); cm = cm[,id]
print(cm)
#knitr::kable(cm,caption = "Confusion Matrix (Clusters Represent Male Sports - k=8)")
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")
R2(x,clusters,2)


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
k=9; kmc = kmeans(female, k);
clusters=kmc$cluster
table(clusters)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm=table(clusters,female_sports);
match.label = minWeightBipartiteMatching(clusters,female_sports); id = c(match.label); cm = cm[,id]
print(cm)
#knitr::kable(cm,caption = "Confusion Matrix (Clusters Represent Female Sports - k=9)")
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")
R2(x,clusters,2)


## ----include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
k=2; kmc = kmeans(x, k);
clusters=kmc$cluster
table(clusters)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm=table(clusters,powend); #  plot(x,pch=19,col=hd.clust)
match.label = minWeightBipartiteMatching(clusters,powend); id = c(match.label); cm = cm[,id]
print(cm)
#knitr::kable(cm,caption = "Confusion Matrix (Clusters Represent Gender)")
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")
R2(x,clusters,2)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
m = matrix(c(   '2.9%', '46%','60%','83%','56%',
                '3.9%','47%','65%','87%','56%',
                '2','2','10','8','9',
                '202','202','202','102','100'),ncol=5,byrow=T)
colnames(m) <- c('Sex','Power/Endurance','Sports','Sports (Male)','Sports (Female)')
rownames(m) <- c('Hierarchical','k-Means','k','n')
tab <- as.table(m)
tab
#knitr::kable(tab,caption="Comparing Methods Error Rate", booktabs=T)


## ----echo=F,out.width="70%",fig.align='center',fig.cap="Pairs Plot"---------------------------------------------------------------------------------------------------------------
colors = c('black', 'red','blue','purple','pink','green','cornflowerblue','darkgoldenrod1','peachpuff3','turquoise3')[unclass(df$sport)]
pairs(df[,9:11],col=colors)


## ----echo=F,out.width="70%",fig.align='center',fig.cap="Pairs Plot"---------------------------------------------------------------------------------------------------------------
pairs(df[,6:9],col=colors_sex)


## ----echo=F,out.width="70%",fig.align='center',fig.cap="Pairs Plot"---------------------------------------------------------------------------------------------------------------
pairs(df[,5:8],col=df$powend)


## ----echo=F,out.width="70%",fig.align='center',fig.cap="Pairs Plot"---------------------------------------------------------------------------------------------------------------
pairs(df[,9:11],col=df$powend)


## ----echo=F-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm=table(clusterse,s);
print(cm)
#knitr::kable(cm,caption = "Confusion Matrix (Clusters Represent Gender)")
cat(" Error rate =",1-sum(diag(cm))/sum(cm),"\n")

