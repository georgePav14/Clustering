

# functions #####

# Entropy
enfunc<-function(x){ if(any(is.nan(x*log(x)))){x[is.nan(x*log(x))]=1}
  sum(sum(x*log(x)))}

#hypothesis testing
hypT<-function(v,cl){
  for(i in 1:length(v)){
    print(paste('hypothesis testing for',names(domit[,v])[i]))
    print('ANOVA')
    print(summary(aov(domit[,i]~cl)))
    print(kruskal.test( domit[,i]~cl))
    print('-----------------------------------------')}
}

# multiple boxplots
boxp<-function(x,v){
  par(ask=TRUE)
  for(i in 1:dim(domit[,v])[2]){
    boxplot(dScale[,v][,i]~x,main=paste('Boxplot of clusters',names(d)[i+1],sep=' '))}
  par(ask=FALSE)}

# function to compute coefficient.
ac <- function(x) {
  agnes(dScale, method = x)$ac
}

#sil plot with different number of clusters
pSil<- function(x,y,z=dScale,meth='ward'){
  par(ask=TRUE)
  temp <- agnes(z, method = meth)
  for(i in y:x){
    temp2 <- cutree(temp, k = i)
    plot(silhouette(temp2, dist(dScale), method = "euclidean"))}
  par(ask=FALSE)}

# plot for average silhouette for different number of clusters
psil2<-function(x,y,z=dis,h=hWard){
  res<-NULL
  temp3<-cutree(h,x:y)
  for (i in (x-1):(y-1)){
    a<-silhouette(temp3[,i],z)
    res<-c(res,mean(a[,3]))}
  plot(x:y,res,type="b",ylim=c(0,0.5),xlab="clusters", main="Average Silhouette")}


### packages ######

install.packages(c('cluster','factoextra','tidyverse','dendextend',"corrgram",
                   "HDclassif","cluster","mclust","FactMixtAnalysis",
                   "nnet","class","tree", "pgmm",'corrplot','car','dbscan'))
install.packages('FactMixtAnalysis')

library(dbscan)
library(tidyverse) # data manipulation
library(dendextend) # for comparing two dendrograms
library(cluster) # clustering algorithms
library(factoextra)
library(dplyr)
library(corrgram)
library(corrplot)
library(HDclassif)
library(cluster)
library(mclust)
library(FactMixtAnalysis)
library(nnet)
library(class)
library(tree)
library(car)

load(file.choose() ) # load the data
dim(alld)
apply(alld,2,function(x) sum(is.na(x))) # Ποσα ΝΑ εχει καθε μεταβλητη

d<-alld[,c(1,seq(131,145,by=2),146,147)] ; View(d) # Επιλεγω τις μεταβλητες που θέλω

head(d)
apply(d,2,function(x) sum(is.na(x))) # Ποσα ΝΑ εχει καθε μεταβλητη


## Data cleaning ####
corrplot(cor(na.omit(d[,-1])),method='ellipse')

dslim<-d[,2:9] ; head(dslim) # remove the higly correlated variables
domit<-na.omit(dslim)
dScale<-scale(domit[,-1])

#### some plots of the data ####
corrplot(cor(dScale),method='ellipse')

r<-cor(dScale,use= 'pairwise.complete.obs')
corrplot(r, method ="number",type = "lower",diag=F, tl.cex=1, cl.cex=1)

pairs(dScale)

# All the variables with omited na lines

# Hierarchical ####

# Dissimilarity matrix
dis <- dist(dScale, method = "euclidean")

# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# agglomerative coefficient, which measures the
# amount of clustering structure found (values
# closer to 1 suggest strong clustering structure).

map_dbl(m, ac)

# average single complete ward
# 0.9406398 0.9201809 0.9572674 0.9922642

hComplete <- agnes(dis, method = "complete")
plot(hComplete)

hAverage <- agnes(dScale, method = "average")
hSingle <- agnes(dScale, method = "single")
hWard <- agnes(dScale, method = "ward")

#### Plot the obtained dendrogram
plot(hComplete,which.plots = 2) # complete not looking good chain effect
plot(hAverage,which.plots = 2) # also chain effect
plot(hSingle,which.plots = 2) # much more chain effect
plot(hWard,which.plots = 2) # better I think 4 or less clusters
rect.hclust(hWard, k = 4, border = 2:5)

## Ward with hclust
hCward <- hclust(dis, method = "ward.D2" )
plot(hCward) # similar but not exacly
rect.hclust(hCward, k = 4, border = 2:5)

# Cut tree into 4 groups ####
sub_grpW <- cutree(hWard, k = 4)
boxp(sub_grpW) # compare cluster in each variable

# Number of members in each cluster ####
table(sub_grpW)

### Silhouette #####
plot(silhouette(sub_grpW, dis))
sil<-silhouette(sub_grpW, dis)
mean(sil[,3])

# try with different number of clusters
pSil(2,4)

# 2 is better
sub_grpW2 <- cutree(hWard, k = 2)
plot(silhouette(sub_grpW2, dis))
sil2<-silhouette(sub_grpW2, dis)
mean(sil2[,3])

######### wilks lambda ##########
clas1<-cutree(hWard,3:5)
m <- manova(dScale~clas1[,3])
plot(m$residuals[,10],type='line')

summary(m,test="Wilks",tol=0) # απορριπτω την υποθεση οτι οι μεσοι των
# cluster ειναι ισοι

### a kind of variable selection ####
allcomb<-t(combn(1:dim(dScale)[2],2)) # all the combination of the ten variables
# in two clusters

res<-NULL
for (i in 1:dim(allcomb)[1]){
  print(i)
  x0<-dScale[,allcomb[i,]]
  mydist<- as.dist(apply(x0, 1, function(i) mahalanobis(x0, i, cov = cov(x0))))
  hc2<-hclust(mydist,method="ward.D2")
  clas1<-cutree(hc2, k=2)
  a<-silhouette(clas1, mydist)
  res<-c(res,mean(a[,3]))
}

### so now we do have the best model with 2 variables
### we can go on with some forward algorithm
### given the 3 variables which is the best one to add to make a quartet?
### and life goes on...
### here your turn starts
plot(res)

### find the best but ... surprise !!!
best<-allcomb[which.max(res),]
names(d[,-1])[best]
apply(d,2,function(x) sum(is.na(x))) # Ποσα ΝΑ εχει καθε μεταβλητη

### the average silhouette of the best combn
x0<-dScale[,best]
mydist<- as.dist(apply(x0, 1, function(i) mahalanobis(x0, i, cov = cov(x0))))
hc2<-hclust(mydist,method="ward.D2")
clas1<-cutree(hc2, k=2)
a<-silhouette(clas1, mydist)
plot(a)
plot(hc2)

### Plot the results ####
co<-cbind(allcomb,res)
co<-co[order(co[,1]),]

ma<-matrix(rep(1,100),ncol=10,nrow=10)
count<-1
for(i in 1:9){
  for(j in min(co[co[,1]==i,2]):max(co[co[,1]==i,2])){
    ma[i,j]<-res[count]
    count<-count+1
  }
}

ma2<-t(ma)*ma
corrplot(cor(dScale),method='ellipse')
corrplot(ma2,method='ellipse')

### Hierarchical Clustering with Mahalanobis distance ####

mydist<- as.dist(apply(dScale, 1, function(i) mahalanobis(dScale, i, cov = cov(dScale))))

hcla<-hclust(mydist, method="ward.D2")
hcla<-hclust(mydist, method="complete")
clas<-cutree(hcla,3)
plot(hcla)
rect.hclust(hclab,3)
plot(silhouette(clas,mydist))

#### K-means #####

sol<-kmeans(dScale,3)
attributes(sol)
sol$withinss
sol$betweenss
fviz_cluster(sol, geom = "point", data = dScale) + ggtitle("k = 3")
### run with 2 up to15 clusters
### and plot the within ss
within<-NULL

for (i in 2:15) {
  within<-c(within,kmeans(dScale[,2:8],i,nstart=20)$tot.withinss) }
plot(2:15,within, type="b",xlab="number of cluster", ylab="total within ss")

######### DBScan ########

knn<-kNNdist(dScale, 2)
mean(knn)

kNNdistplot(dScale, 2)
abline(h=1)

ds <- dbscan(dScale, 1)
plot(silhouette(ds$cluster,dist(dScale)))
ds$cluster

#### Determining Optimal Clusters ####

# Elbow method
fviz_nbclust(dScale, FUN = kmeans, method = "wss",diss=dis,nstart=20) +
  geom_vline(xintercept = 4, linetype = 2) + labs(subtitle = "Elbow method") # 4


fviz_nbclust(dScale, FUN = hcut, method = "wss",diss=dis) +
  geom_vline(xintercept = 4, linetype = 2) + labs(subtitle = "Elbow method") #4

# Silhouette method
fviz_nbclust(dScale, kmeans, method = "silhouette",diss=dis)+
  labs(subtitle = "Silhouette method")

fviz_nbclust(dScale, hcut, method = "silhouette",diss=dis)+
  labs(subtitle = "Silhouette method")

# Mahalanobis distance hierarhical
hc2<-hclust(mydist,method="ward.D2")
clas1<-cutree(hc2, k=4)
a<-silhouette(clas1, mydist)
plot(a)
plot(hc2)

# eyclidian distance hierarhical
hc2<-hclust(dis,method="ward.D2")
clas1<-cutree(hc2, k=2)
a<-silhouette(clas1, dis)
plot(a)
plot(hc2)
mean(a[,3]) # 0.22->4 ,0.44->2

# kmean withinss

within<-NULL
for (i in 2:15) {
  within<-c(within,kmeans(dScale,i,nstart=20)$tot.withinss) }
plot(2:15,within, type="b",xlab="number of cluster", ylab="total within ss")


# Gap statistic
# nboot = 50 to keep the function speedy.
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(dScale, kmeans, nstart = 25, method = "gap_stat", k.max = 10, nboot = 50)+
  labs(subtitle = "Gap statistic method")

gap_stat <- clusGap(dScale, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)


#### MCLUST ####


mc1<-Mclust(dScale, G=2:10, modelNames=c('EII',"VII", "EVE", "EVI", "VEI",'VEV','EVV', "VVI",'VVV'))

#mc1<-Mclust(dScale[,t(combn3)[2,]], G=2:10, modelNames=c('EII',"VII", "EVE", "EVI", "VEI",'VEV','EVV', "VVI",'VVV'))

mc1$G#The optimal number of mixture components.
mc1$BIC#All BIC values. Bigger is better
mc1$bic#Optimal BIC value.
mc1$loglik#The loglikelihood corresponding to the optimal BIC.
mc1$df #The number of estimated parameters.
mc1$parameters#A list with the following components:
mc1$z #posterior probabilities.
mc1$classification #map(z): The classification corresponding to z.
mc1$uncertainty#The uncertainty associated with the classification.
plot(mc1)
summary(mc1)

attributes(mc1)
cbind(round(mc1$z,6), round(mc1$uncertainty,6))
mc1$parameters$pro
table(mc1$classification)/dim(dScale)[1]
# head(cbind('P(i=k1)'=round(mc1$z,6)[,1],'P(i=k2)'=round(mc1$z,6)[,2], uncertainty=round(mc1$uncertainty,6)))


# I try all the combinations of 3 variables and compare them with bic
# I think that for whatever mixture i use the best number of clusters is the same
# look at plot(mc1) BIC. For my needs that is of find the best combination I dont care
# for the best mixture, so in order to run faster i will use the simplest EII
# I dont know whether or not i am correct ???

#t(combn3)
combn3<-combn(1:dim(dScale)[2],3)
g3=bic3=msil3=c()
res3<-matrix(rep(NA,dim(t(combn3))[1]*(3+2)),nrow=dim(t(combn3))[1])

for( i in 1:dim(t(combn3))[1]){
  mc2<-Mclust(dScale[,t(combn3)[i,]], G=2:10, modelNames=c('EII'))
  msil3<-c(msil3,mean(silhouette(mc2$classification,dis)[,3]))
  bic3<-c(bic3,mc2$bic)
  g3<-c(g3,mc2$G)}

res3<-cbind(t(combn3)[1:i,],BIC=bic3,Clusters=g3,meanSil=msil3)

res3<-res3[rev(order(res3[,3+1])),] ; res3 # 349

table(res3[1:10,1:3])
head(dScale[,res3[which.max(res3[,4]),1:3]]) # the 'best' variables res3[which.max(res3[,4]),1:3]

best3<-res3[which.max(res3[,3+1]),1:3]

best3<-res3[5,1:3]

#
mc<-Mclust(dScale[,best3], G=2:10, modelNames=c('VVV'))
plot(silhouette(mc$classification,dis))
mean(silhouette(mc$classification,dis)[,3])
mc$classification
#### Manual Forward Step. Not backward because ideally i want less variables ####
best4=best4

# The combinations of variables with old ones i selected
combn4<-t(combn(1:dim(dScale)[2],4))
var<-c()
for(i in 1: dim(combn4)[1]){
  tVar<-all(best4 %in% combn4[i,])
  var<-c(var,tVar)
}
v4<-combn4[var,]

temp4=g4=bic4=msil4=c()
res4<-matrix(rep(NA,dim(t(combn4))[1]*(4+2)),nrow=dim(t(combn4))[1])
for(i in 1:dim(v4)[1]){
  mc3<-Mclust(dScale[,v4[i,]], G=2:10, modelNames=c('VVV'))
  bic4<-c(bic4,mc3$bic)
  msil4<-c(msil4,mean(silhouette(mc3$classification,dis)[,3]))
  g4<-c(g4,mc3$G)
}
res4<-cbind(v4,BIC=bic4,Clusters=g4,meanSil=msil4)

res4<-res4[rev(order(res4[,4+1])),] ; res4
table(res4[,1:4]) # the most common variables
head(dScale[,res4[which.max(res4[,4]),1:4]]) # the 'best' variables

best4<-res4[which.max(res4[,4+1]),1:4] # variables with max bic
best4<-res4[1,1:4] # select variables i want

#######################TEST###################

combn2<-t(combn(1:dim(dScale)[2],2))
best3<-c(2,4,6);best4<-c(2,3,4,6);best5<-c(1,3,4,6,7);best6<-c(1,3,4,5,6,7);best7<-c(1,2,3,4,5,6,7)


# plot the classification of the hierarchical
hCward <- hclust(dist(dScale[,best4], method = "euclidean"), method = "ward.D2" )
plot(silhouette(cutree(hCward,2 ),dist(dScale[,best4])),main = 'Silhouette plot') 
plot(hCward)
rect.hclust(hCward, k = 2, border = 2:5)
pairs(domit[,best4],col=cutree(hCward,2 )+1)

boxp(cutree(hCward,2 ),best4)

#hypothesis testing Anova
hypT(best5,cutree(hCward,2))

# model based clustering and testing
mc<-Mclust(dScale[,res4[rev(order(res4[,4+1])),][1,1:4]], G=2:10, modelNames=c('VVV'))
plot(silhouette(mc$classification,dist(dScale[,best4])))
mc$classification

# trying different models and clusters
mc<-Mclust(dScale[,best4], G=2:10, modelNames=c("EII", "VII", "EEI", "EVI", "VEI",'EVE','VVE', "VVI",'VVV'))
plot(mc)
summary(mc)
mc$G
mc$BIC

### Entropy ####

totEn<-c()
mc0<-Mclust(dScale[,res4[rev(order(res4[,4+1])),][1,1:4]], G=1, modelNames=c('VVV'))
for(i in 2:10 ){
  mc<-Mclust(dScale[,best4], G=i, modelNames=c('VVV'))
  entropy<-(-1)*sum(apply(mc$z,2,enfunc))/(mc$loglik-mc0$loglik) ; entropy
    totEn<-c(totEn,entropy)
}
plot(2:10,totEn,pch=18,type='b',cex=1.5,ylim=c(0,max(totEn)),ylab='Entropy',xlab='Clusters'); abline(h=0,col='red',lwd=2)

### plot the classification of the model based
mc<-Mclust(dScale[,best4], G=2, modelNames=c('VVV'))
plot(silhouette(mc$classification,dist(dScale[,best4])))
pairs(domit[,best5],col=mc$classification+1)

boxp(mc$classification,best5)

#hypothesis testing Anova
hypT(best5,mc$classification)

# kmeans
sol<-kmeans(dScale[,best4],2)
attributes(sol)

plot(silhouette(sol$cluster,dist(dScale[,best4])),main = 'Silhouette plot k-means')

### plot the classification of the kmeans
pairs(domit[,best4],col=sol$cluster)
boxp(sol$cluster,best4)

#hypothesis testing Anova
hypT(best4,sol$cluster)

#dbScan
knn<-kNNdist(dScale[,best4], 2)
mean(knn)

kNNdistplot(dScale[,best4], 2,minPts=2*length(best4))
abline(h=1)

ds <- dbscan(dScale[,best4], 1,minPts=2*length(best4))
plot(silhouette(ds$cluster,dist(dScale[,best4])))
ds$cluster
### plot the classification of the dbscan
pairs(domit[,best4],col=ds$cluster+1)
boxp(ds$cluster,best4)

#hypothesis testing Anova
hypT(best5,ds$cluster)

##### Adjusted Rand Index####

hCward <- hclust(dist(dScale[,best4], method = "euclidean"), method = "ward.D2" )
clas4<-cutree(hCward,2)

adjustedRandIndex(clas5,clas4) # 5 with 7 close with 0.73

# eyclidian distance hierarhical
hCward <- hclust(dist(dScale[,c(4,9,5)], method = "euclidean"), method = "ward.D2" )
plot(silhouette(cutree(hCward,2),dist(dScale[,best4]))) # best5 1349
plot(hCward)


# Elbow method
fviz_nbclust(dScale[,best4], FUN = kmeans, method = "wss",diss=dis,nstart=20) +
  geom_vline(xintercept = 4, linetype = 2) + labs(subtitle = "Elbow method") # 4


fviz_nbclust(dScale[,best4], FUN = hcut, method = "wss",diss=dis) +
  geom_vline(xintercept = 4, linetype = 2) + labs(subtitle = "Elbow method") #4

# Silhouette method
fviz_nbclust(dScale[,best4], kmeans, method = "silhouette",diss=dis)+
  labs(subtitle = "Silhouette method")

fviz_nbclust(dScale[,best4], hcut, method = "silhouette",diss=dis)+
  labs(subtitle = "Silhouette method")

# kmean withinss
within<-NULL
for (i in 2:15) {
  within<-c(within,kmeans(dScale[,best4],i,nstart=20)$tot.withinss) }
plot(2:15,within, type="b",xlab="number of cluster", ylab="total within ss",main = 'Whithin Variance')

# Gap statistic
set.seed(123)
fviz_nbclust(dScale[,best4], kmeans, nstart = 25, method = "gap_stat", k.max = 10, nboot = 50)+
  labs(subtitle = "Gap statistic method")

gap_stat <- clusGap(dScale[,best4], FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

### Intepretation  ####
fd<-domit[,best4+1] ; names(fd) ; View(fd)
r<-apply(d[,-1],1,function(x) any(is.na(x)))

ad<-alld[!r,1:130]

ad$clusters<-sol$cluster
#ad<-na.omit(ad)

#plots
  t=allt=c()
  for(i in 1:dim(ad[,c(-1,-2)])[2]){
  t<-aov(ad[,c(-1,-2)][,i]~ad$clusters)
  allt<-c(allt,summary(t)[[1]][["Pr(>F)"]][1])}
  
  t=allt=c()
  for(i in 1:dim(ad[,c(-1,-2)])[2]){
    t<-kruskal.test(ad[,c(-1,-2)][,i]~ad$clusters)
    allt<-c(allt,t$p.value)}

    length(allt[allt<0.05])
  
  par(ask=TRUE)
  for(i in 1:dim(ad[,c(-1,-2)][,allt<0.05])[2]){
    boxplot(ad[,c(-1,-2)][,allt<0.05][,i]~ad$clusters,ylab=names(ad[,c(-1,-2)][,allt<0.05])[i],main=paste('Boxplot of clusters',names(ad[,c(-1,-2)][,allt<0.05])[i],sep=' '))
  print(c(i,names(ad[,c(-1,-2)][,allt<0.05])[i]))}
  par(ask=FALSE)
