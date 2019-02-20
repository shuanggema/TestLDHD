rm(list = ls())
load('luad_data.RData')
library(glmnet)
library(survival)
source('TestLDHD.R')

#Identification of gene expressions with independent prognostic power conditional on imaging features

lambda1 = seq(0.1,0.28,0.005)

g1=0.02

lambda2 = seq(0.03,0.3,0.003)
 
g2<-0.6

n=length(y)
w=kmw(y,delta)
w=w/sum(w)*n


#normalize
ymean=sum(w*y)/sum(w)
gmean=colSums(X*(w%*%matrix(1,1,ncol(X))))/sum(w)
zmean=colSums(Z*(w%*%matrix(1,1,ncol(Z))))/sum(w)

y=y-ymean
X=X-matrix(1,n,1)%*%gmean
Z=Z-matrix(1,n,1)%*%zmean
X=scale(X)
Z=scale(Z)

gene=ld.test(X,Z,y,lambda1,lambda2,g1,g2,w)

##survival analysis

idbb=dim(gene)[1]
fit.KM=survfit(my.surv~1)

#plot(fit.KM)

for (ll in idbb){
  
  fit=glm(y~X[,ll]-1,weights = sqrt(w))
  
  
  beta=fit$coefficients
  
  beta_set[ll]<-beta
}
beta=beta_set[idbb]

yg=matrix(0,n,length(idbb))
for (i in 1:length(idbb)) {
  yg[,i]= Z%*%as.numeric(gene[i,5:dim(gene)[2]])
}
log_rank_pg <- apply(yg, 2, function(values1){
  group=ifelse(values1>median(values1),'low','high')
  kmfit2 <- survfit(my.surv~group)
  #plot(kmfit2)
  data.survdiff=survdiff(my.surv~group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
})

y1=matrix(0,n,length(idbb))
for (i in 1:length(idbb)) {
  y1[,i]= X[,gene[i,1]]*as.numeric(gene[i,4])+Z%*%as.numeric(gene[i,5:dim(gene)[2]])
}
log_rank_p <- apply(y1, 2, function(values1){
  group=ifelse(values1>=median(values1),'low','high')
  kmfit2 <- survfit(my.surv~group)
  #plot(kmfit2)
  data.survdiff=survdiff(my.surv~group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
})
sum(log_rank_pg<0.05)
sum(log_rank_p<0.05)

