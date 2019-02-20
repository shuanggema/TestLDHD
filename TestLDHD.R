
ld.test<-function(X,Z,y,lambda1,lambda2,g1,g2,w){
  p=ncol(Z)
  pval=rep(0,ncol(X))
  beta_set=rep(0,ncol(X))
  theta_hat=matrix(0,ncol(X),ncol(Z))
  
  
  for (ll in 1:ncol(X)){
    
    
    XX=cbind(X[,ll],Z)
    
    fit=glmnet(XX,y,penalty.factor = c(0,rep(1,p)),weights = sqrt(w),standardize = FALSE,intercept=FALSE,lambda = lambda1)
    s=rep(0,ncol(fit$beta))
    y_hat<-matrix(0,n,ncol(fit$beta))
    for (i in 1:ncol(fit$beta)) {
      y_hat[,i]=XX%*%fit$beta[,i]
      tmp=0
      for (j in 1:n) {
        tmp=tmp+w[j]*(y[j]-y_hat[j,i])^2
      }
      s[i]=tmp
    }
    ebic=log(s/n)+(log(n)*fit$df+2*g1*log(choose(p,fit$df)))/n
    estimate=as.matrix(fit$beta[,which.min(ebic)])
    
    beta=as.matrix(estimate[1],1,1)
    theta=as.matrix(estimate[2:(p+1)],p,1)
    theta_hat[ll,]=theta
    
    fit2=glmnet(Z,X[,ll],weights = sqrt(w),standardize = FALSE,intercept=FALSE,lambda = lambda2)
    s=rep(0,ncol(fit2$beta))
    y_hat<-matrix(0,n,ncol(fit2$beta))
    for (i in 1:ncol(fit2$beta)) {
      y_hat[,i]=Z%*%fit2$beta[,i]
      tmp=0
      for (j in 1:n) {
        tmp=tmp+w[j]*(X[j,ll]-y_hat[j,i])^2
      }
      s[i]=tmp
    }
    ebic=log(s/n)+(log(n)*fit2$df+2*g2*log(choose(p,fit2$df)))/n
    
    
    tilde_B=as.matrix(fit2$beta[,which.min(ebic)])
    
    
    tilde_X=X[,ll]-Z%*%tilde_B
    
    tilde_beta=solve(t(tilde_X)%*%W%*%X[,ll])%*%t(tilde_X)%*%W%*%(y-Z%*%theta)
    
    sigma_hat<-midSig(y,delta,X[,ll],tilde_X,Z,tilde_beta,theta)
    
    sigma_hat_new<-solve(t(tilde_X)%*%W%*%X[,ll]/n)*sigma_hat*solve(t(X[,ll])%*%W%*%tilde_X/n) 
    se=sqrt(sigma_hat_new/n)
    
    t_stat<-abs(tilde_beta/se)
    pval[ll]<-2*(1-pt(t_stat,n-1))   
    
    
    beta_set[ll]<-tilde_beta
  }
  
  
  
  idbb=which(p.adjust(pval,'BY')<0.05)
  result=cbind(colnames(X[,idbb]),pval[idbb],p.adjust(pval[idbb]),beta_set[idbb],theta_hat[idbb,])
  return(result)
}



gamma0 = function(Y, delta)
{
  g = rep(NA, length(Y))
  k = 1
  for (y in Y)
  {
    index0 = (delta == F) & (Y < y)
    c0 = 0
    n = length(Y)
    for (i in Y[index0])
    {
      c0 = c0 + 1 / n / (1 - 1 / n * sum(Y <= i))
    }
    g[k] = exp(c0)
    k = k + 1
  }
  g = g[!is.na(g)]
  g
}

phi = function(j, Xtilde, Y, X, Z, beta, theta)
{
  resid = Y - crossprod(t(X), beta) - crossprod(t(Z), theta)
  as.numeric(Xtilde[, j] * resid)
}

gamma1 = function(y, Y, delta, j, Phi, Gamma0)
{
  index0 = (delta == T) & (Y > y)
  n = length(Y)
  Phi[index0, j] %*% Gamma0[index0] / (n - sum(Y <= y))
}

gamma2 = function(y, Y, delta, j, Phi, Gamma0)
{
  index1 = (delta == F) & (Y < y)
  c0 = 0
  n = length(Y)
  for (i in Y[index1])
  {
    index2 = (delta == T) & (Y > i)
    c0 = c0 + (Phi[index2, j] %*% Gamma0[index2]) / (n - sum(Y <= i)) ^ 2
  }
  c0
}


zi = function(Y, delta, j, Phi, Gamma0)
{
  n = length(Y)
  vec = 1:n
  for (i in 1:n)
  {
    vec[i] = gamma1(Y[i], Y, delta, j, Phi, Gamma0) * (1 - delta[i]) - 
      gamma2(Y[i], Y, delta, j, Phi, Gamma0)
  }
  vec = vec + Phi[1:n, j] * Gamma0[1:n] * delta
  vec
}

midSig = function(Y, delta, X, Xtilde, Z, beta, theta)
{
  X = as.matrix(X)
  Zi = matrix(nrow = nrow(X), ncol = ncol(X))
  Phi = matrix(nrow = nrow(X), ncol = ncol(X))
  Gamma0 = gamma0(Y, delta)
  for (j in 1:ncol(X))
    Phi[, j] = phi(j, Xtilde, Y, X, Z, beta, theta)
  for (j in 1:ncol(X))
    Zi[, j] = zi(Y, delta, j, Phi, Gamma0)
  Zi = Zi[rowSums(is.nan(Zi)) == 0, ]
  cov(Zi, Zi)
}
kmw<-function(y,delta){
  y_s=y
  delta_s=delta
  kmweight<-c()
  nw<-length(y)
  
  comb<-cbind(y,delta)
  oo=order(y)
  ocomb<-comb[oo,]
  y<-ocomb[,1]
  delta<-ocomb[,2]
  kmweight[1]<-delta[1]/nw
  for(ind in 2:nw){
    tmp<-c()
    for(ind2 in 1:(ind-1)){
      tmp[ind2]<-((nw-ind2)/(nw-ind2+1))^delta[ind2]
    }
    kmweight[ind]<-delta[ind]/(nw-ind+1)*prod(tmp)
  }
  kmweight1=matrix(0,nw,0)
  kmweight1[oo]=kmweight
  return(kmweight=kmweight1)
}