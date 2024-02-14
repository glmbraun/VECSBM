library(gtools)
library(dplyr)
library(RSpectra)
library(mvnfast)
library(Rcsdp)
library(MASS)
library(aricode)
library(expm)
library(KRLS)
library(ggplot2)
library(igraph)
library(xtable)
library(reshape2)
library(ClusterR)
library(rTensor)
library(clusterGeneration)  

############################
#       Scenario 1         #
############################

# Parameters
n=600
K=3
d=3
p=3.5*log(n)/n
q=log(n)/n
# Communities 1 and 2 are not distinguishable
Pi= diag(p,K)
Pi[1,2]=p
Pi[2,1]=p
Pi[3,2]=q
Pi[2,3]=q
Pi[3,1]=q
Pi[1,3]=q

# Communities 2 and 3 are not distinguishable
# mu = array(0,dim=c(K,K,d))
mu = array(1,dim=c(K,K,d))
mu[1,1,] = c(1,1,1)
mu[2,2,] = c(-1,-1,-1)



xp1b<-function(n,p,q,K){
  d=3
  Pi= diag(p,K)
  Pi[1,2]=p
  Pi[2,1]=p
  Pi[3,2]=q
  Pi[2,3]=q
  Pi[3,1]=q
  Pi[1,3]=q
  Z = pure_membership(n, K,rep(1/K,K))
  convZ=convertZ(Z)
  P = Z%*% Pi %*% t(Z)
  A=sample_IER(P)
  
  C = array(0, dim=c(K,K,d,d))
  for(k in 1:K){
    for(kk in 1:K){
      C[k,kk,,] = diag(1,d)
    }
  }
  
  
  
  mu = array(0,dim=c(K,K,d))
  #mu[1,1,] = c(0,0,0)
  mu[2,2,] = c(-2,2,-2)
  
  W=sample_edges_cov(A,Z,mu,C)
  W_tens=as.tensor(W)
  
  
  score = matrix(NA,nrow=1, ncol=5)
  colnames(score) = c('Spec', 'sIR-VEC','rsIR-VEC','IR-VEC','OLMF')
  
  clust_spec = clust_on_mat_eig(A,K)
  score[1,'Spec'] = NMI(clust_spec,convZ)
  
  
  clust_sIR = sIR(A,W,convertClust(clust_spec),4)
  score[1,'sIR-VEC'] = NMI(clust_sIR,convZ)
  
  #score[1,'srIR-VEC'] = NMI(sIR(A,W,pure_membership(n, K,rep(1/K,K)),7),convZ)
  
  clust_IR= IR(A,W,convertClust(clust_sIR),3)[[1]]
  score[1,'IR-VEC'] = NMI(clust_IR,convZ)
  
  x=vector("list", d)
  for(l in 1:d){
    m =max(W)
    x[[l]]= W[,,l]
  }
  
  score[1,'OLMF'] = NMI(lmfo(x,n,K),convZ)
  
  return(score)
}

#xp1b(n,p,q,3)

xp1bb = array(,dim=c(20,5))
for(i in 12:20){
  xp1bb[i,]= xp1b(n,p,q,3)
}


library(foreach)
library(doParallel)
library(parallel)


cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
resultats1bb<-foreach(t= 1:5, .combine = rbind) %dopar% {
  library(gtools)
  library(dplyr)
  library(RSpectra)
  library(mvnfast)
  library(Rcsdp)
  library(MASS)
  library(aricode)
  library(expm)
  library(KRLS)
  library(ggplot2)
  library(igraph)
  library(xtable)
  library(reshape2)
  library(ClusterR)
  library(rTensor)
  library(clusterGeneration)  
  xp1b(n,p,q,3)
  
}
parallel::stopCluster(cl)
colnames(xp1bb)=c('Spec', 'sIR-VEC','rsIR-VEC','IR-VEC','OLMF')
xp1bb= xp1bb[,-3]
write.csv(xp1bb, "xp1bb.csv",row.names=FALSE)
############################
#       Scenario 2         #
############################

C = array(0, dim=c(K,K,d,d))
for(k in 1:K){
  for(kk in 1:K){
    Cov =  genPositiveDefMat(dim = d, covMethod = "eigen", eigenvalue = c(3, 2, 1))$Sigma
    C[k,kk,,] = Cov/norm(Cov,type="2")
  }
}

####################
# Parralelization
library(mvtnorm)


rep=10


xp2<-function(n,Pi,K, mu , C){
  
  Z = pure_membership(n, K,rep(1/K,K))
  convZ=convertZ(Z)
  P = Z%*% Pi %*% t(Z)
  A=sample_IER(P)
  W=sample_edges_cov(A,Z,mu,C)
  W_tens=as.tensor(W)
  
  score = matrix(NA,nrow=1, ncol=5)
  colnames(score) = c('Spec', 'sIR-VEC','rsIR-VEC','IR-VEC','OLMF')
  
  clust_spec = clust_on_mat_eig(A,K)
  score[1,'Spec'] = NMI(clust_spec,convZ)
  
  score[1,'sIR-VEC'] = NMI(sIR(A,W,convertClust(clust_spec),3),convZ)
  
  #score[1,'srIR-VEC'] = NMI(sIR(A,W,pure_membership(n, K,rep(1/K,K)),7),convZ)
  
  score[1,'IR-VEC'] = NMI(IR(A,W,convertClust(clust_spec),3)[[1]],convZ)
  
  x=vector("list", d)
  for(l in 1:d){
    x[[l]]= W[,,l]
  }
  score[1,'OLMF'] = NMI(lmfo(x,n,K),convZ)
  
  return(score)
}



library(foreach)
library(doParallel)
library(parallel)
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
resultats<-foreach(t= 1:20, .combine = rbind) %dopar% {
  library(gtools)
  library(dplyr)
  library(RSpectra)
  library(mvnfast)
  library(Rcsdp)
  library(MASS)
  library(aricode)
  library(expm)
  library(KRLS)
  library(ggplot2)
  library(igraph)
  library(xtable)
  library(reshape2)
  library(ClusterR)
  library(rTensor)
  
  xp2(n,Pi,K,mu, C)
  
}
parallel::stopCluster(cl)

############################
#       Scenario 3         #
############################

# Increasing number of communities

# Parameters
n=1000
K=10
d=3
p=8*log(n)/n
q=4*log(n)/n
Pi= matrix(q, nrow=K, ncol=K)+diag(p-q,K)


mu = array(runif(n=K*K*d, min=-1, max=1),dim=c(K,K,d))


#Isotropic Gaussian covariates
C = array(0, dim=c(K,K,d,d))
for(k in 1:K){
  for(kk in 1:K){
    C[k,kk,,] = diag(1,d)
  }
}


xp3<-function(n,p,q,K){
  
  Pi= matrix(q, nrow=K, ncol=K)+diag(p-q,K)
  Z = pure_membership(n, K,rep(1/K,K))
  convZ=convertZ(Z)
  P = Z%*% Pi %*% t(Z)
  A=sample_IER(P)
  
  C = array(0, dim=c(K,K,d,d))
  for(k in 1:K){
    for(kk in 1:K){
      C[k,kk,,] = diag(1,d)
    }
  }
  
  
  mu = array(runif(n=K*K*d, min=-2, max=2),dim=c(K,K,d))
  W=sample_edges_cov(A,Z,mu,C)
  W_tens=as.tensor(W)
  
  
  score = matrix(NA,nrow=1, ncol=3)
  colnames(score) = c('Spec','sIR-VEC','OLMF')
  
  clust_spec = clust_on_mat_eig(A,K)
  score[1,'Spec'] = NMI(clust_spec,convZ)
  
  score[1,'sIR-VEC'] = NMI(sIR(A,W,convertClust(clust_spec),3),convZ)
  
  x=vector("list", d)
  for(l in 1:d){
    m =max(W)
    x[[l]]= W[,,l]
  }
  
  score[1,'OLMF'] = NMI(lmfo(x,n,K),convZ)
  
  return(score)
}


library(foreach)
library(doParallel)
library(parallel)


cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
resultats2<-foreach(t= 1:20, .combine = rbind) %dopar% {
  library(gtools)
  library(dplyr)
  library(RSpectra)
  library(mvnfast)
  library(Rcsdp)
  library(MASS)
  library(aricode)
  library(expm)
  library(KRLS)
  library(ggplot2)
  library(igraph)
  library(xtable)
  library(reshape2)
  library(ClusterR)
  library(rTensor)
  
  xp3(1000,p,q,2)
  
}
parallel::stopCluster(cl)

cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
resultats4<-foreach(t= 1:20, .combine = rbind) %dopar% {
  library(gtools)
  library(dplyr)
  library(RSpectra)
  library(mvnfast)
  library(Rcsdp)
  library(MASS)
  library(aricode)
  library(expm)
  library(KRLS)
  library(ggplot2)
  library(igraph)
  library(xtable)
  library(reshape2)
  library(ClusterR)
  library(rTensor)
  
  xp3(1000,p,q,4)
  
}
parallel::stopCluster(cl)


cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
resultats6<-foreach(t= 1:20, .combine = rbind) %dopar% {
  library(gtools)
  library(dplyr)
  library(RSpectra)
  library(mvnfast)
  library(Rcsdp)
  library(MASS)
  library(aricode)
  library(expm)
  library(KRLS)
  library(ggplot2)
  library(igraph)
  library(xtable)
  library(reshape2)
  library(ClusterR)
  library(rTensor)
  
  xp3(1000,p,q,6)
  
}
parallel::stopCluster(cl)


cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
resultats8<-foreach(t= 1:20, .combine = rbind) %dopar% {
  library(gtools)
  library(dplyr)
  library(RSpectra)
  library(mvnfast)
  library(Rcsdp)
  library(MASS)
  library(aricode)
  library(expm)
  library(KRLS)
  library(ggplot2)
  library(igraph)
  library(xtable)
  library(reshape2)
  library(ClusterR)
  library(rTensor)
  
  xp3(1000,p,q,8)
  
}
parallel::stopCluster(cl)


cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
resultats10<-foreach(t= 1:20, .combine = rbind) %dopar% {
  library(gtools)
  library(dplyr)
  library(RSpectra)
  library(mvnfast)
  library(Rcsdp)
  library(MASS)
  library(aricode)
  library(expm)
  library(KRLS)
  library(ggplot2)
  library(igraph)
  library(xtable)
  library(reshape2)
  library(ClusterR)
  library(rTensor)
  
  xp3(1000,p,q,10)
  
}
parallel::stopCluster(cl)

####

res_sum = matrix(NA,nrow=10, ncol=6)
colnames(res_sum) = c('Spec', 'sIR-VEC','IR-VEC','OLMF')
res_sum[1,]=apply(resultats2,2,mean)
res_sum[2,]=apply(resultats2,2,sd)
res_sum[3,]=apply(resultats4,2,mean)
res_sum[4,]=apply(resultats4,2,sd)
res_sum[5,]=apply(resultats6,2,mean)
res_sum[6,]=apply(resultats6,2,sd)
res_sum[7,]=apply(resultats8,2,mean)
res_sum[8,]=apply(resultats8,2,sd)
res_sum[9,]=apply(resultats10,2,mean)
res_sum[10,]=apply(resultats10,2,sd)

write.csv2(res_sum, "res_sum.csv",row.names = F)

############################
#  Email EU core dataset   #
############################

mail_edges = read.csv("emailEU.txt", sep=" ")
mail_inst = read.csv("emailEUID.txt", sep=" ")
Z_inst=convertClust(mail_inst$X1)
A= matrix(0, 1004,1004)
for(i in 1:25570){
  A[mail_edges[i,1],mail_edges[i,2]]=1
  A[mail_edges[i,2],mail_edges[i,1]]=1
}


indexes_inst = which(colSums(Z_inst)>50)
indexes_people = which(mail_inst[,2] %in% indexes_inst )

A=A[indexes_people,indexes_people ]
Z_inst=Z_inst[indexes_people,indexes_inst ]

### Generate edge covariates
d = 6
mu2 = array(0,dim=c(K,K,d))
for(k in 1:K){
  for(kk in k:K){
    c=runif(d,0,1)
    mu2[k,kk,]=c/sum(c)
    mu2[kk,k,]=mu2[k,kk,]
  }
}


rowMeans(rmultinom(100, d, mu2[1,1,])/d)

sample_edges_cov2<-function(A, Z, mu){
  library(MASS)
  n= dim(A)[1]
  K=dim(Z)[2]
  d=dim(mu)[3]
  convZ=convertZ(Z)
  W =array(0, dim=c(n,n,d))
  
  indices<-which(A!= 0, arr.ind = T)
  n_ind=dim(indices)[1]
  
  for(ind in 1:n_ind){
    i=indices[ind,1]
    j=indices[ind,2]
    k=convZ[i]
    kk=convZ[j]
    W[i,j,]=rowMeans(rmultinom(20, d, mu[k,kk,])/d)
    
  }
  return(W)
}
W2=sample_edges_cov2(A,Z_inst,mu2)

clust_spec = clust_on_mat_eig(A,K)
NMI(clust_spec,convertZ(Z_inst))

r=IR(A,W2,convertClust(clust_spec)  ,1)
NMI(r[[1]],convertZ(Z_inst))

library(concom)
ind_cc=concomFromMatAdj(A)$components$`1`
A2=A[ind_cc,ind_cc]
Z_inst2=Z_inst[ind_cc,]

clust_spec = clust_on_mat_eig(A2,K)
NMI(clust_spec,convertZ(Z_inst2))

W2=sample_edges_cov2(A2,Z_inst2,mu2)
r=IR(A2,W2,convertClust(clust_spec)  ,15)
NMI(r[[1]],convertZ(Z_inst2))