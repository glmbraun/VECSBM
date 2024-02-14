####################### 
# Auxiliary functions #
####################### 

# Norms

norm_vec2 <- function(x) sum(x^2)

norm_vec2_gamma<-function(vec, Gamma){
  d=dim(vec)[1]
  R=ginv(Gamma)
  return(t(vec)%*%R%*%vec-0.5*det(Gamma))
}



#Convert a membership matrix to a vector 
convertZ<-function(Z){
  n<-dim(Z)[1]
  return(sapply(1:n, function(i) {match(1, Z[i,])}))
}

# Convert a vector of labels to a membership matrix
convertClust<-function(clust){
  n<-length(clust)
  k<-length(unique(clust))
  Z<-matrix(0,nrow=n,ncol=k)
  for(i in 1:n){Z[i,clust[i]]<-1}
  return(Z)
}


hard_threshold = function(M, delta){
  # delta : threshold value
  M[(M) <= delta] = 0
  return(M)
}


normalize_row<-function(M){
  return(diag(1/rowSums(abs(M)))%*%M)
}


laplacian <- function(A, tau = 0) {
  Atau = A+tau
  degree = apply(Atau,1,sum)
  diag(1/sqrt(degree))%*%Atau%*%diag(1/sqrt(degree))
}

laplacian_reg <- function(A, tau = 0) {
  degree = apply(A,1,sum) + tau
  diag(1/sqrt(degree))%*%A%*%diag(1/sqrt(degree))
}


clust_on_mat_eig<-function(A,k){
  return(kmeans(svds(A,k)$u,k)$cluster)
}

spectral.norm <- function( x )
{
  if ( !is.numeric( x ) ) {
    stop( "argument x is not numeric" )
  }
  if ( is.vector( x ) ) {
    return( sqrt( sum( x * x ) ) )
  }
  if ( !is.matrix( x ) ) {
    return( "argument x is not a matrix" )
  }
  A <- t(x) %*% x
  eigenA <- eigen( A )
  lambdaA <- eigenA$values
  maxLambdaA <- lambdaA[1]
  if ( maxLambdaA < 0 ) {
    stop( "t(x) %*% x is negative definite" )
  }
  return( sqrt( maxLambdaA ) )
}

############################
##### Graph generation #####
############################

# Generate an inhomogeneous Erdös-Renyi graph from a matrix P
sample_IER<-function(P){
  A = apply(P,MARGIN = c(1,2),function(u) rbinom(1,1,prob = min(u,1)))  
  A[lower.tri(A)] = 0
  A =  A + t(A)
  diag(A) = 0
  return(A)
}

# Generate a membership matrix with no overlap
pure_membership<-function(n, K, alpha){
  return(t(rmultinom(n,1,alpha)))
}

sample_edges_cov<-function(A, Z, mu, C){
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
    W[i,j,]=mvrnorm(1,mu=mu[k,kk,],Sigma=C[k,kk,,])
    
  }
  return(W)
}

###########################
##### Main algorithms #####
###########################

sIR<-function(A,W,Z_init,iter){
  library(rTensor)
  n= dim(A)[1]
  K=dim(Z_init)[2]
  d=dim(W)[3]
  convZ=convertZ(Z_init)
  Z=Z_init
  W_tens=as.tensor(W)
  
  for(t in 1:iter){
    # Estimate model parameters
    Z=convertClust(convZ)
    Zn= Z%*% diag(1/colSums(Z))
    Pi = t(Zn)%*% A %*% Zn
    
    mu = ttm(W_tens,t(Z),m=1)
    mu = ttm(mu, t(Z),m=2)
    
    
    for(k in 1:K){
      for(kk in 1:K){
        index_k= which(convZ==k)
        index_kk= which(convZ==kk)
        mu[k,kk,]=mu[k,kk,]/sum(A[index_k,index_kk])
      }
    }
    
    
    # Refinement of the partition
    for(i in 1:n){
      map = 1:K
      logPi=log(Pi)
      logmPi=log(1-Pi)
      for(k in 1:K){
        
        #f<-function(j){A[i,j]*log(Pi[k, convZ[j]])+(1-A[i,j])*log(1-Pi[k, convZ[j]])-A[i,j]*norm_vec2(W[i,j,]-vec(mu[k,convZ[j],]))}
        f<-function(j){A[i,j]*logPi[k, convZ[j]]+(1-A[i,j])*logmPi[k, convZ[j]]-A[i,j]*norm_vec2(W[i,j,]-vec(mu[k,convZ[j],]))}
        j_values = 1:n
        j_values= j_values[j_values!= i]
        map[k] =sapply(j_values,f)%>%sum()
      }
      convZ[i] = which.max(map)
    }
  }
  return(convZ) 
}


IR<-function(A,W,Z_init,iter){
  library(rTensor)
  n= dim(A)[1]
  K=dim(Z_init)[2]
  d=dim(W)[3]
  convZ=convertZ(Z_init)
  Z=Z_init
  W_tens=as.tensor(W)
  
  for(t in 1:iter){
    # Estimate model parameters
    Z=convertClust(convZ)
    Zn= Z%*% diag(1/colSums(Z))
    Pi = t(Zn)%*% A %*% Zn
    
    mu = ttm(W_tens,t(Z),m=1)
    mu = ttm(mu, t(Z),m=2)
    
    Gamma = array(0,dim=c(K,K,d,d))
    
    for(k in 1:K){
      for(kk in 1:K){
        index_k= which(convZ==k)
        index_kk= which(convZ==kk)
        mu[k,kk,]=mu[k,kk,]/sum(A[index_k,index_kk])
      }
    }
    
    for(k in 1:K){
      for(kk in 1:K){
        index_k= which(convZ==k)
        index_kk= which(convZ==kk)
        v= array(0,dim=c(d,d))
        for(i in index_k){
          for(j in index_kk){
            w=rep(0,d)
            for(dd in 1:d){
              #print(W[i,j,dd])
              w[dd] = W[i,j,dd]-mu@data[k,kk,dd]}
            v= v+A[i,j]*w%*%t(w)
          }
        }
        
        Gamma[k,kk,,]=v/sum(A[index_k,index_kk])
      }
    }
    
    # Refinement of the partition
    for(i in 1:n){
      map = 1:K
      logPi=log(Pi)
      logmPi=log(1-Pi)
      for(k in 1:K){
        
        #f<-function(j){A[i,j]*log(Pi[k, convZ[j]])+(1-A[i,j])*log(1-Pi[k, convZ[j]])-A[i,j]*norm_vec2(W[i,j,]-vec(mu[k,convZ[j],]))}
        f<-function(j){A[i,j]*logPi[k, convZ[j]]+(1-A[i,j])*logmPi[k, convZ[j]]-A[i,j]*norm_vec2_gamma(W[i,j,]-mu@data[k,convZ[j],], Gamma[k,convZ[j],,])}
        j_values = 1:n
        j_values= j_values[j_values!= i]
        map[k] =sapply(j_values,f)%>%sum()
      }
      convZ[i] = which.max(map)
    }
  }
  return(list(convZ,Gamma)) 
}


################# Implementation of orthogonal LMF used in Paul and Chen, Annals of Statistics 2020.  ##################

### Create the list of Laplacian matrices

laplacian<-function(x,weights=rep(1,length(x))){
  M=length(x)
  laplist<-lapply(1:M,function(m){
    d=rowSums(x[[m]])
    d=d+mean(d)
    deg<-diag(1/sqrt(d))
    lap=deg%*%x[[m]]%*%deg
    return(lap)
  })
  return(laplist)
}


## The objective function 

lmffunctiono<-function(param,laplist,n,k){
  M=length(laplist)
  ustar<-matrix(param[1:(n*k)],n,k)
  lambda<-lapply(1:M,function(m){return(matrix(param[(n*k+(m-1)*k^2+1):(n*k+m*k^2)],k,k))})
  objloop<- sum(unlist(lapply(1:M,function(m){
    specobj<-norm(laplist[[m]]-ustar%*%lambda[[m]]%*%t(ustar),type="F")^2
    return(specobj)
  })))
  obj=objloop
  return(obj)
}


##  The gradients

lmfdero<-function(param,laplist,n,k){
  M=length(laplist)
  ustar<-matrix(param[1:(n*k)],n,k)
  lambda<-lapply(1:M,function(m){return(matrix(param[(n*k+(m-1)*k^2+1):(n*k+m*k^2)],k,k))})
  derlist1<-lapply(1:M,function(m){
    specobj= -(diag(n)-ustar%*%t(ustar))%*%laplist[[m]]%*%ustar%*%lambda[[m]]
    return(specobj)
  })
  derlist2<-lapply(1:M,function(m){
    specobj= -t(ustar)%*%(laplist[[m]]-ustar%*%lambda[[m]]%*%t(ustar))%*%ustar
    return(specobj)
  })
  der1<-Reduce("+",derlist1)
  der2<-unlist(derlist2)
  return(c(as.vector(der1),as.vector(der2)))
}


## The main function with BFGS optimization

lmfo<-function(x,n,k){
  M=length(x)
  #laplist<-laplacian(x)
  laplist<-x
  # Initialize with mean laplacian
  lapmean<-Reduce("+",laplist)
  spectra<-eigen(lapmean)
  ustar<-spectra$vectors[,1:k]
  lambda<-lapply(1:M,function(m){return(diag(spectra$values[1:k]))})
  param<-c(as.vector(ustar),as.vector(unlist(lambda)))
  optimized <-optim(par=param,fn=lmffunctiono,gr=lmfdero,method="BFGS",control=list(reltol=0.0001,maxit=200),laplist=laplist,n=n,k=k)
  param<-optimized$par
  
  ustar<-matrix(param[1:(n*k)],n,k)
  lambda<-lapply(1:M,function(m){return(matrix(param[(n*k+(m-1)*k^2+1):(n*k+m*k^2)],k,k))})
  
  specstar<-kmeans(ustar,k)
  specclus<-specstar$cluster
  return(specclus)
  
}
