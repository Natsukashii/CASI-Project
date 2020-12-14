gen_X <- function(m=300,n=200,r=50,p=0.8,sig_e=1){
  A <- matrix(rnorm(m*r), nrow=m, ncol=r)
  B <- matrix(rnorm(n*r), nrow=n, ncol=r)
  eps <- matrix(rnorm(m*n,sd=sig_e), nrow=m, ncol=n)
  ABt <- A%*%t(B)
  X_complete <- ABt + eps
  X_vector <- as.vector(X_complete)
  X_vector[sample(m*n, size = p*m*n)]=NA
  X_obs <- matrix(X_vector, nrow = m, ncol=n)
  return(list("X_complete"=X_complete, "X_obs"=X_obs, "A_true"=A, "B_true"=B))
}
ALS_obj <- function(X,A,B,lambda){
  0.5*sum((X-A%*%t(B))^2,na.rm = TRUE)+0.5*lambda*(norm(A,"F")^2+norm(B,"F")^2)
}
ALS <- function(X_obs, r=25, itermax=300, lambda=10, tol=1e-6){
  m <- nrow(X_obs)
  n <- ncol(X_obs)
  A_int <- matrix(rnorm(m*r), m, r)
  B_int <- matrix(rnorm(n*r), n, r)
  A <- A_int; B <- B_int;
  ind_mat <- !is.na(X_obs)
  obj_prev <- ALS_obj(X_obs,A_int,B_int,lambda)
  rel_obj <- matrix(NA, nrow = itermax)
  time_vec <- matrix(NA, nrow = itermax)
  start_time <- proc.time()[[3]]
  for (iter in 1:itermax){
    for (i in 1:m){
      nz_cols <- c(1:n)[ind_mat[i,]]
      B_submat <- B[nz_cols,]
      BtB <- t(B_submat)%*%B_submat
      A[i,] <- solve(BtB+lambda*diag(r))%*%(t(B_submat)%*%X_obs[i,nz_cols])
    }
    for (j in 1:n){
      nz_rows <- c(1:m)[ind_mat[,j]]
      A_submat <- A[nz_rows,]
      AtA <- t(A_submat)%*%A_submat
      B[j,] <- solve(AtA+lambda*diag(r))%*%(t(A_submat)%*%X_obs[nz_rows,j]) 
    }
    obj <- ALS_obj(X_obs,A,B,lambda)
    rel_obj[iter] <- log((obj_prev)/obj)
    if (rel_obj[iter]<tol){
      cat("ALS converges in", iter, "iterations")
      break
    }
    obj_prev <- obj
    time_vec[iter] <- proc.time()[[3]]-start_time
  }
  cat("\n")
  return(list("A_est"=A, "B_est"=B,"rel_obj"=rel_obj,"obj"=obj,"time_vec"=time_vec))
}
## Define soft-thredsholding operator
soft_thred <- function(x,lambda){
  y <- x-lambda
  y[y<0]=0
  return(y)
}
softImpute_obj <- function(X,M,S,lambda){
  0.5*sum((X-M)^2,na.rm = TRUE)+lambda*sum(S)
}

softImpute <- function(X_obs, itermax=500, lambda=10, tol=1e-6){
  m <- nrow(X_obs)
  n <- ncol(X_obs)
  na_mat <- is.na(X_obs)
  Mhat <- matrix(0, m, n)
  X0 <- X_obs
  X0[is.na(X0)] = 0
  obj_prev <- softImpute_obj(X_obs,Mhat,0,lambda)
  rel_obj <- matrix(NA, nrow = itermax)
  time_vec <- matrix(NA, nrow = itermax)
  start_time <- proc.time()[[3]]
  for (iter in 1:itermax){
    ## Replace the missing entries in X
    Xhat <- X0 + Mhat*na_mat
    ## Update Mhat
    svd_res <- svd(Xhat)
    S <- soft_thred(svd_res$d,lambda)
    Mhat <- sweep(svd_res$u, MARGIN = 2, S, "*") %*% t(svd_res$v)
    obj <- softImpute_obj(X_obs, Mhat, S, lambda)
    rel_obj[iter] <- log((obj_prev)/obj)
    if (rel_obj[iter]<tol){
      cat("softImpute converges in", iter, "iterations")
      break
    }
    obj_prev <- obj
    time_vec[iter] <- proc.time()[[3]]-start_time
  }
  cat("\n")
  return(list("Slambda"=S, "Mhat"=Mhat, "rel_obj"=rel_obj,"obj"=obj,"time_vec"=time_vec))
}
softImpute_ALS_obj <- function(X,ABt,Dsq,lambda){
  0.5*sum((X-ABt)^2,na.rm = TRUE)+lambda*sum(Dsq)
}
softImpute_ALS <- function(X_obs, r=25,itermax=1000,lambda=10,tol=1e-6){
  m <- nrow(X_obs)
  n <- ncol(X_obs)
  ## Intialization
  U <- matrix(rnorm(m*n), m, n)
  U <- U_prev <- svd(U, nu = r)$u
  Dsq <- Dsq_prev <- rep(1,r)
  A <- U%*%diag(sqrt(Dsq))
  V <- V_prev <- matrix(0, n, r)
  B <- V%*%diag(sqrt(Dsq))
  ABt <- A%*%t(B)
  na_mat <- is.na(X_obs)
  obj_prev <- softImpute_ALS_obj(X_obs,ABt,Dsq,lambda)
  rel_obj <- matrix(NA, nrow = itermax)
  time_vec <- matrix(NA, nrow = itermax)
  start_time <- proc.time()[[3]]
  for (iter in 1:itermax){
    ## Update B
    X0 <- X_obs - ABt
    X0[na_mat]=0
    Xstar <- X0 + ABt
    Btildet <- t(U)%*%Xstar*Dsq/(Dsq+lambda)
    Btilde_svd <- svd(t(Btildet))
    V <- Btilde_svd$u
    Dsq <- Btilde_svd$d
    U <- U%*%Btilde_svd$v
    ABt <- U%*% (Dsq*t(V))
    ## Update A
    BAt <- t(ABt)
    X0 <- t(X_obs) - BAt
    X0[is.na(X0)]=0
    Xstar <- X0 + BAt
    Atildet <- t(V)%*%Xstar*Dsq/(Dsq+lambda)
    Atilde_svd <- svd(t(Atildet))
    U <- Atilde_svd$u
    Dsq <- Atilde_svd$d
    V <- V%*%Atilde_svd$v
    ABt <- U%*% (Dsq*t(V))
    ## Check for convergence
    obj <- softImpute_ALS_obj(X_obs,ABt,Dsq,lambda)
    rel_obj[iter] <- log(obj_prev/obj)
    if (rel_obj[iter]<tol){
      cat("softImpute-ALS converges in", iter, "iterations")
      break
    }
    # cat("iteration:",iter,"relative obj:",rel_obj,"\n")
    obj_prev <- obj
    time_vec[iter] <- proc.time()[[3]]-start_time
  }
  X0 <- X_obs - ABt
  X0[na_mat]=0
  Xstar <- X0 + ABt
  Msvd <- svd(Xstar%*%V)
  U <- Msvd$u
  V <- V%*%Msvd$v
  Dsigma <- Msvd$d
  Dsigmalambda <- soft_thred(Dsigma,lambda)
  # Dsqlambda <- soft_thred(Dsq,lambda)
  cat("\n")
  res <- list("U"=U, "V"=V, "Dsq"=Dsq,"Dsigmalambda"=Dsigmalambda, "rel_obj"=rel_obj,"time_vec"=time_vec)
}

mse <- function(Omega, X_est, X_true){
  sum(X_true[Omega]-X_est[Omega])^2/sum(X_true[Omega]^2)
}
  
  
  
  