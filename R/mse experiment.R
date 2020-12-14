# setwd("/Users/poetofquality/Nutstore Files/Nutstore/HKUST/Computer Age Statistical Inference with Applications/Project")
source("functions.R")
library(ggplot2)
set.seed(123)
m=300;n=200;r=50;itermax=500;tol=1e-6
lambda_list <- seq(10,150, by=4)
mse_mat <- matrix(0, length(lambda_list), 2)
sim_res <- gen_X(m=m, n=n, r=r,p=0.8,sig_e=1)
ABt_true <- sim_res$A_true%*%t(sim_res$B_true)
X_true <- sim_res$X_complete
X_obs <- sim_res$X_obs
na_mat <- is.na(X_obs)
for (i in 1:length(lambda_list)){
  # als_res <- ALS(X_obs, r=r, itermax=itermax, lambda=lambda_list[i], tol=tol)
  # sI_res <- softImpute(X_obs, itermax=itermax, lambda=lambda_list[i], tol=tol)
  # X_als <- als_res$A_est%*%t(als_res$B_est)
  # X_sI <- sI_res$Mhat
  sI_als_res <- softImpute_ALS(X_obs, r=r, itermax=itermax, lambda=lambda_list[i], tol=tol)
  X_sI_als <- sI_als_res$U%*%diag(sI_als_res$Dsigmalambda)%*%t(sI_als_res$V)
  ## Training error
  mse_mat[i, 1] <- mse(!na_mat, X_sI_als, X_true)
  ## Test error
  mse_mat[i, 2] <- mse(na_mat, X_sI_als, ABt_true)
  cat("iteration", i,"\n")
}
rmse <- sqrt(mse_mat)
rmse_df <- data.frame(x=lambda_list,
                     y=matrix(rmse, ncol = 1),
                     method=factor(c(rep("train",length(lambda_list)),
                                     rep("test",length(lambda_list)))))
plt_df <- ggplot(rmse_df,aes(x=x,y=y,color=method))+
  ylab("RMSE")+
  xlab("Lambda")+
  geom_line(lwd=0.5) + geom_point(aes(shape = method), size =3, show.legend = FALSE) +
  theme_bw() +
  theme(legend.position = c(.84,.21))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(2)),
        legend.text = element_text(size=rel(2.5)),
        legend.title = element_blank(),
        axis.title = element_text(hjust=0.5,size=rel(2.5)),
        axis.text = element_text(size=rel(2.5)))
plt_df
rmse_df2 <- data.frame(x=rmse[,1],y=rmse[,2])
plt_df2 <- ggplot(rmse_df2,aes(x=x,y=y,color="red"))+
  ylab("Test RMSE")+
  xlab("Train RMSE")+
  geom_line(lwd=0.5) + geom_point(size =3, show.legend = FALSE) +
  theme_bw() +
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5,size=rel(2)),
        legend.text = element_text(size=rel(2.5)),
        legend.title = element_blank(),
        axis.title = element_text(hjust=0.5,size=rel(2.5)),
        axis.text = element_text(size=rel(2.5)))
plt_df2
