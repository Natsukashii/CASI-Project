# setwd("/Users/poetofquality/Nutstore Files/Nutstore/HKUST/Computer Age Statistical Inference with Applications/Project")
source("functions.R")
library(ggplot2)
set.seed(123)
m=800;n=600;r=50;itermax=500;lambda=100;tol=1e-6
sim_res <- gen_X(m=m, n=n, r=r,p=0.8,sig_e=0.1)
X_obs <- sim_res$X_obs
als_res <- ALS(X_obs, r=r, itermax=itermax, lambda=lambda, tol=tol)
sI_res <- softImpute(X_obs, itermax=itermax, lambda=lambda, tol=tol)
sI_als_res <- softImpute_ALS(X_obs, r=r, itermax=itermax, lambda=lambda, tol=tol)
sI_rank <- length((sI_res$Slambda[sI_res$Slambda>1e-10]))
sI_als_rank <- length(sI_als_res$Dsqlambda[sI_als_res$Dsqlambda>1e-10])
als_time <- als_res$time_vec[!is.na(als_res$time_vec)]
sI_time <- sI_res$time_vec[!is.na(sI_res$time_vec)]
sI_als_time <- sI_als_res$time_vec[!is.na(sI_als_res$time_vec)]
k <- min(length(als_time),length(sI_time),length(sI_als_time))
method_name <- c("ALS","softImpute","softImpute-ALS")
m <- 5 
res_df <- data.frame(x=c(als_time[m:length(als_time)],
                         sI_time[m:length(sI_time)],
                         sI_als_time[m:length(sI_als_time)]),
                     y=c(als_res$rel_obj[m:length(als_time)],
                         sI_res$rel_obj[m:length(sI_time)],
                         sI_als_res$rel_obj[m:length(sI_als_time)]),
                     method=factor(c(rep(method_name[1],length(als_time)-m+1),
                                     rep(method_name[2],length(sI_time)-m+1),
                                     rep(method_name[3],length(sI_als_time)-m+1))))

plt_df <- ggplot(res_df,aes(x=x,y=y,color=method))+
  ylab("Relative Objective (log scale)")+
  xlab("Time in seconds")+
  coord_cartesian(ylim = c(1e-6, 1e-1)) +
  scale_y_log10()+
  geom_point(lwd=1) + geom_point(aes(shape = method), size =2, show.legend = FALSE) +
  theme_bw() +
  theme(legend.position = c(.84,.86))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(2)),
        legend.text = element_text(size=rel(2.5)),
        legend.title = element_blank(),
        axis.title = element_text(hjust=0.5,size=rel(2.5)),
        axis.text = element_text(size=rel(2.5)))
plt_df


