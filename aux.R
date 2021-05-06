#Quadratic programming algorithm  - WE ONLY ALLOW FOR POSITIVE DEFINITE DESIGN MATRICES
library(quadprog)
library(parallel)
library(fixest)

#Synthetic control estimator
synth_control_est <-function(pred_before, y_after, ridge = 0, Vmat = NULL)
{
  y= pred_before[,1]
  X =pred_before[,-1]
  
  if(is.null(Vmat))
  {
    Dmat = t(X)%*%X + ridge*diag(ncol(X))
    dvec = t(X)%*%y
  } else {
    Dmat = t(X)%*%Vmat%*%X + ridge*diag(ncol(X))
    dvec = t(X)%*%t(Vmat)%*%y
  }
  
  Amat = t(rbind(rep(1,ncol(X)),diag(ncol(X))))
  bvec = c(1, rep(0,ncol(X)))
  
  synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  w = synth_model$solution
  effects = y_after%*%c(1,-w)
  
  return(list("w" = w, "effects" = effects, "value_min" = (2*synth_model$value +y%*%y) ))
  
}


#Synthetic control estimator (demeaned)
synth_control_est_demean <-function(pred_before, y_after, ridge = 0, Vmat = NULL)
{
  y= pred_before[,1]
  X =cbind(1,pred_before[,-1])
  
  penalty = ridge*diag(ncol(X))
  penalty[1,1] = 0
  
  if(is.null(Vmat))
  {
    Dmat = t(X)%*%X + penalty 
    dvec = t(X)%*%y
  } else {
    Dmat = t(X)%*%Vmat%*%X + penalty
    dvec = t(X)%*%t(Vmat)%*%y
  }
  
  Amat = t(rbind(c(0,rep(1,ncol(X)-1)),diag(ncol(X))[-1,]))
  bvec = c(1, rep(0,ncol(X)-1))
  
  synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  intercept = synth_model$solution[1]
  w = synth_model$solution[-1]
  effects = -intercept + y_after%*%c(1,-w)
  
  return(list("w" = w, "intercept" = intercept, "effects" = effects, "value_min" = (2*synth_model$value +y%*%y) ))
  
}

#Computes SC ridge penalty objective function, using the CV approach of Doudchenko and Imbens (2016)
#option = "standard" or "demean" for choice of estimator
#You should minimize this function (e.g. using optimize of Base R) to find penalty
#Paralellization currently works on Unix systems only
ridge_objective <- function(lambda, pre_treat, post_treat, type = "standard", cores = detectCores()){
  val = mean(unlist(mclapply(2:ncol(pre_treat), function(j)
  {
    train = cbind(pre_treat[,j],pre_treat[,c(-1,-j)])
    eval= cbind(post_treat[,j],post_treat[,c(-1,-j)])
    
    if(type=="standard")
      effects = synth_control_est(train, eval, ridge = lambda)$effects else effects = synth_control_est_demean(train, eval, ridge = lambda)$effects
        
    mean((effects)^2)
  }, mc.cores = cores)))
  return(val)
}


#Synthetic diff in diff estimator of Athey et al (2020)
synthetic_diff_in_diff <- function(pre_treat, post_treat)
{
  mat_data = rbind(pre_treat, post_treat)
  
  xi = sd(as.vector(apply(pre_treat[,-1], 2, diff)))
  
  cross_section_ridge = synth_control_est_demean(pre_treat,pre_treat, ridge = (xi^2)*nrow(pre_treat))
  
  Time_post = colMeans(post_treat[,-1])
  Time_pre = t(pre_treat[,-1])
  Time_mat = cbind(Time_post, Time_pre)
  
  time_lasso = synth_control_est_demean(Time_mat,Time_mat, ridge = 0)
  
  weight_cross_section = c(1,cross_section_ridge$w)
  weight_cross_section[weight_cross_section<0] = 0
  
  weight_time = c(time_lasso$w, rep(1/nrow(post_treat), nrow(post_treat)))
  weight_time[weight_time<0] = 0
  
  base_final = data.frame("outcome"= as.vector(mat_data), "unit" = rep(1:ncol(mat_data), each = nrow(mat_data)),
                          "time" = rep(1:nrow(mat_data), times = ncol(mat_data)),
                          "w_i" = rep(weight_cross_section, each = nrow(mat_data)),
                          "lambda_t" =  rep(weight_time, times = ncol(mat_data)))
  
  base_final$treatment = 1*((base_final$unit==1)&(base_final$time>nrow(pre_treat)))
  
  #Diff-in-diff
  dd = feols(outcome~treatment|as.factor(unit)+as.factor(time), data = base_final)
  
  sc_did = feols(outcome~treatment|as.factor(unit)+as.factor(time), 
              data = base_final, weights=base_final$w_i*base_final$lambda_t)
  
  names(weight_cross_section) = colnames(mat_data)
  names(weight_time) = rownames(mat_data)
  
  
  return(list("sc_did" = sc_did$coefficients["treatment"],
              "did" = dd$coefficients["treatment"], "unit_weights" = weight_cross_section,
              "time_weights" = weight_time))
}