
# TMB version 1.9.16
# stringr version 1.5.1
# ggplot2 version 3.5.1
# lattice version 0.22-6

library(TMB)
library(stringr)
library(ggplot2)
library(lattice)

load("./Data/pop_dy_inputs.RData")

setwd("./src/")

# compile the cpp file for the population dynamics
dyn.unload("state_space_pop_dy")
compile("state_space_pop_dy.cpp")
dyn.load("state_space_pop_dy")

# set-up model settings (these will later turn on or off how data and parameters are defined)
model_run<-list(
  rec_covar_switch = 0, # Covariate option for recruitment: 0 no enviro #1 climate (not possible with this version of the code)
  rec_nll_switch = 2, # Likelihood function for recruitment: 0 is dnorm, 1 is ar1, 2 is rw
  M_RE_switch = 1, # Model M as a random effect: 0 is off, 1 is on
  M_nll_switch = 2, # Likelihood function for M: 0 means dnorm, 1 is ar1, 2 is rw (different by age)
  M_covar_switch = 1, # Covariate option for M: 0 no enviro #1 condition index
  YE_switch = 0 # Year effects for the survey index: 0 is no YE, 1 is AR(1), 2 is dnorm
)

covar_name<-"cond"

# define the data object for the model
tmb.data<-list(
  weight = model_inputs$stock_weights,
  mat = model_inputs$new_mat,
  midy_weight = model_inputs$new_midy_wgt,
  C = model_inputs$updated_catch[1:60,]/1000,  
  index = log(model_inputs$index_df$index/100),
  iyear = model_inputs$index_df$iyear,
  iage = model_inputs$index_df$iage, 
  isurvey = model_inputs$index_df$isurvey,
  log_landings = log((c(model_inputs$landings$Total[1:25], model_inputs$landings$STACFISa[26:60])/1000)),
  fs = model_inputs$fs,
  A = 15,
  Y = 60,
  is_year = model_inputs$index_df$is_year,
  Ns = length(unique(model_inputs$index_df$is_year)),
  isurvey1 = c(rep(0, 29), rep(1,31)),
  isd = model_inputs$index_df$isd,
  ssb_start = rep(100,60),
  std_landings = 0.1,
  rec_covar_switch = model_run$rec_covar_switch, 
  rec_nll_switch = model_run$rec_nll_switch, 
  M_RE_switch = model_run$M_RE_switch,
  M_nll_switch = model_run$M_nll_switch, 
  M_covar_switch = model_run$M_covar_switch,
  YE_switch = model_run$YE_switch, 
  mean_midy_weight = model_inputs$mean_midy_weight,
  crl = model_inputs$crl,
  log_condM_dat = model_inputs$condM_dat,
  log_sdcond = log(0.2)
)

# define parameter starting values for the model
parameters<-list( 
  log_No = rep(9.2,14),  
  log_Rec_mean = c(log(5),log(3)), 
  log_std_log_R = rep(log(2),2),
  m_q = model_inputs$updated_mq,
  log_cv_index = rep(0.5,30),
  logit_ar_logF_year = log(0.99/0.01),        
  logit_ar_logF_age = log(0.99/0.01), 
  log_M_dev =  rep(log(0.01), 60),
  log_sd_M = log(0.01),
  log_Rec_dev = rep(0.5,60),   
  log_F = matrix(log(0.1),nrow=60,ncol=11),
  log_alpha=log(50),
  log_beta=log(0.015),
  log_init_rec = log(20),
  log_intercept_R=rep(0,1),
  log_phi_R=rep(0,1),
  intercept_M=rep(0,1),
  phi_M=rep(0,1),
  logit_log_R = log(0.99/0.01),
  logit_ar_logM_year = log(0.99/0.01),        
  logit_ar_logM_age = log(0.99/0.01),
  logit_ar_iye_year_spr = log(0.99/0.01),
  logit_ar_iye_year_fall = log(0.99/0.01),
  log_std_iye_spr = log(0.1),
  log_std_iye_fall = log(0.1),
  iye_spr = rep(0, 60),
  iye_fall = rep(0, 60),
  log_std_crl = matrix(log(0.2),nrow=60,ncol=10),
  logit_ar_crl = c(2.2,-10),
  log_sdrwM = log(sd(tmb.data$log_condM_dat)),
  log_condM = tmb.data$log_condM_dat,
  init_log_condM = log(0.5)
)

# define maps to estimate similar parameters for some ages/years
log_std_crl_map<-factor(c(rep("5-6",60*2), rep("7-11",60*5), rep("12-14",60*3)))

mapmq_mat<-matrix(model_inputs$m_q, nrow=56, ncol=14)
new_mapmq_mat<-rbind(mapmq_mat[1:26,],mapmq_mat[25:26,],mapmq_mat[27:56,],mapmq_mat[55:56,])

map_stdcrl_mat<-matrix(model_inputs$log_std_crl,nrow=58, ncol=10)
new_mapstdcrl_mat<-rbind(map_stdcrl_mat,map_stdcrl_mat[57:58,])

map <- list(
  m_q = factor(new_mapmq_mat),
  log_cv_index = model_inputs$log_cv_index,
  log_std_crl = factor(new_mapstdcrl_mat),
  logit_ar_logM_age = factor(NA),
  intercept_M=factor(NA),
  log_alpha=factor(NA),
  log_beta=factor(NA)
)

# switches to turn off estimation of some parameters, depending on the model format
if(model_run$rec_covar_switch==0){
  map <- append(map, list(  log_intercept_R=factor(NA),
                            log_phi_R=factor(NA) )  )
}

if(model_run$M_covar_switch==0){
  map <- append(map, list(  intercept_M=factor(NA),
                            phi_M=factor(NA) )  )
}

if(model_run$rec_nll_switch==0|model_run$rec_nll_switch==2){
  map <- append(map, list(  logit_log_R = factor(NA) )  )
}

if(model_run$M_nll_switch==0|model_run$M_nll_switch==2){
  map <- append(map, list( logit_ar_logM_year = factor(NA),       
                           logit_ar_logM_age = factor(NA) )  )
}

if(model_run$M_RE_switch==0){
  map <- append(map, list(log_M_dev =  factor(rep(NA, 60)),   
                          log_sd_M = as.factor(c(rep(NA,1)))  )  )
}

if(model_run$YE_switch==0){
  map <- append(map, list(   logit_ar_iye_year_spr = factor(NA),
                             logit_ar_iye_year_fall = factor(NA),
                             log_std_iye_spr = factor(NA),
                             log_std_iye_fall = factor(NA),
                             iye_spr = as.factor(c(rep(factor(NA),60))),
                             iye_fall = as.factor(c(rep(factor(NA),60))) )  )
}

if(model_run$YE_switch==2){
  map <- append(map, list(   logit_ar_iye_year_spr = factor(NA),
                             logit_ar_iye_year_fall = factor(NA))  )
}

if(model_run$M_covar_switch<1){
  
  map <- append(map, list(  intercept_M=factor(NA)))
  
  map<-append(map, list(phi_M = factor(NA)))
}  

map<-append(map, list(log_condM = factor(c(NA, seq(from=2, to=30, by=1)))))

# define random effects
rnames<-c("log_Rec_dev","log_F", "log_condM")

#if not estimating M, but estimating YE
if(model_run$YE_switch>0){
  rnames<-append(rnames, c("iye_spr", "iye_fall"))
}

if(model_run$M_RE_switch>0){
  rnames<-append(rnames, c("log_M_dev"))
}


# prepare the model
obj <- MakeADFun(tmb.data,parameters,map,
                 random=rnames,
                 DLL="state_space_pop_dy",inner.control=list(maxit=500,trace=F))

# estimate parameters
opt<-nlminb(obj$par,obj$fn,obj$gr,
            control = list(trace=10,eval.max=2000,iter.max=1000)) 

obj$gr()

# check gradient
max(obj$gr())

# save output
rep = obj$report()

AIC=2*opt$objective+ 2*length(opt$par) #AIC
BIC=2*opt$objective + log(length(tmb.data$iage))*length(obj$par) #BIC
AIC
BIC

# save sd estimates
sdrep<-sdreport(obj)
