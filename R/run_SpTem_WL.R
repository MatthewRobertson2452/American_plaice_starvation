
# Matrix version 1.7-1
# ggplot2 version 3.5.1
# patchwork version 1.3.0
# TMB version 1.9.16
# reshape2 version 1.4.4
# lattice version 0.22-6
# latticeExtra version 0.6-30
# xtable version 1.8-4
# corrplot version 0.95

library(Matrix)
library(ggplot2)
library(patchwork)
library(TMB)  
library(reshape2) 
library(lattice)
library(latticeExtra)
library(xtable)             
library(corrplot)

load("./Data/distance_mat_strat")
load("./Data/strat_area")
load("./Data/lw_mod_inputs.RData")


##################################################

setwd("./src/")

TMBmodel <- "SpTem_WL"

# Build TMB object
TMB::compile( file=paste0(TMBmodel,".cpp"),
              flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign") 

dyn.load( dynlib(TMBmodel) )

#########################################################
nyears = 30

prec_mat = distance_mat_strat

rownum = sum(prec_mat>0)/2
row_non0 = rep(0,rownum)
col_non0 = rep(0,rownum)

i_now = 1
for( i in 1:(length(prec_mat[1,])-1) ){
  for( j in (i+1):length(prec_mat[1,]) ){
    if(prec_mat[i,j]>0){
      row_non0[i_now] = i
      col_non0[i_now] = j
      i_now = i_now + 1
    }#if
  }#j
}#i

strat_area0 <- c()
for(si in 1:ncol(prec_mat)){
  strat_area0[si] <- strat_area$norm_area[strat_area$strat%in%(colnames(prec_mat)[si])]
}
strat_area <- strat_area0

#=====================================================
SIGrhoPOW<-matrix(NA,
                  ncol=length(unique(model_inputs$JD_CAT_reff$JDX_cat)),
                  nrow=length(unique(model_inputs$JD_CAT_reff$JDX_cat)))

for(i in 1:ncol(SIGrhoPOW)){
  for(j in 1:nrow(SIGrhoPOW)){  
    dist<-abs(i-j)
    if(dist <= (ncol(SIGrhoPOW)-1)/2){
      SIGrhoPOW[i,j]<- dist;
    }else{
      SIGrhoPOW[i,j]<- ncol(SIGrhoPOW)-dist;
    }
  }
}
#=====================================================



NUMfn<-function(x){return(ifelse(floor(x/10)==0,paste0("0",x,"_"),paste0(x,"_")))}

fit_fun <- function(use_X,Mnum,BCtrue=FALSE,INITs){

  
  ad_obj<-NULL; Fitted_AIC<-NULL;  Fitted_BIC<-NULL;
  obj<-NULL; opt<-NULL; rep.nowt<-NULL; WL_A_B<-NULL;  Ew_y<-NULL; 
  
  #########################################################
  tmb.data = list(
    use_Ag     = ifelse(use_X[ 1]=="+",1,0), 
    use_At     = ifelse(use_X[ 2]=="+",1,0), 
    use_Al     = ifelse(use_X[ 3]=="+",1,0), 
    use_Agt    = ifelse(use_X[ 4]=="+",1,0), 
    use_Agl    = ifelse(use_X[ 5]=="+",1,0), 
    use_Atl    = ifelse(use_X[ 6]=="+",1,0), 
    use_Bc     = ifelse(use_X[ 7]=="+",1,0), 
    use_JL     = ifelse(use_X[ 8]=="+",1,0), 
    K_crit     = -0.18,
    nyears     = nyears,
    nstrat     = length(prec_mat[1,]),
    nlen       = length(min(model_inputs$new_data_length_weight$length):max(model_inputs$new_data_length_weight$length)),
    nlen_cat   = length(min(model_inputs$new_data_length_weight$ilen_cat):max(model_inputs$new_data_length_weight$ilen_cat)),
    log_wei    = log(model_inputs$new_data_length_weight$weight),
    log_len    = log(model_inputs$new_data_length_weight$length),# - mean( log(data_length_weight$length) ),
    min_len    = min(model_inputs$new_data_length_weight$length),# - mean( log(data_length_weight$length) ),
    min_JD     = min(model_inputs$new_data_length_weight$julian_date),
    iyears     = model_inputs$new_data_length_weight$year - min(model_inputs$new_data_length_weight$year),
    istrat     = model_inputs$new_data_length_weight$istrat - 1,
    ilen       = model_inputs$new_data_length_weight$ilen - 1,
    ilen_cat   = model_inputs$new_data_length_weight$ilen_cat - 1,
    iJD_cat    = model_inputs$new_data_length_weight$iJD - 1,
    spatiocorr = prec_mat,
    diag_Mat   = as(diag(rep(1,length(prec_mat[1,]))),"dgTMatrix"),
    SIGrhoPOW  = SIGrhoPOW,
    nofJD      = length(unique(model_inputs$JD_CAT_reff$JDX_cat)), 
    row_non0   = row_non0-1,
    col_non0   = col_non0-1,
    area           = strat_area,
    len_v_sub      = seq(20,60,by=10),
    JD_v_sub       = c(4,10) 
  )
  
  for(i in 1:length(tmb.data$len_v_sub)){
    tmb.data$len_cat_v_sub[i] <- model_inputs$LEN_CAT_reff$LenghtX_cat[model_inputs$LEN_CAT_reff$LenghtX==tmb.data$len_v_sub[i]] -1
  }
  for(i in 1:length(tmb.data$JD_v_sub)){
    tmb.data$JD_cat_v_sub[i] <- model_inputs$JD_CAT_reff$JDX_cat[model_inputs$JD_CAT_reff$JDX==tmb.data$JD_v_sub[i]] -1
  }
  
  
  ModelName = paste0(c(ifelse(tmb.data$use_Ag==1,    "_Ag","_xx"),
                       ifelse(tmb.data$use_At==1,    "_At","_xx"),
                       ifelse(tmb.data$use_Al==1,    "_Al","_xx"),
                       ifelse(tmb.data$use_Agt==1,  "_Agt","_xxx"),
                       ifelse(tmb.data$use_Agl==1,  "_Agl","_xxx"),
                       ifelse(tmb.data$use_Atl==1,  "_Atl","_xxx"),
                       ifelse(tmb.data$use_Bc==1,    "_Bc","_xx"),
                       ifelse(tmb.data$use_JL==1,    "_JL","_xx")),collapse = "")
  
  print("#####################################################################")
  print(ModelName)
  print("#####################################################################")
  
  write(paste0(ModelName),file="track.txt", append = T)
  
  rname_temp = na.omit(c(ifelse(tmb.data$use_Ag==1,  "Delta_Ag" ,NA),
                         ifelse(tmb.data$use_At==1,  "Delta_At" ,NA),
                         ifelse(tmb.data$use_Al==1,  "Delta_Al" ,NA),
                         ifelse(tmb.data$use_Agt==1, "Delta_Agt",NA),
                         ifelse(tmb.data$use_Agl==1, "Delta_Agl",NA),
                         ifelse(tmb.data$use_Atl==1, "Delta_Atl",NA),
                         ifelse(tmb.data$use_JL==1,  "Delta_JL" ,NA)))
  rname<-NULL
  if(length(rname_temp)!=0)rname<-rname_temp
  #-------------------------------------------------------------------------------------
  map_tmb <- list(  )
  
  if(tmb.data$use_Bc==0){
    map_tmb <- append(map_tmb, list( log_beta  = factor(NA) )  )
  }
  if(tmb.data$use_Ag==0){
    map_tmb <- append(map_tmb,
                      list(
                        log_q_Ag      = factor(NA),
                        log_omeg_Ag   = factor(NA),
                        Delta_Ag      = factor(rep(NA,tmb.data$nstrat))
                      )  )
  }
  if(tmb.data$use_At==0){
    map_tmb <- append(map_tmb,
                      list(
                        logit_phi_At  = factor(NA),
                        log_sigma_At  = factor(NA),
                        Delta_At      = factor(rep(NA,tmb.data$nyears))
                      )  )
  }
  if(tmb.data$use_Al==0){
    map_tmb <- append(map_tmb,
                      list(
                        logit_phi_Al  = factor(NA),
                        log_sigma_Al  = factor(NA),
                        Delta_Al      = factor(rep(NA,tmb.data$nlen_cat))
                      )  )
  }
  if(tmb.data$use_Agt==0){
    map_tmb <- append(map_tmb,
                      list(
                        log_q_Agt     = factor(NA),
                        log_omeg_Agt  = factor(NA),
                        logit_phi_Agt = factor(NA), 
                        Delta_Agt     = factor(matrix(NA,tmb.data$nstrat,tmb.data$nyears))
                      )  )
  }
  if(tmb.data$use_Agl==0){
    map_tmb <- append(map_tmb,
                      list(
                        log_q_Agl     = factor(NA),
                        log_omeg_Agl  = factor(NA),
                        logit_phi_Agl = factor(NA), 
                        Delta_Agl     = factor(matrix(NA,tmb.data$nstrat,tmb.data$nlen_cat))
                      )  )
  }
  if(tmb.data$use_Atl==0){
    map_tmb <- append(map_tmb,
                      list(
                        logit_phi_Atl_t  = factor(NA),
                        logit_phi_Atl_l  = factor(NA),
                        log_sigma_Atl    = factor(NA),
                        Delta_Atl        = factor(matrix(NA,tmb.data$nyears,tmb.data$nlen_cat))
                      )  )
  }
  if(tmb.data$use_JL==0){
    map_tmb <- append(map_tmb,
                      list(
                        logit_phi_JL_J  = factor(NA),
                        logit_phi_JL_L  = factor(NA),
                        log_sigma_JL    = factor(NA),
                        Delta_JL        = factor(matrix(NA,
                                                        tmb.data$nofJD,
                                                        diff(range(tmb.data$ilen_cat))+1))
                      )  )
  }
  if(tmb.data$use_Agt+tmb.data$use_Agl+tmb.data$use_Ag==0){
    map_tmb <- append(map_tmb, list( tao = factor(NA) )  )
  }
  #-------------------------------------------------------------------------------------
  
  
  map_tmb <- append(map_tmb, list( logit_phi_JL_J = factor(NA) )  )
  
  # ====================================================================
  FITfn<-function(phi_0_tem, a_w_0_tem, tao_0_tem, sig_0_tem){
    
    
    parameters_tem <- list(
      a_w             = a_w_0_tem,
      log_c_w         = log(0.2),
      d_w             = -4,
      e_w             = 0,
      log_beta        = log(3),
      log_q_Ag        = log(10),
      log_q_Agt       = log(10),
      log_q_Agl       = log(10),
      log_omeg_Ag     = log(0.1),
      log_omeg_Agt    = log(0.1),
      log_omeg_Agl    = log(0.1),
      logit_phi_At    = log( phi_0_tem / (1-phi_0_tem) ), # (0,1) -> (-Inf,Inf)
      logit_phi_Al    = log( phi_0_tem / (1-phi_0_tem) ), # (0,1) -> (-Inf,Inf)
      logit_phi_Agt   = log( phi_0_tem / (1-phi_0_tem) ), # (0,1) -> (-Inf,Inf)
      logit_phi_Agl   = log( phi_0_tem / (1-phi_0_tem) ), # (0,1) -> (-Inf,Inf)
      logit_phi_Atl_t = log( phi_0_tem / (1-phi_0_tem) ), # (0,1) -> (-Inf,Inf)
      logit_phi_Atl_l = log( phi_0_tem / (1-phi_0_tem) ), # (0,1) -> (-Inf,Inf)
      logit_phi_JL_J  = log(0.99/(1-0.99)),#log( phi_0_tem / (1-phi_0_tem) ), # (0,1) -> (-Inf,Inf)
      logit_phi_JL_L  = log( phi_0_tem / (1-phi_0_tem) ), # (0,1) -> (-Inf,Inf)
      log_sigma_At    = log(sig_0_tem),
      log_sigma_Al    = log(sig_0_tem),
      log_sigma_Atl   = log(sig_0_tem),
      log_sigma_JL    = log(sig_0_tem),
      tao             = tao_0_tem,
      ################ random effects ################
      Delta_Ag   = rep(0,tmb.data$nstrat),
      Delta_At   = rep(0,tmb.data$nyears),
      Delta_Al   = rep(0,diff(range(tmb.data$ilen_cat))+1),
      Delta_Agt  = matrix(0,tmb.data$nstrat,tmb.data$nyears),
      Delta_Agl  = matrix(0,tmb.data$nstrat,diff(range(tmb.data$ilen_cat))+1),
      Delta_Atl  = matrix(0,tmb.data$nyears,diff(range(tmb.data$ilen_cat))+1),
      Delta_JL   = matrix(0,
                          tmb.data$nofJD,
                          diff(range(tmb.data$ilen_cat))+1)
    ) 
    obj_tem<-NULL
    opt_tem<-NULL
    error_r_tem=try( ( obj_tem <- MakeADFun(tmb.data,
                                            parameters_tem,
                                            random=rname,
                                            map=map_tmb,
                                            random.start=expression(last.par[random]), 
                                            DLL=TMBmodel, 
                                            inner.control=list(maxit=10000,trace=F)) ), 
                     silent=TRUE )
    # ====================================================================
    OUT<-NULL
    no_of_runs_tem<-NULL
    
    tryCatch({
      
      obj_tem$env$tracemgc <- FALSE
      
      obj_tem$gr(obj_tem$par)
      obj_tem$fn(obj_tem$par)
      opt_tem<-nlminb(obj_tem$par,obj_tem$fn,obj_tem$gr,
                      control = list(trace=1,iter.max=10000,
                                     eval.max=10000,sing.tol=1e-20))
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # Run until it converges
      no_of_runs_tem<-1
      NANvals<-1
      if(opt_tem$convergence==0 & is.finite(opt_tem$objective)){
        ad_obj=sdreport(obj_tem)
        NANvals<-sum(is.na(summary(ad_obj, "fixed")[,2]))
      }
      
      while(!(opt_tem$convergence==0 & is.finite(opt_tem$objective) & NANvals==0) & no_of_runs_tem<2 ){
        print(paste("========================",no_of_runs_tem,"========================"))
        no_of_runs_tem <- no_of_runs_tem + 1
        
        opt_tem<-nlminb(opt_tem$par,obj_tem$fn,obj_tem$gr,#lower=lower,upper=upper,
                        control = list(trace=1,iter.max=10000,
                                       eval.max=10000,sing.tol=1e-20))
        
        
        if(opt_tem$convergence==0 & is.finite(opt_tem$objective)){
          ad_obj=sdreport(obj_tem)
          NANvals<-sum(is.na(summary(ad_obj, "fixed")[,2]))
        }
      }
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      OUT<-list(opt=opt_tem,
                obj=obj_tem,
                parameters=parameters_tem,
                phi_0=phi_0_tem,
                no_of_runs=no_of_runs_tem,
                error_r=error_r_tem,
                NANvals=NANvals)
      
      if(!is.null(opt_tem)){
        if(opt_tem$convergence!=0){OUT<-NULL}
        if(!is.finite(opt_tem$objective)){OUT<-NULL}
      }
      
    },error=function(e){}); 
    
    return(OUT)
  }
  # ====================================================================

  # #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  LOOP_stop <- FALSE
  NLLvals<-rep(NA,length=nrow(INITs))
  covNANs<-rep(NA,length=nrow(INITs))
  BestFITind <- NULL
  
  fi <- 1
  
  while( (fi<=nrow(INITs)) & (!LOOP_stop) ){
    TEM_fit<-NULL
    tryCatch({ TEM_fit<-FITfn(INITs$phi_0s[fi],
                              INITs$a_w_0s[fi],
                              INITs$tao_0s[fi],
                              INITs$sig_0s[fi]) },error=function(e){});
    if(!is.null(TEM_fit)){
      NLLvals[fi]<-TEM_fit$opt$objective
      covNANs[fi]<-TEM_fit$NANvals
      if(covNANs[fi]==0){
        BestFITind <- fi
        LOOP_stop <- TRUE
        FiTeD <- TEM_fit
      }
    }
    
    print(data.frame(cbind(INITs,NLLvals,covNANs))[1:fi,])
    write(paste0(fi," : ",NLLvals[fi]," - ",covNANs[fi]," - covNANs"),file="track.txt", append = T)
    
    fi <- fi + 1
  }
  
  if(is.null(BestFITind)){
    BestFITind <- which(NLLvals==min(na.omit(NLLvals[covNANs==min(covNANs,na.rm = T)])))[1]
    
    FiTeD<-NULL
    tryCatch({ FiTeD<-FITfn(INITs$phi_0s[BestFITind],
                            INITs$a_w_0s[BestFITind],
                            INITs$tao_0s[BestFITind],
                            INITs$sig_0s[BestFITind]) },error=function(e){});
  }
  
  obj        <- FiTeD$obj
  opt        <- FiTeD$opt
  parameters <- FiTeD$parameters
  phi_0      <- FiTeD$phi_0
  no_of_runs <- FiTeD$no_of_runs
  error_r    <- FiTeD$error_r
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  rep.nowt = obj$report()
  
  if(!BCtrue){ ad_obj=sdreport(obj) }
  if(BCtrue){ad_obj=sdreport(obj,bias.correct=TRUE,bias.correct.control=list(sd = TRUE))}
  
  ad_obj$cov <- NULL
  
  Fitted_AIC<-2*opt$objective + 2*length(obj$par)
  Fitted_BIC<-2*opt$objective + log(length(tmb.data$log_wei))*length(obj$par)
  
  #===========================================================================================
  if(!dir.exists("check"))dir.create("check")
  FNme <- file(paste0("check/M",Mnum,ModelName,gsub(" ","",gsub("-","",gsub(":","",Sys.time()))),".txt"))
  sink(FNme, append=TRUE)
  print(str(tmb.data))
  print(paste("phi_0 =",phi_0))
  print(str(parameters))
  print(str(rname))
  print(data.frame(cbind(INITs,NLLvals,covNANs))[1:(fi-1),])
  print(opt)
  print(obj$gr(opt$par))
  print(paste("AIC:",Fitted_AIC))
  print(paste("BIC:",Fitted_BIC))
  print(ad_obj)
  sink() 
  #===========================================================================================
  
  if(!dir.exists("Model_comp"))dir.create("Model_comp")
  save(ad_obj, model_inputs, error_r, Fitted_AIC, 
       Fitted_BIC, map_tmb, ModelName, no_of_runs, obj, opt, parameters, 
       prec_mat, rep.nowt, rname, strat_area, tmb.data, TMBmodel, 
       LEN_CAT_reff, JD_CAT_reff,
       file=paste0("Model_comp/",Mnum, TMBmodel, ModelName,"rev", ".RData"))
}


INITs0<-expand.grid(
  phi_0s = c(0.01,0.4,0.5,0.9),
  a_w_0s = c(-13,-10,-12),
  tao_0s = c(0.01,1),
  sig_0s = c(0.1,0.5)
)


ni=0


#uncomment out whichever models you want to try

#----------------------------------- c("g","t","l","gt","gl","tl","b","JL")
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","+","+", "+", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+"," "," ", " ", " ", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c(" ","+"," ", " ", " ", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c(" ","+"," ", " ", " ", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c(" "," ","+", " ", " ", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c(" "," ","+", " ", " ", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+"," ", " ", " ", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+"," ", " ", " ", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+"," ","+", " ", " ", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+"," ","+", " ", " ", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c(" ","+","+", " ", " ", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c(" ","+","+", " ", " ", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+"," ", "+", " ", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+"," ", "+", " ", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+"," ","+", " ", "+", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+"," ","+", " ", "+", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c(" ","+","+", " ", " ", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c(" ","+","+", " ", " ", "+","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", " ", " ", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", " ", " ", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", "+", " ", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", "+", " ", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", " ", "+", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", " ", "+", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", " ", " ", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", " ", " ", "+","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", "+", "+", " ","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", "+", "+", " ","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", "+", " ", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", "+", " ", "+","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", " ", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", " ", "+", "+","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", "+", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","+", "+", "+", "+","+", " "),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});

ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","","+", "+", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","+","", "+", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","","+", "+", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","+","", "+", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("+","","", "+", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  tryCatch({fit_fun(use_X=c("","","", "+", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){}); #best model listed in the manuscript

ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","","", "", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","","", "+", "", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","","", "+", "+", "","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","","", "+", "+", "+","+", ""),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","","", "", "", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","","", "", "+", "","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","","", "", "+", "+","+", ""),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});
ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","","", "", "", "","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});

ni<-ni+1 ;  #tryCatch({fit_fun(use_X=c("","+","+", "+", "+", "+","+", "+"),Mnum=NUMfn(ni),BCtrue=FALSE,INITs=INITs0) },error=function(e){});


