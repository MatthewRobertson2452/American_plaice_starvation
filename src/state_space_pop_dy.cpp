
#include <TMB.hpp> 
#include "pnorm4.h" //Atomic functions for censored likelihoods
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  
  // input data;  
  DATA_MATRIX(weight); 
  DATA_MATRIX(mat); 
  DATA_MATRIX(midy_weight);
  DATA_MATRIX(C);         
  DATA_VECTOR(index);
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(iage);
  DATA_IVECTOR(isurvey); 
  DATA_VECTOR(log_landings);
  DATA_VECTOR(fs);
  DATA_INTEGER(A);
  DATA_INTEGER(Y);
  DATA_IVECTOR(is_year);
  DATA_INTEGER(Ns);
  DATA_IVECTOR(isurvey1);
  DATA_IVECTOR(isd);
  DATA_VECTOR(ssb_start);
  DATA_SCALAR(std_landings);
  DATA_SCALAR(rec_covar_switch);
  DATA_SCALAR(rec_nll_switch);
  DATA_SCALAR(M_RE_switch);
  DATA_SCALAR(M_nll_switch);
  DATA_SCALAR(M_covar_switch);
  DATA_SCALAR(YE_switch);
  DATA_VECTOR(mean_midy_weight);
  DATA_ARRAY(crl);
  DATA_VECTOR(log_condM_dat);
  DATA_SCALAR(log_sdcond);
  
  int n = index.size();  
  Type one = 1.0;
  Type zero = 0.0;
  
  matrix<Type> log_C = C.array().log();
  vector<Type> log_index = log(index);
  
  //define parameters;  
  PARAMETER_VECTOR(log_No);  
  PARAMETER_VECTOR(log_Rec_mean);
  PARAMETER_VECTOR(log_std_log_R);
  PARAMETER_ARRAY(m_q);
  PARAMETER_VECTOR(log_cv_index); 
  PARAMETER(logit_ar_logF_year);        
  PARAMETER(logit_ar_logF_age); 
  PARAMETER_VECTOR(log_M_dev);
  PARAMETER(log_sd_M);
  PARAMETER_VECTOR(log_Rec_dev);   
  PARAMETER_ARRAY(log_F); 
  PARAMETER(log_alpha);
  PARAMETER(log_beta);
  PARAMETER(log_init_rec);
  PARAMETER_VECTOR(log_intercept_R);
  PARAMETER_VECTOR(log_phi_R);
  PARAMETER_VECTOR(intercept_M);
  PARAMETER_VECTOR(phi_M);
  PARAMETER(logit_log_R);
  PARAMETER(logit_ar_logM_year);        
  PARAMETER(logit_ar_logM_age); 
  PARAMETER(logit_ar_iye_year_spr); 
  PARAMETER(logit_ar_iye_year_fall); 
  PARAMETER(log_std_iye_spr);
  PARAMETER(log_std_iye_fall);
  PARAMETER_VECTOR(iye_spr);
  PARAMETER_VECTOR(iye_fall);
  
  PARAMETER_ARRAY(log_std_crl);
  PARAMETER_VECTOR(logit_ar_crl);
  PARAMETER(log_sdrwM); 
  PARAMETER_VECTOR(log_condM);
  PARAMETER(init_log_condM);
  
  Type sdrwM = exp(log_sdrwM);
  Type sdcond = exp(log_sdcond);
  
  Type phi_log_R = exp(logit_log_R)/(one + exp(logit_log_R));
  
  Type ar_iye_year_spr = exp(logit_ar_iye_year_spr)/(one + exp(logit_ar_iye_year_spr));
  Type ar_iye_year_fall = exp(logit_ar_iye_year_fall)/(one + exp(logit_ar_iye_year_fall));
  
  Type std_iye_spr = exp(log_std_iye_spr);
  Type std_iye_fall = exp(log_std_iye_fall);
  
  vector<Type> intercept_R = log_intercept_R;
  vector<Type> phi_R = log_phi_R;
  
  Type init_rec = exp(log_init_rec);
  vector<Type> Rec_dev = exp(log_Rec_dev);
  Type sd_M = exp(log_sd_M);
  
  vector<Type> std_log_R = exp(log_std_log_R); 
  
  vector<Type> cv_index = exp(log_cv_index);
  
  Type ar_logF_age = exp(logit_ar_logF_age)/(one + exp(logit_ar_logF_age));    
  Type ar_logF_year = exp(logit_ar_logF_year)/(one + exp(logit_ar_logF_year)); 
  Type ar_logM_age = exp(logit_ar_logM_age)/(one + exp(logit_ar_logM_age));    
  Type ar_logM_year = exp(logit_ar_logM_year)/(one + exp(logit_ar_logM_year)); 
  
  matrix<Type> log_N(Y,A); 
  matrix<Type> N(Y,A);     
  matrix<Type> F(Y,A-4);    
  matrix<Type> Z(Y,A);  
  matrix<Type> EC(Y,A-4);  
  matrix<Type> log_ECW(Y,A-4); 
  matrix<Type> ECW(Y,A-4); 
  matrix<Type> C_resid(Y,A); 
  matrix<Type> C_resid_std(Y,A);     
  matrix<Type> std_C_resid(Y,A);
  matrix<Type> B_matrix(Y,A);           
  matrix<Type> SSB_matrix(Y,A); 
  
  vector<Type> Elog_index(n);
  vector<Type> resid_index(n); 
  vector<Type> std_resid_index(n); 
  
  Type ar_crl_age = invlogit(logit_ar_crl(0));
  Type ar_crl_year = invlogit(logit_ar_crl(1));
  
  //**********  SD report objects ***************;
  
  vector<Type> biomass(Y); 
  vector<Type> log_biomass(Y);  
  vector<Type> ssb(Y);
  ssb(0)=ssb_start(0); 
  vector<Type> log_ssb = log(ssb);
  vector<Type> aveF_69(Y);
  vector<Type> log_aveF_69(Y);   
  vector<Type> aveF_914(Y);
  vector<Type> log_aveF_914(Y); 
  
  vector<Type> log_condMest(Y-30);
  vector<Type> condM;
  
  Type alpha = exp(log_alpha);
  Type beta = exp(log_beta);
  
  array<Type> M(Y,A);
  
  //**********  start the engine ***************;
  
  using namespace density;
  Type nll = 0.0;  
  vector<Type> jnll(26);
  jnll.setZero();
  
 //prior for sd_M
  if(M_RE_switch==1){
    jnll(22)-=dnorm(sd_M, Type(0.0), Type(0.05), true);
    nll += jnll(22);
  }
  
  
  log_condMest(0) = init_log_condM;
  jnll(23) -= dnorm(log_condM(1),init_log_condM,sdrwM, true);
  log_condMest(1) = log_condM(1);
  for(int y = 2;y < Y-30;++y){ 
    jnll(24) -= dnorm(log_condM(y),log_condM(y-1),sdrwM, true);
    log_condMest(y) = log_condM(y);
  }
  
  nll += jnll(23);
  nll += jnll(24);

  condM = exp(log_condMest);
  
  vector<Type> enviro_impact(Y);
  vector<Type> log_enviro_impact(Y);
  
  // The cohort model;
  vector<Type> Rec(Y);
  vector<Type> log_Rec(Y);
  
  log_Rec(0) = init_rec;
  log_N(0,0) = init_rec;
  for(int a = 1;a < A;++a){
    log_N(0,a) = log_No(a-1);
  }
  
  for(int a = 0;a < A;++a){
    N(0,a) = exp(log_N(0,a));
    
    B_matrix(0,a) = weight(0,a)*N(0,a);
    SSB_matrix(0,a) = mat(0,a)*B_matrix(0,a);
    
    ssb(0) += SSB_matrix(0,a);
  }
  
  
  
  for(int a = 0;a < A;++a){ 
    for(int y = 0;y < Y;++y){
      
      if(M_RE_switch==0){     
        if(a<3){
          M(y,a) = 0.5;
          Z(y,a) = M(y,a);
        }
        if(a==3){
          M(y,a) = 0.3;
          Z(y,a) = M(y,a);
        }
      }
      
      if(M_RE_switch==1){
        
        if(M_covar_switch==0){
          if(a==0){
            M(y,a) = 0.2 * exp(log_M_dev(y));
            Z(y,a) = M(y,a);
          }
          if(a<3 & a>0){
            M(y,a) = 0.2 * exp(log_M_dev(y));
            Z(y,a) = M(y,a);
          }
          if(a==3){
            M(y,a) = 0.2 * exp(log_M_dev(y));
            Z(y,a) = M(y,a);
          }
        }
        if(M_covar_switch==1){
          if(a==0){
            if(y<=29){
              M(y,a) = 0.2 * exp(log_M_dev(y));
              Z(y,a) = M(y,a);
            }
            if(y>29){
              M(y,a) = 0.2 * exp(log_M_dev(y)) + exp(phi_M(0)*log_condMest(y-30));
              Z(y,a) = M(y,a);
            }
          }
          if(a<3 & a>0){
            if(y<=29){
              M(y,a) = 0.2 * exp(log_M_dev(y));
              Z(y,a) = M(y,a);
            }
            if(y>29){
              M(y,a) = 0.2 * exp(log_M_dev(y)) + exp(phi_M(0)*log_condMest(y-30));
              Z(y,a) = M(y,a);
            }
          }
          if(a==3){
            if(y<=29){
              M(y,a) = 0.2 * exp(log_M_dev(y));
              Z(y,a) = M(y,a);
            }
            if(y>29){
              M(y,a) = 0.2 * exp(log_M_dev(y)) + exp(phi_M(0)*log_condMest(y-30));
              Z(y,a) = M(y,a);
            }
          }
        }
      }
      
      if(a>=4){
        F(y,a-4) = exp(log_F(y,a-4));
        
        if(M_RE_switch==0){
          
          if(M_covar_switch==0){
            M(y,a) = 0.2;
            Z(y,a) = F(y,a-4) + M(y,a);
          }
          
          if(M_covar_switch==2){
            if(y<=29){
              M(y,a) = 0.2;
              Z(y,a) =  F(y,a-4) + M(y,a);}
            if(y>29){
              M(y,a) = 0.2 * exp(intercept_M(0) +phi_M(0)*condM(y-30));
              Z(y,a) =  F(y,a-4) + M(y,a);}
          }
        }
        
        if(M_RE_switch==1){
          if(M_covar_switch==0){
            M(y,a) = 0.2 * exp(log_M_dev(y));
            Z(y,a) =  F(y,a-4) + M(y,a);
          }
          
          //cond_M
          if(M_covar_switch==1){
            if(y<=29){
              M(y,a) = 0.2 * exp(log_M_dev(y));
              Z(y,a) =  F(y,a-4) + M(y,a);}
            if(y>29){
              M(y,a) = 0.2 * exp(log_M_dev(y)) + exp(phi_M(0)*log_condMest(y-30));
              Z(y,a) =  F(y,a-4) + M(y,a);}
          }
          
        }
        
        
      }
    }
  }
  
  
  for(int y = 1;y < Y;++y){
    
    //no enviro
    if(rec_covar_switch==0){
      if(y<30){
        log_Rec(y) = log_Rec_mean(0) + log_Rec_dev(y-1);}
      if(y>29){
        log_Rec(y) = log_Rec_mean(1) + log_Rec_dev(y-1);}
    }
    
    log_N(y,0) = log_Rec(y);  
    for(int a = 1;a < A;++a){
      log_N(y,a) = (log_N(y-1,a-1) - Z(y-1,a-1));
      if(a == A-1){//Plus group
        log_N(y,a) = log(exp(log_N(y-1,a-1)-Z(y-1,a-1))+exp(log_N(y-1,a)-Z(y-1,a)));
      }
    }
    
    for(int a = 0;a < A;++a){
      N(y,a) = exp(log_N(y,a));
      
      B_matrix(y,a) = weight(y,a)*N(y,a);
      SSB_matrix(y,a) = mat(y,a)*B_matrix(y,a);
      
      ssb(y) += SSB_matrix(y,a);
    }
    
  }
  
  Rec = exp(log_Rec);
  
  
  vector<Type> landings_pred(Y);
  vector<Type> C_tot(Y);
  vector<Type> CW_tot(Y);
  
  //Baranov catch equation predictions and residuals
  for(int y = 0;y < Y;++y){
    C_tot(y) = zero;
    CW_tot(y) = zero;
    for(int a = 0;a < A-4;++a){
      EC(y,a) = N(y,a+4)*((one - exp(-one*Z(y,a+4)))*F(y,a)/Z(y,a+4));
      ECW(y,a) = EC(y,a)*midy_weight(y,a);
      C_tot(y) += EC(y,a);
      CW_tot(y) += ECW(y,a);
    } }
  
  vector<Type> log_landings_pred = log(CW_tot);
  landings_pred = CW_tot;
  vector<Type> landings_resids = log_landings - log_landings_pred;
  vector<Type> std_landings_resid = landings_resids/std_landings;
  
  matrix<Type> p_ya(A-4,Y);
  matrix<Type> pya_sum(A-4,Y);
  matrix<Type> pi_ya(A-4,Y);
  matrix<Type> pred_crl(Y,A-5);
  matrix<Type> resid_crl(Y,A-5);
  matrix<Type> std_resid_crl(Y,A-5);
  array<Type> std_crl(Y,A-5);
  
  int i,j;
  for(i = 0;i < Y;++i){
    for(j = 0;j < A-5;++j){std_crl(i,j) = exp(log_std_crl(i,j));}
  }
  
  //age composition catch
  for(int i = 0;i < Y;++i){
    for(int j = 0;j < A-4;++j){
      p_ya(j,i) = EC(i,j)/C_tot(i);
    }}
  
  Type total;
  for(int i = 0;i <Y;++i){
    total=zero;
    pya_sum(0,i) = one;
    for(int j = 1;j < A-4;++j){
      total+=p_ya(j-1,i);
      pya_sum(j,i) = (one-total);
    }}
  
  for(int i = 0;i <Y;++i){
    for(int j = 0;j < A-4;++j){
      pi_ya(j,i) = p_ya(j,i)/pya_sum(j,i);}}
  
  
  for(int i = 0;i<Y;i++){
    for(int j = 0;j<A-5;j++){
      pred_crl(i,j)= log(pi_ya(j,i)/(one - pi_ya(j,i)));
      resid_crl(i,j) = crl(i,j) - pred_crl(i,j);
      std_resid_crl(i,j) = resid_crl(i,j)/std_crl(i,j);}}
  
  matrix<Type> q(Ns,A);
  
  int ia,iy,iy1,icv,is;   
  
  for(int i = 0;i < Ns;++i){
    for(int j = 1;j < A-1;++j){
      m_q(i,j) = exp(m_q(i,j));
    }}
  
  
  for(int i = 0;i < Ns;++i){
    q(i,0) = m_q(i,0);
    
    for(int j=1;j<A-1;++j){
      q(i,j) = q(i,j-1)+m_q(i,j);
    }
    int j=A-1;
    q(i,j) = q(i,j-1);}
  
  
  vector<Type> std_index_vec(n);
  matrix<Type> mresid_index(Ns,A);
  matrix<Type> msd_index(Ns,A);
  vector<Type> Eindex(n);
  
  //  Survey index predictions, and residuals;
  for(int i = 0;i < n;++i){
    ia = iage(i);
    iy = iyear(i);
    is = isurvey(i);
    iy1 = is_year(i);
    icv = isd(i);
    
    Elog_index(i) = q(iy1,ia) + log(N(iy,ia)) - fs(i)*Z(iy,ia);
    Eindex(i) = exp(Elog_index(i));
    std_index_vec(i) = cv_index(icv);
    if(YE_switch==0){resid_index(i) = index(i) - Elog_index(i);}
    if(YE_switch>0){
      if(is==0){resid_index(i) = index(i) - Elog_index(i)+iye_spr(iy);} 
      if(is==1){resid_index(i) = index(i) - Elog_index(i)+iye_fall(iy);} 
    }
    std_resid_index(i) = resid_index(i)/std_index_vec(i);
    
    mresid_index(iy1,ia) = resid_index(i);
    msd_index(iy1,ia) = std_index_vec(i);
  }
  
  // End of model, now fit functions (nll - negative loglikelihoods);

  
  ////SURVEY
  
  vector<Type> del(A);
  vector<Type> sd_del(A);
  
  for(int i = 0;i < n;++i){
    jnll(0) -= dnorm(resid_index(i),zero,std_index_vec(i),true);
  }
  
  nll += jnll(0);
  
  // //landings nll
  vector<Type> std_landings_resids(Y);
  for(int y = 0;y < Y;++y){
    if(y<18){
      std_landings_resids(y) = landings_resids(y)/(std_landings);
      nll -= dnorm(landings_resids(y),zero,std_landings,true);
    }
    if(y<34 & y>17){
      std_landings_resids(y) = landings_resids(y)/(std_landings);
      nll -= dnorm(landings_resids(y),zero,std_landings,true);
    }
    if(y>33){
      std_landings_resids(y) = landings_resids(y)/(std_landings);
      nll -= dnorm(landings_resids(y),zero,std_landings,true);
    }
  }
  
  
  array<Type> temp(Y,A-5);
  temp = resid_crl.array();
  jnll(3)  += VECSCALE(SEPARABLE(AR1(ar_crl_age),AR1(ar_crl_year)),std_crl)(temp);
  
  
  nll += jnll(3);
  
  //index year effect nll;
  if(YE_switch==1){
    jnll(5) += SCALE(AR1(ar_iye_year_spr),std_iye_spr)(iye_spr);
    jnll(6) += SCALE(AR1(ar_iye_year_fall),std_iye_fall)(iye_fall);
  }
  
  nll += jnll(5);
  nll += jnll(6);
  
  //index year effect nll;
  if(YE_switch==2){
    for(int y = 0;y < Y;++y){
      jnll(7) -= dnorm(iye_spr(y), zero, std_iye_spr, true);
      jnll(8) -= dnorm(iye_fall(y), zero, std_iye_fall, true);
    }
  }
  
  nll += jnll(7);
  nll += jnll(8);
  
  vector<Type> resid_condM(Y-30);
  for(int y = 30;y < Y;++y){
    resid_condM(y-30) = log_condM_dat(y-30)-log_condMest(y-30);
    jnll(25) -= dnorm(resid_condM(y-30), zero, sdcond, true);
  }
  nll += jnll(25);
  
  /////RECRUITMENT
  
  //Normal
  if(rec_nll_switch==0){
    for(int y = 0;y < Y;++y){
      jnll(9) -= dnorm(log_Rec_dev(y), zero, std_log_R(0), true);
    }
  }
  
  nll += jnll(9);
  
  //AR1(Y)
  if(rec_nll_switch==1){
    jnll(10) += SCALE(AR1(phi_log_R),std_log_R(0))(log_Rec_dev);
  }
  
  nll += jnll(10);
  
  //RW
  if(rec_nll_switch==2){
    jnll(11) -= dnorm(log_Rec_dev(0), zero, std_log_R(0), true);
    for(int y = 1;y < Y;++y){
      if(y<30){
        jnll(12) -= dnorm(log_Rec_dev(y)-log_Rec_dev(y-1), zero, std_log_R(0), true);}
      if(y>29){
        jnll(12) -= dnorm(log_Rec_dev(y)-log_Rec_dev(y-1), zero, std_log_R(1), true);}
    }
  }
  
  nll += jnll(11);
  nll += jnll(12);
  
  /////NATURAL MORTALITY
  
  //NORMAL
  if(M_RE_switch==1){
    if(M_nll_switch==0){
      for(int y = 0;y < Y;++y){
        jnll(13) -= dnorm(log_M_dev(y), zero, sd_M, true);
      }
    }
  }
  
  nll += jnll(13);
  
  //AR1
  if(M_nll_switch==1){
    //year x age correlation on first+1:last ages;
    jnll(16) += SCALE(AR1(ar_logM_year), sd_M)(log_M_dev);
  }
  
  nll += jnll(16);
  
  if(M_nll_switch==2){
    jnll(17) -= dnorm(log_M_dev(0),Type(0.0),sd_M, true);
    for(int y = 1;y < Y;++y){ 
      jnll(17) -= dnorm(log_M_dev(y),log_M_dev(y-1),sd_M, true);
    }
  }
  
  
  nll += jnll(17);
  
  // //RW on first age;
  vector<Type> delF = log_F.col(0);
  jnll(19) -= dnorm(delF(0),Type(-10.0),one, true);
  for(int y = 1;y < Y;++y){ 
    jnll(19) -= dnorm(delF(y),delF(y-1),one, true);
  }
  array<Type> log_F1(Y,A-5);
  for(int a = 1;a < A-4;++a){
    log_F1.col(a-1) = log_F.col(a);
  } 
   
  // //year x age correlation on first+1:last ages;
  jnll(19) += SEPARABLE(AR1(ar_logF_age),AR1(ar_logF_year))(log_F1);

  nll += jnll(19);
  
  // create some REPORT output;  
  
  for(int y = 0;y < Y;++y){
    biomass(y) = zero;     
    ssb(y) = zero;
    for(int a = 0;a < A;++a){
      biomass(y) += B_matrix(y,a); 
      ssb(y) += SSB_matrix(y,a);}
  }          
  log_biomass = log(biomass);
  log_ssb = log(ssb);       
  
  //pop size weighted ave F;  
  
  Type tni;
  
  for(int y = 0;y < Y;++y){
    aveF_69(y) = zero; 
    tni = zero;
    for(int a = 1;a <= 4;++a){
      aveF_69(y) += F(y,a)*N(y,a); 
      tni += N(y,a);
    }
    aveF_69(y) = aveF_69(y)/tni;  
    aveF_914(y) = zero; 
    tni = zero;
    for(int a = 4;a <= 9;++a){
      aveF_914(y) += F(y,a)*N(y,a); 
      tni += N(y,a);
    }
    aveF_914(y) = aveF_914(y)/tni;
  }
  
  log_aveF_69 = log(aveF_69);
  log_aveF_914 = log(aveF_914);
  
  vector<Type> N_total = N.rowwise().sum() ;
  
  REPORT(log_Rec_mean);
  REPORT(std_crl);
  REPORT(jnll);
  REPORT(iye_spr);
  REPORT(iye_fall);
  REPORT(M);
  REPORT(phi_R);
  REPORT(intercept_R);
  REPORT(phi_M);
  REPORT(intercept_M);
  REPORT(M);
  REPORT(sd_M);
  REPORT(std_log_R);  
  REPORT(ar_logF_age);           
  REPORT(ar_logF_year);
  REPORT(log_M_dev);
  REPORT(N);          
  REPORT(F);                  
  REPORT(Z);                  
  REPORT(B_matrix);             
  REPORT(SSB_matrix);                
  REPORT(biomass);                  
  REPORT(ssb);            
  REPORT(aveF_69);
  REPORT(aveF_914);
  REPORT(EC);                          
  REPORT(C_resid);             
  REPORT(std_C_resid);  
  REPORT(ECW);        
  REPORT(landings_pred); 
  REPORT(C_tot);
  REPORT(CW_tot);
  REPORT(resid_crl);
  REPORT(std_resid_crl);
  REPORT(Elog_index);                 
  REPORT(resid_index);             
  REPORT(std_resid_index); 
  REPORT(std_landings_resids);
  REPORT(Rec);
  REPORT(log_Rec);
  REPORT(log_Rec_dev);      
  REPORT(log_F);      
  REPORT(alpha);
  REPORT(beta);
  REPORT(q);
  REPORT(condM);
  REPORT(sdrwM);
  REPORT(sdcond);
  REPORT(std_index_vec);
  REPORT(cv_index);
  REPORT(Eindex);
  REPORT(mresid_index);
  REPORT(msd_index);
  REPORT(del);
  REPORT(sd_del);
  REPORT(landings_resids);
  REPORT(pred_crl);
  REPORT(resid_condM);
  REPORT(log_condM);
  REPORT(log_condMest);
  
  
  ADREPORT(intercept_M);
  ADREPORT(phi_M);
  ADREPORT(intercept_R);
  ADREPORT(phi_R);
  ADREPORT(M);
  ADREPORT(log_M_dev);
  ADREPORT(log_Rec_dev);
  ADREPORT(log_init_rec);
  ADREPORT(N_total); 
  ADREPORT(log_landings_pred);     
  ADREPORT(log_biomass);     
  ADREPORT(log_Rec);
  ADREPORT(log_ssb);   
  ADREPORT(log_aveF_914);  
  ADREPORT(ar_logF_age);           
  ADREPORT(ar_logF_year);  
  
  return nll;
}

