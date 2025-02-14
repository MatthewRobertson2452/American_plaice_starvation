#include <TMB.hpp> 
#include <iostream>

/** \brief  Approximate inverse normal cumulative distribution function, similar to R's qnorm (one-argument case only).
 * \details
 To be replaced by more accurate version based on Rmath library.
 */


template<class Type>
Type objective_function<Type>::operator() ()
{  
  DATA_INTEGER(use_Ag);
  DATA_INTEGER(use_At);
  DATA_INTEGER(use_Al);
  DATA_INTEGER(use_Agt);
  DATA_INTEGER(use_Agl);
  DATA_INTEGER(use_Atl);
  DATA_INTEGER(use_Bc);
  DATA_INTEGER(use_JL);
  DATA_SCALAR(K_crit);
  DATA_INTEGER(nyears);
  DATA_INTEGER(nstrat); 
  DATA_INTEGER(nlen); 
  DATA_INTEGER(nlen_cat); 
  DATA_VECTOR(log_wei);        
  DATA_VECTOR(log_len);  
  DATA_INTEGER(min_len);  
  DATA_INTEGER(min_JD);
  DATA_IVECTOR(iyears);
  DATA_IVECTOR(istrat); 
  DATA_IVECTOR(ilen);  
  DATA_IVECTOR(ilen_cat);  
  DATA_IVECTOR(iJD_cat);  
  DATA_MATRIX(spatiocorr);
  DATA_SPARSE_MATRIX(diag_Mat);
  DATA_MATRIX(SIGrhoPOW);
  DATA_INTEGER(nofJD);
  DATA_IVECTOR(row_non0); 
  DATA_IVECTOR(col_non0); 
  DATA_VECTOR(area); 
  DATA_IVECTOR(len_v_sub); 
  DATA_IVECTOR(JD_v_sub); 
  DATA_IVECTOR(len_cat_v_sub); 
  DATA_IVECTOR(JD_cat_v_sub); 
  
  PARAMETER(a_w);  
  PARAMETER(log_c_w);
  PARAMETER(d_w);
  PARAMETER(e_w);
  PARAMETER(log_beta);
  PARAMETER(log_q_Ag);  
  PARAMETER(log_q_Agt);  
  PARAMETER(log_q_Agl);  
  PARAMETER(log_omeg_Ag);
  PARAMETER(log_omeg_Agt);
  PARAMETER(log_omeg_Agl);
  PARAMETER(logit_phi_At);
  PARAMETER(logit_phi_Al);
  PARAMETER(logit_phi_Agt);
  PARAMETER(logit_phi_Agl);
  PARAMETER(logit_phi_Atl_t);
  PARAMETER(logit_phi_Atl_l);
  PARAMETER(logit_phi_JL_J);
  PARAMETER(logit_phi_JL_L);
  PARAMETER(log_sigma_At);
  PARAMETER(log_sigma_Al);
  PARAMETER(log_sigma_Atl);
  PARAMETER(log_sigma_JL);
  PARAMETER(tao);
  PARAMETER_VECTOR(Delta_Ag);
  PARAMETER_VECTOR(Delta_At);
  PARAMETER_VECTOR(Delta_Al);
  PARAMETER_ARRAY(Delta_Agt);
  PARAMETER_ARRAY(Delta_Agl);
  PARAMETER_ARRAY(Delta_Atl);
  PARAMETER_ARRAY(Delta_JL);
  
  using namespace density;   
  using namespace Eigen; 
  using namespace tmbutils;
  
  Type beta = exp(log_beta);
  Type c_w = exp(log_c_w);
  
  int nobs = log_wei.size(); 
  int nnon0 = row_non0.size();
  //============================================================================
  matrix<Type> spatiocorr_now = spatiocorr;
  matrix<Type> grid_nei_now = spatiocorr;
  grid_nei_now.setZero();
  
  for(int i=0; i<nnon0; i++){
    spatiocorr_now(row_non0(i),col_non0(i))  = pow(spatiocorr(row_non0(i), col_non0(i)), tao);
    spatiocorr_now(col_non0(i),row_non0(i))  = spatiocorr_now(row_non0(i), col_non0(i));
    grid_nei_now(  row_non0(i),row_non0(i)) += spatiocorr_now(row_non0(i), col_non0(i));
    grid_nei_now(  col_non0(i),col_non0(i)) += spatiocorr_now(col_non0(i), row_non0(i));
  }
  
  SparseMatrix<Type> spatiocorr_sparse = asSparseMatrix(spatiocorr_now);
  SparseMatrix<Type> grid_nei = asSparseMatrix(grid_nei_now);
  //============================================================================
  Type zero = 0.0;
  
  array<Type> A_gtl(nstrat,nyears,nlen_cat);
  
  for(int i=0; i<nstrat; i++){ 
    for(int j=0; j<nyears; j++){
      for(int k=0; k<nlen_cat; k++){
        
        A_gtl(i,j,k)  =  a_w ;
        
        if(use_Ag==1){
          A_gtl(i,j,k) = A_gtl(i,j,k) + Delta_Ag(i) ;
        }
        if(use_At==1){
          A_gtl(i,j,k) = A_gtl(i,j,k) + Delta_At(j) ;
        }
        if(use_Al==1){
          A_gtl(i,j,k) = A_gtl(i,j,k) + Delta_Al(k) ;
        }
        if(use_Agt==1){
          A_gtl(i,j,k) = A_gtl(i,j,k) + Delta_Agt(i,j) ;
        }
        if(use_Agl==1){
          A_gtl(i,j,k) = A_gtl(i,j,k) + Delta_Agl(i,k) ;
        }
        if(use_Atl==1){
          A_gtl(i,j,k) = A_gtl(i,j,k) + Delta_Atl(j,k) ;
        }
      }
    }
  }
  //
  Type asym_sig = c_w + exp( d_w - e_w * 1000.0);
  
  array<Type> M_gtlj(nstrat,nyears,nlen_cat,nofJD);
  
  for(int i=0; i<nstrat; i++){ 
    for(int j=0; j<nyears; j++){
      for(int k=0; k<nlen_cat; k++){
        for(int m=0; m<nofJD; m++){
          Type A_gtlj = A_gtl(i,j,k) + Delta_JL(m,k) ;
          Type P_gtlj = pnorm((K_crit - A_gtlj + a_w) / asym_sig);
          M_gtlj(i,j,k,m) = -log(1 - P_gtlj);
        }
      }
    }
  }
  
  array<Type> M_gtl(nstrat,nyears,nlen_cat);
  
  for(int i=0; i<nstrat; i++){ 
    for(int j=0; j<nyears; j++){
      for(int k=0; k<nlen_cat; k++){
        M_gtl(i,j,k) = zero;
        for(int m=0; m<nofJD; m++){
          M_gtl(i,j,k) = M_gtl(i,j,k) + M_gtlj(i,j,k,m);
        }
      }
    }
  }
  
  Type area_sum = zero;
  for(int i=0; i<nstrat; i++){ 
    area_sum = area_sum + area(i);
  }
  
  array<Type> M_tl(nyears,nlen_cat);
  array<Type> log_M_tl(nyears,nlen_cat);
  
  for(int j=0; j<nyears; j++){
    for(int k=0; k<nlen_cat; k++){
      M_tl(j,k) = zero;
      for(int i=0; i<nstrat; i++){ 
        M_tl(j,k) = M_tl(j,k) + M_gtl(i,j,k) * area(i);
      }
      M_tl(j,k) = M_tl(j,k)/area_sum;
      log_M_tl(j,k) = log(M_tl(j,k));
    }
  }
  //
  
  // ===========================================================
  array<Type> A_g_bar(nstrat);
  array<Type> A_t_bar(nyears);
  array<Type> A_l_bar(nlen_cat);
  array<Type> exp_A_g_bar(nstrat);
  array<Type> exp_A_t_bar(nyears);
  array<Type> exp_A_l_bar(nlen_cat);
  
  for(int i=0; i<nstrat; i++){ 
    A_g_bar(i) = zero;
    exp_A_g_bar(i) = zero;
    for(int j=0; j<nyears; j++){
      for(int k=0; k<nlen_cat; k++){ 
        A_g_bar(i)     += A_gtl(i,j,k);
        exp_A_g_bar(i) += exp(A_gtl(i,j,k));
      }
    }
    A_g_bar(i) = A_g_bar(i) / (nyears*nlen_cat);
    exp_A_g_bar(i) = exp_A_g_bar(i) / (nyears*nlen_cat);
  }
  
  for(int j=0; j<nyears; j++){
    A_t_bar(j) = zero;
    exp_A_t_bar(j) = zero;
    for(int i=0; i<nstrat; i++){ 
      for(int k=0; k<nlen_cat; k++){ 
        A_t_bar(j)     += A_gtl(i,j,k);
        exp_A_t_bar(j) += exp(A_gtl(i,j,k));
      }
    }
    A_t_bar(j) = A_t_bar(j) / (nstrat*nlen_cat);
    exp_A_t_bar(j) = exp_A_t_bar(j) / (nstrat*nlen_cat);
  }
  
  for(int k=0; k<nlen_cat; k++){ 
    A_l_bar(k) = zero;
    exp_A_l_bar(k) = zero;
    for(int i=0; i<nstrat; i++){ 
      for(int j=0; j<nyears; j++){
        A_l_bar(k)     += A_gtl(i,j,k);
        exp_A_l_bar(k) += exp(A_gtl(i,j,k));
      }
    }
    A_l_bar(k) = A_l_bar(k) / (nstrat*nyears);
    exp_A_l_bar(k) = exp_A_l_bar(k) / (nstrat*nyears);
  }
  // ===========================================================
  
  
  
  
  // ===========================================================
  array<Type> Delta_Agt_gBar(nstrat);
  array<Type> Delta_Agt_tBar(nyears);
  array<Type> Delta_Agl_gBar(nstrat);
  array<Type> Delta_Agl_lBar(nlen_cat);
  array<Type> Delta_Atl_tBar(nyears);
  array<Type> Delta_Atl_lBar(nlen_cat);
  
  for(int i=0; i<nstrat; i++){ 
    Delta_Agt_gBar(i) = zero;
    for(int j=0; j<nyears; j++){
      Delta_Agt_gBar(i) += Delta_Agt(i,j);
    }
    Delta_Agt_gBar(i) = Delta_Agt_gBar(i) / (nyears);
  }
  for(int j=0; j<nyears; j++){
    Delta_Agt_tBar(j) = zero;
    for(int i=0; i<nstrat; i++){ 
      Delta_Agt_tBar(j) += Delta_Agt(i,j);
    }
    Delta_Agt_tBar(j) = Delta_Agt_tBar(j) / (nstrat);
  }
  
  for(int i=0; i<nstrat; i++){ 
    Delta_Agl_gBar(i) = zero;
    for(int k=0; k<nlen_cat; k++){
      Delta_Agl_gBar(i) += Delta_Agl(i,k);
    }
    Delta_Agl_gBar(i) = Delta_Agl_gBar(i) / (nlen_cat);
  }
  for(int k=0; k<nlen_cat; k++){
    Delta_Agl_lBar(k) = zero;
    for(int i=0; i<nstrat; i++){ 
      Delta_Agl_lBar(k) += Delta_Agl(i,k);
    }
    Delta_Agl_lBar(k) = Delta_Agl_lBar(k) / (nstrat);
  }
  
  for(int j=0; j<nyears; j++){ 
    Delta_Atl_tBar(j) = zero;
    for(int k=0; k<nlen_cat; k++){
      Delta_Atl_tBar(j) += Delta_Atl(j,k);
    }
    Delta_Atl_tBar(j) = Delta_Atl_tBar(j) / (nlen_cat);
  }
  for(int k=0; k<nlen_cat; k++){
    Delta_Atl_lBar(k) = zero;
    for(int j=0; j<nyears; j++){ 
      Delta_Atl_lBar(k) += Delta_Atl(j,k);
    }
    Delta_Atl_lBar(k) = Delta_Atl_lBar(k) / (nyears);
  }
  // ===========================================================
  array<Type> Delta_JL_jBar(nofJD);
  array<Type> Delta_JL_lBar(nlen_cat);
  
  if(use_JL==1){
    for(int l0=0; l0<nlen_cat; l0++){
      Delta_JL_lBar(l0) = zero;
      for(int j0=0; j0<nofJD; j0++){
        Delta_JL_lBar(l0) += Delta_JL(j0,l0);
      }
      Delta_JL_lBar(l0) = Delta_JL_lBar(l0) / (nofJD);
    }
    
    for(int j0=0; j0<nofJD; j0++){
      Delta_JL_jBar(j0) = zero;
      for(int l0=0; l0<nlen_cat; l0++){
        Delta_JL_jBar(j0) += Delta_JL(j0,l0);
      }
      Delta_JL_jBar(j0) = Delta_JL_jBar(j0) / (nlen_cat);
    }
  }
  // ===========================================================
  
  Type nll = zero;
  
  if(use_Ag==1){
    Type q_Ag = exp(log_q_Ag); 
    Type omeg_Ag = exp(log_omeg_Ag); 
    SparseMatrix<Type> prec_Ag = q_Ag*(grid_nei - spatiocorr_sparse + omeg_Ag * diag_Mat);
    nll += GMRF(prec_Ag)(Delta_Ag);
  }
  
  if(use_At==1){
    Type sigma_At = exp(log_sigma_At); 
    Type phi_At = Type(1)/(Type(1)+exp(-logit_phi_At)); // (-Inf,Inf) -> (0,1)
    nll += SCALE(AR1(phi_At),sigma_At)(Delta_At) ;
  }
  
  if(use_Al==1){
    Type sigma_Al = exp(log_sigma_Al);
    Type phi_Al = Type(1)/(Type(1)+exp(-logit_phi_Al)); // (-Inf,Inf) -> (0,1)
    nll += SCALE(AR1(phi_Al),sigma_Al)(Delta_Al) ;
  }
  
  if(use_Agt==1){
    Type q_Agt = exp(log_q_Agt); 
    Type omeg_Agt = exp(log_omeg_Agt); 
    Type phi_Agt = Type(1)/(Type(1)+exp(-logit_phi_Agt)); // (-Inf,Inf) -> (0,1)
    SparseMatrix<Type> prec_Agt = q_Agt*(grid_nei - spatiocorr_sparse + omeg_Agt * diag_Mat);
    nll += SEPARABLE(AR1(phi_Agt), GMRF(prec_Agt))(Delta_Agt);
  }
  
  if(use_Agl==1){
    Type q_Agl = exp(log_q_Agl); 
    Type omeg_Agl = exp(log_omeg_Agl); 
    Type phi_Agl = Type(1)/(Type(1)+exp(-logit_phi_Agl)); // (-Inf,Inf) -> (0,1)
    SparseMatrix<Type> prec_Agl = q_Agl*(grid_nei - spatiocorr_sparse + omeg_Agl * diag_Mat);
    nll += SEPARABLE(AR1(phi_Agl), GMRF(prec_Agl))(Delta_Agl);
  }
  
  if(use_Atl==1){
    Type sigma_Atl = exp(log_sigma_Atl);
    Type phi_Atl_t = Type(1)/(Type(1)+exp(-logit_phi_Atl_t)); // (-Inf,Inf) -> (0,1)
    Type phi_Atl_l = Type(1)/(Type(1)+exp(-logit_phi_Atl_l)); // (-Inf,Inf) -> (0,1)
    nll += SCALE(SEPARABLE(AR1(phi_Atl_l),AR1(phi_Atl_t)),sigma_Atl)(Delta_Atl) ;
  }
  

    Type sigma_JL = exp(log_sigma_JL);
    Type phi_JL_J = Type(1)/(Type(1)+exp(-logit_phi_JL_J)); // (-Inf,Inf) -> (0,1)
    Type phi_JL_L = Type(1)/(Type(1)+exp(-logit_phi_JL_L)); // (-Inf,Inf) -> (0,1)
    //============================================================================
    matrix<Type> SIG = SIGrhoPOW * 0.0;
    if(use_JL==1){
    for(int ijd=0; ijd<nofJD; ijd++){
      for(int jjd=0; jjd<nofJD; jjd++){
        SIG(ijd,jjd) = pow(phi_JL_J, SIGrhoPOW(ijd,jjd)) / (1 - pow(phi_JL_J,2));
      }
    }
    //============================================================================
    nll += SCALE(SEPARABLE(AR1(phi_JL_L),MVNORM(SIG)),sigma_JL)(Delta_JL) ;
  }
  
  
  
  // ==================================================================
  // Area weighted mean weight at length
  // ==================================================================
  int nlength_sub = len_v_sub.size();
  int niJD_cat = JD_v_sub.size();
  array<Type> EcountA_sub(nyears,nlength_sub,niJD_cat);
  EcountA_sub.setZero();
  array<Type> Ecount_weightA_sub(nyears,nlength_sub,niJD_cat);
  Ecount_weightA_sub.setZero();
  array<Type> Eweight_yearA_sub(nyears,nlength_sub,niJD_cat);
  Eweight_yearA_sub.setZero();

  for(int l=0;l<niJD_cat;l++){
    for(int j=0;j<nyears;j++){
      for(int k=0;k<nlength_sub;k++){
        for(int i=0;i<nstrat;i++){
          int l_now = len_v_sub(k);

          Type weight_now = exp( A_gtl(i,j,len_cat_v_sub(k)) );
          if(use_Bc==1){
            weight_now = weight_now * exp(beta*log(l_now)) ;
          }

          if(use_JL==1){
            weight_now = weight_now * exp(Delta_JL(JD_cat_v_sub(l),len_cat_v_sub(k))) ;
          }

          EcountA_sub(j,k,l)        += area(i);
          Ecount_weightA_sub(j,k,l) += area(i)*weight_now;
        }
        Eweight_yearA_sub(j,k,l) = Ecount_weightA_sub(j,k,l)/EcountA_sub(j,k,l);
      }
    }
  }
  // ==================================================================
  
  vector<Type> res(nobs);
  vector<Type> sig_w(nobs);
  vector<Type> pred_log_w(nobs);
  
  for(int i = 0;i < nobs;++i){  
    
    Type A_gtl_i = A_gtl(istrat(i), iyears(i), ilen_cat(i));
    
    pred_log_w(i) = A_gtl_i ; 
    
    if(use_Bc==1){
      pred_log_w(i) += beta*log_len(i) ; 
    }
    
    if(use_JL==1){
      pred_log_w(i) += Delta_JL(iJD_cat(i),ilen_cat(i));  
    }
    
    sig_w(i) = c_w + exp( d_w - e_w * log_len(i));
    
    res(i) = log_wei(i) - pred_log_w(i);
    nll   -= dnorm(res(i), zero, sig_w(i), true);
  }
  
  vector<Type> res_st = res/sig_w;
  
  REPORT(Eweight_yearA_sub);
  ADREPORT(Eweight_yearA_sub);
  
  Type q_Agt_inv = 1.0/exp(log_q_Agt);
  Type q_Agl_inv = 1.0/exp(log_q_Agl);
  Type q_Ag_inv  = 1.0/exp(log_q_Ag);
  Type root_q_Agt_inv = pow(q_Agt_inv,0.5);
  Type root_q_Agl_inv = pow(q_Agl_inv,0.5);
  Type root_q_Ag_inv  = pow(q_Ag_inv,0.5);
  
  REPORT(asym_sig);
  REPORT(a_w);
  REPORT(A_gtl);
  REPORT(A_g_bar);
  REPORT(A_t_bar);
  REPORT(A_l_bar);
  REPORT(exp_A_g_bar);
  REPORT(exp_A_t_bar);
  REPORT(exp_A_l_bar);
  ADREPORT(A_g_bar);
  ADREPORT(A_t_bar);
  ADREPORT(A_l_bar);
  REPORT(M_gtlj);
  REPORT(log_M_tl);
  ADREPORT(log_M_tl);
  
  if(use_Ag==1){
    REPORT(Delta_Ag);
    ADREPORT(Delta_Ag);
    REPORT(q_Ag_inv);
    ADREPORT(q_Ag_inv);
    REPORT(root_q_Ag_inv);
    ADREPORT(root_q_Ag_inv);
  }
  if(use_At==1){
    REPORT(Delta_At);
    ADREPORT(Delta_At);
  }
  if(use_Al==1){
    REPORT(Delta_Al);
    ADREPORT(Delta_Al);
  }
  if(use_Agt==1){
    REPORT(Delta_Agt);
    REPORT(Delta_Agt_gBar);
    ADREPORT(Delta_Agt_gBar);
    REPORT(Delta_Agt_tBar);
    ADREPORT(Delta_Agt_tBar);
    REPORT(q_Agt_inv);
    ADREPORT(q_Agt_inv);
    REPORT(root_q_Agt_inv);
    ADREPORT(root_q_Agt_inv);
  }
  if(use_Agl==1){
    REPORT(Delta_Agl);
    REPORT(Delta_Agl_gBar);
    ADREPORT(Delta_Agl_gBar);
    REPORT(Delta_Agl_lBar);
    ADREPORT(Delta_Agl_lBar);
    REPORT(q_Agl_inv);
    ADREPORT(q_Agl_inv);
    REPORT(root_q_Agl_inv);
    ADREPORT(root_q_Agl_inv);
  }
  if(use_Atl==1){
    REPORT(Delta_Atl);
    REPORT(Delta_Atl_tBar);
    ADREPORT(Delta_Atl_tBar);
    REPORT(Delta_Atl_lBar);
    ADREPORT(Delta_Atl_lBar);
  }
  if(use_JL==1){
    REPORT(Delta_JL);
    REPORT(SIG);
    REPORT(Delta_JL_jBar);
    REPORT(Delta_JL_lBar);
    REPORT(phi_JL_L);
    REPORT(phi_JL_J);
    REPORT(M_gtlj);
    REPORT(M_gtl);
    ADREPORT(Delta_JL);
    ADREPORT(Delta_JL_jBar);
    ADREPORT(Delta_JL_lBar);
  }
  
  REPORT(res_st);
  REPORT(res);
  REPORT(pred_log_w);
  REPORT(sig_w);
  ADREPORT(sig_w);
  
  return nll;
}

