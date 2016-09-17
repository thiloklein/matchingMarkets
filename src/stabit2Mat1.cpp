#define ARMA_NO_DEBUG
// Uncomment the above line if you want to disable all run-time checks.
// This will result in faster code, but you first need to make sure that your code runs correctly!

// [[Rcpp::depends(RcppArmadillo,RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>

using namespace Rcpp;
arma::colvec mvrnormArma(arma::colvec mu, arma::mat sigma, int ncols);
double norm_rs(double a, double b);
double half_norm_rs(double a, double b);
double unif_rs(double a, double b);
double exp_rs(double a, double b);
double truncn2(double mu, double sigma, double lower, double upper);

//' @export
// [[Rcpp::export]]
List stabit2Mat1(Rcpp::List Cr, Rcpp::List Cmatchr, Rcpp::List Sr, Rcpp::List Smatchr, 
  Rcpp::List Dr, Rcpp::List dr, Rcpp::List Mr, Rcpp::List Hr, arma::colvec nCollegesr, 
  arma::colvec nStudentsr, Rcpp::List CCr, Rcpp::List SSr, Rcpp::List CCmatchr, Rcpp::List SSmatchr,
  Rcpp::List Lr, Rcpp::List studentIdsr, Rcpp::List collegeIdr, arma::colvec nEquilibsr,
  Rcpp::List sopt2equr, Rcpp::List equ2soptr, Rcpp::List copt2equr, Rcpp::List equ2coptr, arma::colvec coptidr,
  int n, int N, int niter, int T, int thin, bool display_progress = true) {
  
  // ---------------------------------------------
  // Independent variables
  // ---------------------------------------------
    
  arma::field<arma::mat> C(T), Cmatch(T), S(T), Smatch(T), CC(T), SS(T), CCmatch(T), SSmatch(T), collegeId(T); 
  arma::field<arma::uvec> L(T), d(T);
  arma::field<arma::umat> M(T), studentIds(N), sopt2equ(T), equ2sopt(T), copt2equ(T), equ2copt(T);
  arma::field<arma::colvec> D(T); 
  arma::field<arma::cube> H(T);

  for(int i=0; i<T; i++){
    C(i) = Rcpp::as<arma::mat>(Cr[i]);  
    Cmatch(i) = Rcpp::as<arma::mat>(Cmatchr[i]);  
    S(i) = Rcpp::as<arma::mat>(Sr[i]);  
    Smatch(i) = Rcpp::as<arma::mat>(Smatchr[i]);  
    D(i) = Rcpp::as<arma::colvec>(Dr[i]); 
    d(i) = Rcpp::as<arma::uvec>(dr[i]);  // uvec is used for logical indexing
    M(i) = Rcpp::as<arma::umat>(Mr[i]);  // umat is used for logical indexing
    H(i) = Rcpp::as<arma::cube>(Hr[i]);  
    CC(i) = Rcpp::as<arma::mat>(CCr[i]);  
    SS(i) = Rcpp::as<arma::mat>(SSr[i]); 
    CCmatch(i) = Rcpp::as<arma::mat>(CCmatchr[i]); 
    SSmatch(i) = Rcpp::as<arma::mat>(SSmatchr[i]); 
    L(i) = Rcpp::as<arma::uvec>(Lr[i]);  // uvec is used for logical indexing
    collegeId(i) = Rcpp::as<arma::mat>(collegeIdr[i]); 
    sopt2equ(i) = Rcpp::as<arma::umat>(sopt2equr[i]);
    equ2sopt(i) = Rcpp::as<arma::umat>(equ2soptr[i]);
    copt2equ(i) = Rcpp::as<arma::umat>(copt2equr[i]);
    equ2copt(i) = Rcpp::as<arma::umat>(equ2coptr[i]);
  }    
  
  for(int i=0; i<N; i++){
    studentIds(i) = Rcpp::as<arma::umat>(studentIdsr[i]);  // umat is used for logical indexing
  }
  
  arma::colvec nColleges = nCollegesr;
  arma::colvec nStudents = nStudentsr;
  arma::colvec nEquilibs = nEquilibsr;
  arma::colvec coptid = coptidr;
  
  // ---------------------------------------------
  // Outcome variables.
  // ---------------------------------------------

  arma::field<arma::colvec> Vc(T);
  for(int i=0; i<T; i++){
    Vc(i) = arma::zeros(C(i).n_rows,1);              // latent match valuation of colleges
  }
  arma::field<arma::colvec> Vs(T);
  for(int i=0; i<T; i++){
    Vs(i) = arma::zeros(S(i).n_rows,1);              // latent match valuation of students
  }
  
  // ---------------------------------------------
  // Priors.
  // ---------------------------------------------
  
  const int kC = C(0).n_cols; // number of parameters
  const int kS = S(0).n_cols; // number of parameters

  // beta
  arma::colvec betabar = arma::zeros(kC,1);
  arma::mat sigmabarbeta = 10*arma::eye(kC,kC);
  arma::mat sigmabarbetainverse = arma::inv(sigmabarbeta);    
  arma::colvec betabaroversigma = sigmabarbetainverse*betabar;

  // gamma
  arma::colvec gammabar = arma::zeros(kS,1);
  arma::mat sigmabargamma = 10*arma::eye(kS,kS);
  arma::mat sigmabargammainverse = arma::inv(sigmabargamma);    
  arma::colvec gammabaroversigma = sigmabargammainverse*gammabar;
  
  // ---------------------------------------------
  // Prior modes for parameters.
  // ---------------------------------------------

  // beta
  arma::colvec beta = betabar;
  arma::colvec betahat = arma::zeros(kC,1);
  arma::mat sigmahatbeta = arma::zeros(kC,kC);
  
  // gamma
  arma::colvec gamma = gammabar;
  arma::colvec gammahat = arma::zeros(kS,1);
  arma::mat sigmahatgamma = arma::zeros(kS,kS);
  
  // sums
  arma::mat sum3    = arma::zeros(kC,kC); 
  arma::colvec sum4 = arma::zeros(kC,1);
  arma::mat sum5    = arma::zeros(kS,kS); 
  arma::colvec sum6 = arma::zeros(kS,1);
  arma::colvec sum7 = arma::zeros(1,1); 
  arma::colvec sum8 = arma::zeros(1,1);
  
  // eta and delta
  arma::colvec eta(n);
  arma::colvec delta(n);
    
  // ---------------------------------------------
  // Matrices for parameter draws.
  // ---------------------------------------------
  
  arma::mat betadraws(kC,(niter-(niter % thin))/thin);
  arma::mat gammadraws(kS,(niter-(niter % thin))/thin);
  arma::mat etadraws(n,(niter-(niter % thin))/thin);
  arma::mat deltadraws(n,(niter-(niter % thin))/thin);

  // ---------------------------------------------
  // Main loop.
  // ---------------------------------------------  
  
  Rcpp::Rcout << "Drawing " << niter << " MCMC samples..." << std::endl;
  
  // Initiate Progress Bar
  Progress prog(niter, display_progress);
  
  for(int iter = 0; iter < niter; iter++){
    
    // update Progress Bar
    prog.increment(); 
    
    // ---------------------------------------------
    // Simulate each of the Vstar's conditional on all other Vstar's, data, and parameters.
    // ---------------------------------------------  
    
    double Vhat, Vcupperbar, Vsupperbar, Vclowerbar, Vslowerbar;
    arma::uvec iuvec;
    arma::vec ivec = arma::zeros(1,1);
    
    for(int t=0; t < T; t++){
      for(int i=0; i < nColleges(t); i++){
        // convert int i to uvec iuvec for non-contiguous indexing in umatrix M
        ivec(0) = i;  
        iuvec = arma::conv_to<arma::uvec>::from(ivec); 
        
        for(unsigned int j=0; j < nStudents(t); j++){
          
          // --- Calculation of equilibrium bounds ---
          
          Vcupperbar = arma::datum::inf;
          Vsupperbar = arma::datum::inf;
          Vclowerbar = -arma::datum::inf;
          Vslowerbar = -arma::datum::inf;
          
          // --- Draw of u_j,i: student j's valuation over college i (truncated by equilibrium bounds) ---
          
          for(int k=0; k < nEquilibs(t); k++){
          
            if(H(t)(i,j,k) == 0){  // non-equilibrium matches
              
              // Eqn (A4)
              // if student's valuation is higher than minimum of college's students (so college prefers student)
              if( Vs(t)(M(t)(i,j)) > arma::min(Vs(t)(M(t)(iuvec, studentIds(L(t)(i)).col(k) ))) ){  
                // then college's valuation has to be lower than that of student's equilibrium college
                Vcupperbar = std::min( Vcupperbar, Vc(t)(M(t)(collegeId(t)(j,k),j)) ); 
              }
              
            } else{  // equilibrium matches
              
              if(k == 0){ // student-optimal equilibrium
                
                for( int iprime=0; iprime < nColleges(t); iprime++ ){
                  // convert int iprime to uvec iprimeuvec for non-contiguous indexing in umatrix M
                  ivec(0) = iprime;  
                  arma::uvec iprimeuvec = arma::conv_to<arma::uvec>::from(ivec); 
              
                  // if student not matched to college iprime
                  if( H(t)(iprime,j,k) == 0 ){  
                    // Eqn (A6)
                    // if student's valuation is higher than minimum of imprime's students (so college prefers student)
                    if( Vs(t)(M(t)(iprime,j)) > arma::min(Vs(t)(M(t)(iprimeuvec,studentIds(L(t)(iprime)).col(k)))) ){  
                      // Eqn (A5), term 1
                      // then student must value college i higher than iprime
                      Vclowerbar = std::max( Vclowerbar, std::max(Vclowerbar,Vc(t)(M(t)(iprime,j))) ); 
                    }  
                  }     
                } 
                
                for(int l=1; l < nEquilibs(t); l++){
                  if( collegeId(t)(j,l) != collegeId(t)(j,0) ){ // only if mu'(j) != mu^s(j)
                    // Eqn (A5), term 2
                    Vclowerbar = std::max( Vclowerbar, Vc(t)(M(t)(collegeId(t)(j,l),j)) );
                  }
                }
                
              } else if(k == coptid[t]){ // college-optimal equilibrium
                
                for( int iprime=0; iprime < nColleges(t); iprime++ ){
                  // convert int iprime to uvec iprimeuvec for non-contiguous indexing in umatrix M
                  ivec(0) = iprime;  
                  arma::uvec iprimeuvec = arma::conv_to<arma::uvec>::from(ivec); 
                  
                  // if student not matched to college iprime
                  if( H(t)(iprime,j,k) == 0 ){  
                    // Eqn (A6)
                    // if student's valuation is higher than minimum of imprime's students (so college prefers student)
                    if( Vs(t)(M(t)(iprime,j)) > arma::min(Vs(t)(M(t)(iprimeuvec,studentIds(L(t)(iprime)).col(k)))) ){  
                      // Eqn (A10_b)
                      // then student must value college i higher than iprime
                      Vclowerbar = std::max( Vclowerbar, std::max(Vclowerbar,Vc(t)(M(t)(iprime,j))) ); 
                    }  
                  }     
                } 
                
                for(int l=0; l < nEquilibs(t); l++){
                  if( (collegeId(t)(j,l) != collegeId(t)(j,k)) ){ // only if mu'(j) != mu^s(j)
                    // Eqn (A11_b)
                    Vcupperbar = std::min( Vcupperbar, Vc(t)(M(t)(collegeId(t)(j,l),j)) );
                  }
                }

              } else{ // other equilibrium
                
                for( int iprime=0; iprime < nColleges(t); iprime++ ){
                  // convert int iprime to uvec iprimeuvec for non-contiguous indexing in umatrix M
                  ivec(0) = iprime;  
                  arma::uvec iprimeuvec = arma::conv_to<arma::uvec>::from(ivec); 
                
                  // if student not matched to college iprime
                  if( H(t)(iprime,j,k) == 0 ){  
                    // Eqn (A6)
                    // if student's valuation is higher than minimum of imprime's students (so college prefers student)
                    if( Vs(t)(M(t)(iprime,j)) > arma::min(Vs(t)(M(t)(iprimeuvec,studentIds(L(t)(iprime)).col(k)))) ){  
                      // Eqn (A11) = Eqn (A15_b), term 1
                      // then student must value college i higher than iprime
                      Vclowerbar = std::max( Vclowerbar, std::max(Vclowerbar,Vc(t)(M(t)(iprime,j))) ); 
                    }  
                  }     
                }
                
                if( collegeId(t)(j,k) != collegeId(t)(j,coptid[t]) ){ // only if mu'(j) != mu^s(j)
                  // Eqn (A15_b), term 2
                  Vclowerbar = std::max( Vclowerbar, Vc(t)(M(t)(collegeId(t)(j,coptid[t]),j)) );
                }
                
                if( collegeId(t)(j,k) != collegeId(t)(j,0) ){ // only if mu'(j) != mu^s(j)
                  // Eqn (A12) = Eqn (A16)
                  Vcupperbar = std::min( Vcupperbar, Vc(t)(M(t)(collegeId(t)(j,0),j)) );
                }
                
              } // student-optimal vs. other equilibrium         
            } // non-equilibrium vs. equilibrium 
          } // k = nEquilibs
          
          if((Vclowerbar != 0) | (Vcupperbar != 0)){ 
            
            Vhat = arma::as_scalar( C(t).row(M(t)(i,j))*beta ); 
            Vc(t)(M(t)(i,j)) = truncn2(Vhat, 1.0, Vclowerbar, Vcupperbar); // mu, sigma, lower, upper
          }
          
          // --- Draw of u_i,j: college i's valuation over student j (truncated by equilibrium bounds) ---
          
          for(int k=0; k < nEquilibs(t); k++){
            
            if(H(t)(i,j,k) == 0){  // non-equilibrium matches
              
              // Eqn (A3)            
              // if student prefers college to current match
              if( Vc(t)(M(t)(i,j)) > Vc(t)(M(t)(collegeId(t)(j,k),j)) ){ 
                // then students' valuation has to be lower than that of worst student attending the college
                Vsupperbar = std::min( Vsupperbar, arma::min(Vs(t)(M(t)( iuvec,studentIds(L(t)(i)).col(k) ))) ); 
              }
              
            } else{  // equilibrium matches
              
              if(k==0){ // student-optimal equilibrium
                
                for( int jprime=0; jprime < nStudents(t); jprime++ ){             
                  // if college not matched to student jprime
                  if( H(t)(i,jprime,k)==0 ){  
                    // Eqn (A9)
                    // but student jprime values i over his equilibrium school
                    if( Vc(t)(M(t)(i,jprime)) > Vc(t)(M(t)(collegeId(t)(jprime,k),jprime)) ){ 
                      // Eqn (A7)
                      // then school must value student j higher than jprime  
                      Vslowerbar = std::max( Vslowerbar, Vs(t)(M(t)(i,jprime)) );
                    } 
                  } 
                }

                for(int l=1; l < nEquilibs(t); l++){
                  if( j != sopt2equ(t)(j,l) ){ // only if mu^s(i) != mu'(i)
                    // Eqn (A8)
                    Vsupperbar = std::min( Vsupperbar, Vs(t)(M(t)(i,sopt2equ(t)(j,l))) );
                  }
                }
                
              } else if(k == coptid[t]){ // college-optimal equilibrium

                for( int jprime=0; jprime < nStudents(t); jprime++ ){             
                  // if college not matched to student jprime
                  if( H(t)(i,jprime,k)==0 ){  
                    // Eqn (A9)
                    // but student jprime values i over his equilibrium school
                    if( Vc(t)(M(t)(i,jprime)) > Vc(t)(M(t)(collegeId(t)(jprime,k),jprime)) ){ 
                      // Eqn (A12_b), term 1
                      // then school must value student j higher than jprime  
                      Vslowerbar = std::max( Vslowerbar, Vs(t)(M(t)(i,jprime)) );
                    } 
                  } 
                }
                
                for(int l=0; l < nEquilibs(t); l++){ 
                  if( (j != copt2equ(t)(j,l)) & (l != k) ){ // only if mu^s(i) != mu'(i)
                    // Eqn (A12_b), term 2
                    Vslowerbar = std::max( Vslowerbar, Vs(t)(M(t)(i,copt2equ(t)(j,l))) );
                  }
                }
                
              } else{ // other equilibrium
                
                for( int jprime=0; jprime < nStudents(t); jprime++ ){             
                  // if college not matched to student jprime
                  if( H(t)(i,jprime,k)==0 ){  
                    // Eqn (A9)
                    // but student jprime values i over his equilibrium school
                    if( Vc(t)(M(t)(i,jprime)) > Vc(t)(M(t)(collegeId(t)(jprime,k),jprime)) ){ 
                      // Eqn (A10), term 1
                      // then school must value student j higher than jprime  
                      Vslowerbar = std::max( Vslowerbar, Vs(t)(M(t)(i,jprime)) );
                    } 
                  } 
                }

                if( j != equ2sopt(t)(j,k) ){ // only if mu'(i) != mu^s(i)
                  // Eqn (A10), term 2
                  Vslowerbar = std::max( Vslowerbar, Vs(t)(M(t)(i,equ2sopt(t)(j,k))) );
                }
                
                if( j != equ2copt(t)(j,k) ){ // only if mu'(i) != mu^s(i)
                  // Eqn (A14_b)
                  Vsupperbar = std::min( Vsupperbar, Vs(t)(M(t)(i,equ2copt(t)(j,k))) );
                }
                
              } // student-optimal vs. other equilibrium         
            } // non-equilibrium vs. equilibrium 
          } // k = nEquilibs
          
          if((Vslowerbar != 0) | (Vsupperbar != 0)){ 
            
            Vhat = arma::as_scalar( S(t).row(M(t)(i,j))*gamma ); 
            Vs(t)(M(t)(i,j)) = truncn2(Vhat, 1.0, Vslowerbar, Vsupperbar); // mu, sigma, lower, upper
          }

        } // j = nStudents
      } // i = nColleges
    } // t = T
    
    // ---------------------------------------------
    // Simulate each group of parameters conditional on latents, data, and all other parameters.
    // ---------------------------------------------

    // ---------------------------------------------
    // beta.
    // ---------------------------------------------

    sum3.zeros(kC,kC); // reset to zero
    sum4.zeros(kC,1);  // reset to zero

    for(int t=0; t<T; t++){ 
      sum3 += CC(t);
      sum4 += -trans(C(t))*Vc(t);
    } 
    sigmahatbeta = arma::inv(sigmabarbetainverse + sum3); 
    betahat = -sigmahatbeta * (-betabaroversigma + sum4); 
    beta = mvrnormArma(betahat,sigmahatbeta,kC);
    
    // ---------------------------------------------
    // gamma.
    // ---------------------------------------------
     
    sum5.zeros(kS,kS); // reset to zero
    sum6.zeros(kS,1);  // reset to zero
    
    for(int t=0; t<T; t++){ 
      sum5 += SS(t);
      sum6 += -trans(S(t))*Vs(t);
    } 
    sigmahatgamma = arma::inv(sigmabargammainverse + sum5); 
    gammahat = -sigmahatgamma * (-gammabaroversigma + sum6); 
    gamma = mvrnormArma(gammahat,sigmahatgamma,kS);
    
    // ---------------------------------------------
    // eta.
    // ---------------------------------------------

    double etacount = 0;
    for(int t=0; t<T; t++){
      eta.rows(etacount, etacount + nStudents(t) - 1) = Vc(t)(d(t)) - Cmatch(t)*beta;
      etacount = etacount + nStudents(t);
    }
    
    // ---------------------------------------------
    // delta.
    // ---------------------------------------------
    
    double deltacount = 0;
    for(int t=0; t<T; t++){
      delta.rows(deltacount, deltacount + nStudents(t) - 1) = Vs(t)(d(t)) - Smatch(t)*gamma;
      deltacount = deltacount + nStudents(t);
    }

    // ---------------------------------------------  
    // save this iterations draws.
    // ---------------------------------------------
    
    if((iter % thin == 0) & (iter<(niter-thin+1))){
      betadraws.col(iter/thin) = beta;
      gammadraws.col(iter/thin) = gamma;
      etadraws.col(iter/thin) = eta;
      deltadraws.col(iter/thin) = delta;
    }
  }
  
  // ---------------------------------------------  
  // Return the parameter draws.
  // ---------------------------------------------  

  return List::create(  
    // parameter draws
    Named("betadraws") = betadraws,
    Named("gammadraws") = gammadraws,
    Named("etadraws") = etadraws,
    Named("deltadraws") = deltadraws
  );  
}




