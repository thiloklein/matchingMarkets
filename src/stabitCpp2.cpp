//#define ARMA_NO_DEBUG
// Uncomment the above line if you want to disable all run-time checks.
// This will result in faster code, but you first need to make sure that your code runs correctly!

// [[Rcpp::depends(RcppArmadillo,RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>

using namespace Rcpp;
arma::colvec mvrnormArma(arma::colvec mu, arma::mat sigma, int ncols);
double truncn(double bound, bool lb, double mu, double sigma);
double truncn2(double a, double b, double mu, double sigma);
double f(double y);

//' @export
// [[Rcpp::export]]
List stabitCpp2(Rcpp::List Yr, Rcpp::List Xmatchr, Rcpp::List Cr, 
  Rcpp::List Cmatchr, Rcpp::List Dr, Rcpp::List dr, Rcpp::List Mr, Rcpp::List Hr, 
  arma::colvec nCollegesr, arma::colvec nStudentsr, Rcpp::List XXmatchr,
  Rcpp::List CCr, Rcpp::List CCmatchr, 
  Rcpp::List Lr, Rcpp::List studentIdsr, Rcpp::List collegeIdr, int n, int N,
  bool binary, int niter, int T, int censored, bool display_progress = true) {
  
  //-- 
  // bool gPrior, arma::mat sigmabarbetainverse, arma::mat sigmabaralphainverse, arma::mat sigmabargammainverse,
  //int censored = 1;
  
  // Enable/Disable verbose debug tracing.
  //bool DEBUG = FALSE;

  // ---------------------------------------------
  // Independent variables
  // ---------------------------------------------
    
  arma::field<arma::mat> Xmatch(T), C(T), Cmatch(T), H(T), XXmatch(T), CC(T), CCmatch(T); 
  arma::field<arma::uvec> L(T), studentIds(N), d(T);
  arma::field<arma::umat> M(T);
  arma::field<arma::colvec> D(T), collegeId(T);

  for(int i=0; i<T; i++){
    Xmatch(i) = Rcpp::as<arma::mat>(Xmatchr[i]);  
    C(i) = Rcpp::as<arma::mat>(Cr[i]);  
    Cmatch(i) = Rcpp::as<arma::mat>(Cmatchr[i]);  
    D(i) = Rcpp::as<arma::colvec>(Dr[i]); 
    d(i) = Rcpp::as<arma::uvec>(dr[i]);  // uvec is used for logical indexing
    M(i) = Rcpp::as<arma::umat>(Mr[i]);  // umat is used for logical indexing
    H(i) = Rcpp::as<arma::mat>(Hr[i]);  
    XXmatch(i) = Rcpp::as<arma::mat>(XXmatchr[i]); 
    CC(i) = Rcpp::as<arma::mat>(CCr[i]);  
    CCmatch(i) = Rcpp::as<arma::mat>(CCmatchr[i]); 
    L(i) = Rcpp::as<arma::uvec>(Lr[i]);  // uvec is used for logical indexing
    collegeId(i) = Rcpp::as<arma::colvec>(collegeIdr[i]); 
  }    
  
  for(int i=0; i<N; i++){
    studentIds(i) = Rcpp::as<arma::uvec>(studentIdsr[i]);  // uvec is used for logical indexing
  }
  
  arma::colvec nColleges = nCollegesr;
  arma::colvec nStudents = nStudentsr;
  
  // ---------------------------------------------
  // Outcome variables.
  // ---------------------------------------------
  
  arma::field<arma::colvec> R(T), Y(T);
  for(int i=0; i<T; i++){
    R(i)      = Rcpp::as<arma::colvec>(Yr[i]);      // dependent variable
    if(binary == TRUE){
      Y(i)    = arma::zeros(R(i).n_rows,1); // latent outcome (set to zero)
    } else{
      Y(i)    = R(i);          
    }
  }
  
  arma::field<arma::colvec> Vc(T);
  for(int i=0; i<T; i++){
    Vc(i) = arma::zeros(C(i).n_rows,1);              // latent match valuation of colleges
  }
  
  // ---------------------------------------------
  // Priors.
  // ---------------------------------------------
  
  const int kX = Xmatch(0).n_cols; // number of parameters
  const int kC = C(0).n_cols; // number of parameters
  
  // alpha
  arma::colvec alphabar = arma::zeros(kX,1);
  //if(gPrior == FALSE){
    arma::mat sigmabaralpha = 10*arma::eye(kX,kX);
    arma::mat sigmabaralphainverse = arma::inv(sigmabaralpha);
  //} else{
  //  sigmabaralphainverse = sigmabaralphainverse;
  //}
  arma::mat alphabaroversigma = sigmabaralphainverse*alphabar;  

  // beta
  arma::colvec betabar = arma::zeros(kC,1);
  //if(gPrior == FALSE){
    arma::mat sigmabarbeta = 10*arma::eye(kC,kC);
    arma::mat sigmabarbetainverse = arma::inv(sigmabarbeta);    
  //} else{
  //  sigmabarbetainverse = sigmabarbetainverse;
  //}
  arma::colvec betabaroversigma = sigmabarbetainverse*betabar;
  
  // kappa
  double kappabar = 0;
  double sigmabarsquarekappa = 10;
    
  // ---------------------------------------------
  // Prior modes for parameters.
  // ---------------------------------------------

  // alpha
  arma::colvec alpha = alphabar;
  arma::colvec alphahat = arma::zeros(kX,1);
  arma::mat sigmahatalpha = arma::zeros(kX,kX);

  // beta
  arma::colvec beta = betabar;
  arma::colvec betahat = arma::zeros(kC,1);
  arma::mat sigmahatbeta = arma::zeros(kC,kC);
    
  // kappa
  double kappa = kappabar;
  double kappahat = 0;
  double sigmahatsquarekappa = 0;
    
  // nu
  double a = 2, b = 1, ahat, bhat;
  double sigmasquarenuinverse = b*(a-1);
  double sigmasquarenu = 1/sigmasquarenuinverse;
  
  // sums
  arma::mat sum1    = arma::zeros(kX,kX); 
  arma::colvec sum2 = arma::zeros(kX,1);
  arma::mat sum3    = arma::zeros(kC,kC); 
  arma::colvec sum4 = arma::zeros(kC,1);
  arma::colvec sum7 = arma::zeros(1,1); 
  arma::colvec sum8 = arma::zeros(1,1);
  
  // eta and delta
  arma::colvec eta(n);
  //arma::colvec delta(n);
    
  // ---------------------------------------------
  // Matrices for parameter draws.
  // ---------------------------------------------
  
  arma::mat alphadraws(kX,niter);
  arma::mat betadraws(kC,niter);
  arma::mat kappadraws(1,niter);
  arma::mat etadraws(n,niter);
  //arma::mat deltadraws(n,niter);
  arma::mat sigmasquarenudraws(1,niter);
  
  // ---------------------------------------------
  // Main loop.
  // ---------------------------------------------  
  
  Rcout << "Drawing " << niter << " MCMC samples..." << std::endl;
  
  // Initiate Progress Bar
  Progress prog(niter, display_progress);
  
  for(int iter = 0; iter < niter; iter++){
    
    //if(iter % 1000 == 999){
    //  Rcout << iter+1 << " of " << niter << std::endl;
    //}
    
    // update Progress Bar
    prog.increment(); 
    
    // ---------------------------------------------
    // Simulate the latent outcomes conditional on the latent group valuations, data, and parameters.
    // ---------------------------------------------
    
    double Yhat; 
    
    if(binary == TRUE){
      
      for(int t=0; t < T; t++){ // two-group-market identifiers.        
        for(int ij=0; ij < nStudents(t); ij++){
          Yhat = arma::as_scalar( Xmatch(t).row(ij)*alpha + ( Vc(t)(ij) - Cmatch(t).row(ij)*beta )*kappa );
          if(R(t)(ij) == 1){
            Y(t)(ij) = truncn(0, TRUE, Yhat, 1.0); // bound, lb, mu, sigma
          } else{
            Y(t)(ij) = truncn(0, FALSE, Yhat, 1.0); // bound, lb, mu, sigma
          }            
        }
      }                    
    }
    
    // ---------------------------------------------
    // Simulate each of the Vstar's conditional on all other Vstar's, data, and parameters.
    // ---------------------------------------------
      
    double Vhat, sigmahatsquareV, Vclowerbar, Vcupperbar;
    arma::uvec iuvec, iprimeuvec;
    arma::vec ivec = arma::zeros(1,1);
    
    for(int t=0; t < T; t++){
      
      for(int i=0; i < nColleges(t); i++){
        // convert int i to uvec iuvec for non-contiguous indexing in umatrix M
        ivec(0) = i;  
        iuvec = arma::conv_to<arma::uvec>::from(ivec); 
        
        for(int j=0; j < nStudents(t); j++){
          
          // --- Calculation of equilibrium bounds ---
          
          //Vcupperbar = arma::datum::inf;
          //Vsupperbar = arma::datum::inf;
          Vclowerbar = -arma::datum::inf;
          //Vslowerbar = -arma::datum::inf;
          
          if(H(t)(i,j) == 0){  // non-equilibrium matches
            
            // --- Draw valuation (truncated by upper equilibrium bound) ---
                   
            Vhat = arma::as_scalar( C(t).row(M(t)(i,j))*beta );            
            // Eqn (5)
            Vcupperbar = std::max( Vc(t)(M(t)(collegeId(t)(j),j)), arma::min(Vc(t)(M(t)(iuvec,studentIds(L(t)(i))))) ); 
            Vc(t)(M(t)(i,j)) = truncn(Vcupperbar, FALSE, Vhat, 1.0); // bound, lb, mu, sigma

          } else{  // equilibrium matches
                          
            // --- Draw valuation (truncated by lower equilibrium bounds) ---
                        
            // Eqn (7)
            for( int jprime=0; jprime < nStudents(t); jprime++ ){              
              if( H(t)(i,jprime)==0 ){  // college not matched to student jprime
                if( Vc(t)(M(t)(i,jprime)) > Vc(t)(M(t)(collegeId(t)(jprime),jprime)) ){
                  Vclowerbar = std::max(Vclowerbar, Vc(t)(M(t)(i,jprime)));                  
                } 
              } 
            } 
            
            // Eqn (8)
            for( int iprime=0; iprime < nColleges(t); iprime++ ){
              // convert int iprime to uvec iprimeuvec for non-contiguous indexing in umatrix M
              ivec(0) = iprime;  
              iprimeuvec = arma::conv_to<arma::uvec>::from(ivec); 
              
              if( H(t)(iprime,j) == 0 ){  // student not matched to college iprime
                if( Vc(t)(M(t)(iprime,j)) > arma::min(Vc(t)(M(t)(iprimeuvec,studentIds(L(t)(iprime))))) ){  // student's valuation is higher than minimum of imprime's students (so college prefers student)
                  Vclowerbar = std::max(Vclowerbar, Vc(t)(M(t)(iprime,j)));
                }  
              }     
            } 
            
            sum7 = Y(t)(M(t)(i,j)) - Xmatch(t).row(M(t)(i,j))*alpha;
            Vhat = arma::as_scalar( C(t).row(M(t)(i,j))*beta + kappa*sum7/(sigmasquarenu+pow(kappa,2)) );
            sigmahatsquareV = sigmasquarenu/(sigmasquarenu+pow(kappa,2));
            Vc(t)(M(t)(i,j)) = truncn(Vclowerbar, TRUE, Vhat, sqrt(sigmahatsquareV)); // bound, lb, mu, sigma
            
          }  
          
        }
      } 
    }
    
    // ---------------------------------------------
    // Simulate each group of parameters conditional on latents, data, and all other parameters.
    // ---------------------------------------------
  
    // ---------------------------------------------
    // alpha.
    // ---------------------------------------------
    
    sum1.zeros(kX,kX); // reset to zero
    sum2.zeros(kX,1);  // reset to zero
    
    for(int t = 0; t < T; t++){
      sum1 += XXmatch(t);
      sum2 += trans(Xmatch(t))*( Y(t) - kappa*(Vc(t)(d(t)) - Cmatch(t)*beta));
    }          
    sigmahatalpha = arma::inv(sigmabaralphainverse + (1/sigmasquarenu)*sum1);
    alphahat = -sigmahatalpha * (-alphabaroversigma - (1/sigmasquarenu)*sum2);
    alpha = mvrnormArma(alphahat,sigmahatalpha,kX); 
        
    // ---------------------------------------------
    // beta.
    // ---------------------------------------------
    
    sum3.zeros(kC,kC); // reset to zero
    sum4.zeros(kC,1);  // reset to zero
    
    for(int t=0; t<T; t++){ 
      sum3 += CC(t) + (pow(kappa,2)/sigmasquarenu)*CCmatch(t);
      sum4 += -trans(C(t))*Vc(t) + (kappa/sigmasquarenu)*(trans(Cmatch(t))*(Y(t) - Xmatch(t)*alpha - kappa*Vc(t)(d(t))) );
    } 
    sigmahatbeta = arma::inv(sigmabarbetainverse + sum3); 
    betahat = -sigmahatbeta * (-betabaroversigma + sum4); 
    beta = mvrnormArma(betahat,sigmahatbeta,kC);
    
    // ---------------------------------------------
    // kappa.
    // ---------------------------------------------
      
    sum7.zeros(1,1); // reset to zero
    sum8.zeros(1,1); // reset to zero
    
    for(int t=0; t<T; t++){ 
      //sum7 += sum( pow(Vc(t)(d(t)) - Cmatch(t)*beta, 2) );      
      //sum8 += sum( (Y(t) - Xmatch(t)*alpha - lambda*(Vs(t)(d(t)) - Smatch(t)*gamma)) * (Vc(t)(d(t)) - Cmatch(t)*beta) );
      sum7 += trans(Vc(t)(d(t)) - Cmatch(t)*beta) * (Vc(t)(d(t)) - Cmatch(t)*beta);      
      sum8 += trans(Y(t) - Xmatch(t)*alpha) * (Vc(t)(d(t)) - Cmatch(t)*beta);
    }      
    sigmahatsquarekappa = 1/(1/sigmabarsquarekappa + (1/sigmasquarenu)*arma::as_scalar(sum7));
    kappahat = -sigmahatsquarekappa*(-kappabar/sigmabarsquarekappa - (1/sigmasquarenu)*arma::as_scalar(sum8));      
    if(censored == 1){ // from below, positive covariation of residuals
      kappa = truncn(0, TRUE, kappahat, sqrt(sigmahatsquarekappa) ); // bound, lb, mu, sigma
    } else if(censored == 2){ // from above, negative covariation of residuals
      kappa = truncn(0, FALSE, kappahat, sqrt(sigmahatsquarekappa) ); // bound, lb, mu, sigma
    } else{ // not censored
      kappa = ::Rf_rnorm(kappahat, sqrt(sigmahatsquarekappa));
    }
    
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
    //for(int t=0; t<T; t++){
    //  delta.rows(2*t,2*t+1) = V(t).rows(0,1) - W(t).rows(0,1)*alpha;
    //}
    
    // ---------------------------------------------
    // sigma nu
    // ---------------------------------------------
    
    if(binary == FALSE){
      
      sum7.zeros(1,1); // reset to zero
      
      for(int t = 0; t < T; t++){
        sum7 += sum( pow(Y(t) - Xmatch(t)*alpha - kappa*(Vc(t)(d(t)) - Cmatch(t)*beta), 2) );
      }    
      
      ahat = a + n/2;          // shape
      bhat = 1/(1/b + arma::as_scalar(sum7)/2); // scale = 1/rate
      sigmasquarenuinverse = ::Rf_rgamma(ahat, bhat); // shape, scale
      sigmasquarenu = 1/sigmasquarenuinverse;
      
    }
    
    // ---------------------------------------------  
    // save this iterations draws.
    // ---------------------------------------------
    
    alphadraws.col(iter) = alpha;
    betadraws.col(iter) = beta;
    kappadraws.col(iter) = kappa;
    etadraws.col(iter) = eta;
    //deltadraws.col(iter) = delta;
    
    if(binary == FALSE){
      sigmasquarenudraws.col(iter) = sigmasquarenu;
    }
  }
  
  // print to screen
  //Rcout << "done." << std::endl;
  //Rcout << std::endl;
  
  // ---------------------------------------------  
  // The last half of all draws are used in approximating the posterior means and the posterior standard deviations.
  // ---------------------------------------------  
  
  int startiter = niter/2;
  
  if(binary == TRUE){
    return List::create(  
      // parameter draws
      Named("alphadraws") = betadraws,
      Named("betadraws") = alphadraws,
      Named("deltadraws") = kappadraws,
      // posterior means
      Named("eta") = mean(etadraws.cols(startiter,niter-1),1),
      //Named("delta") = mean(deltadraws.cols(startiter,niter-1),1),
      Named("alpha") = join_rows( mean(betadraws.cols(startiter,niter-1),1), stddev(betadraws.cols(startiter,niter-1),0,1) ),
      Named("beta") = join_rows( mean(alphadraws.cols(startiter,niter-1),1), stddev(alphadraws.cols(startiter,niter-1),0,1) ),
      Named("delta") = join_rows( mean(kappadraws.cols(startiter,niter-1),1), stddev(kappadraws.cols(startiter,niter-1),0,1) ),
      Named("sigmasquarexi") = join_rows( arma::ones(1,1) , arma::zeros(1,1) ),
      // vcov
      Named("alphavcov") = cov(trans(betadraws.cols(startiter,niter-1))),
      Named("betavcov") = cov(trans(alphadraws.cols(startiter,niter-1)))
    );  
  } else if(binary == FALSE){
    return List::create(  
      // parameter draws
      Named("alphadraws") = betadraws,
      Named("betadraws") = alphadraws,
      Named("deltadraws") = kappadraws,
      Named("sigmasquarexidraws") = sigmasquarenudraws,
      // posterior means
      Named("eta") = mean(etadraws.cols(startiter,niter-1),1),
      //Named("delta") = mean(deltadraws.cols(startiter,niter-1),1),
      Named("alpha") = join_rows( mean(betadraws.cols(startiter,niter-1),1), stddev(betadraws.cols(startiter,niter-1),0,1) ),
      Named("beta") = join_rows( mean(alphadraws.cols(startiter,niter-1),1), stddev(alphadraws.cols(startiter,niter-1),0,1) ),
      Named("delta") = join_rows( mean(kappadraws.cols(startiter,niter-1),1), stddev(kappadraws.cols(startiter,niter-1),0,1) ),
      Named("sigmasquarexi") = join_rows( mean(sigmasquarenudraws.cols(startiter,niter-1),1), stddev(sigmasquarenudraws.cols(startiter,niter-1),0,1) ),
      // vcov
      Named("alphavcov") = cov(trans(betadraws.cols(startiter,niter-1))),
      Named("betavcov") = cov(trans(alphadraws.cols(startiter,niter-1)))
    );
  } else{
      return 0;
  }
}


// ---------------------------------------------  
// Returns one draw from truncated normal distribution (mu,sigma^2) with range
// (a,b) from http://athens.src.uchicago.edu/jenni/econ319_2003/lecture.html
// ---------------------------------------------  

double f(double y){
  return exp(-0.5*pow(y,2));
}

double truncn2(double a, double b, double mu, double sigma){
  
  arma::colvec u = arma::zeros(2,1);
  double c, c1, c2, x, cdel, f1, f2, z, dtruncn;
  double eps=1e-20, t1=0.375, t2=2.18, t3=0.725, t4=0.45;
  bool lflip;
  //double c, z, w;
  
  // 1. standardised cut-off c for truncation from below or above 
  c1 = (a-mu)/sigma;
  c2 = (b-mu)/sigma;
  lflip = FALSE;
    
  // 2. standardised draw using Geweke's (1991) normal rejection sampling
  if(c1*c2 < 0){
    if((f(c1)>t1) & (f(c2)>t1)){
      cdel = c2 - c1;
      
      u(0) = ::Rf_runif(0.0,1.0); // u = runif(2);
      u(1) = ::Rf_runif(0.0,1.0); 
      x = c1 + cdel*u(0);
      while(-(u(1) < f(x))){
        u(0) = ::Rf_runif(0.0,1.0); // u = runif(2);
        u(1) = ::Rf_runif(0.0,1.0); 
        x = c1 + cdel*u(0);
      }
    } else{
      
      x = ::Rf_rnorm(0.0,1.0);
      while(-((x>c1) & (x<c2))){
        x = ::Rf_rnorm(0.0,1.0);
      }
    }
  } else{
    
    if(c1<0){
      c  =  c1;
      c1 = -c2;
      c2 = -c;
      lflip = TRUE;
    }
    f1 = f(c1);
    f2 = f(c2);
    
    if((f2<eps) | (f1/f2>t2)){
      if(c1>t3){  // exponential rejection sampling
        
        c = c2 - c1;
        
        u(0) = ::Rf_runif(0.0,1.0);  // u = runif(2);
        u(1) = ::Rf_runif(0.0,1.0);
        z = -log(u(0))/c1;
        while(-((z<c) & (u(1)<f(z)))){
          u(0) = ::Rf_runif(0.0,1.0);  // u = runif(2);
          u(1) = ::Rf_runif(0.0,1.0);
          z = -log(u(0))/c1;
        }
        x = c1 + z;
        
      } else{  // half-normal rejection sampling
        
        x = ::Rf_rnorm(0.0,1.0);
        x = std::abs(x);
        while(-((x>c1) & (x<c2))){
          x = ::Rf_rnorm(0.0,1.0);
          x = std::abs(x);
        }
      }
    } else{  // uniform rejection sampling
      cdel = c2 - c1;
      
      u(0) = ::Rf_runif(0.0,1.0);  // u = runif(2);
      u(1) = ::Rf_runif(0.0,1.0);
      x = c1 + cdel*u(0);
      while(-(u(1) < (f(x)/f1))){
        u(0) = ::Rf_runif(0.0,1.0);  // u = runif(2);
        u(1) = ::Rf_runif(0.0,1.0);
        x = c1 + cdel*u(0);
      }
    }
  }
  
  // 3. reverse standardisation
  if(lflip){
    return dtruncn = mu - sigma*x;
  } else{
    return dtruncn = mu + sigma*x;
  }
}


