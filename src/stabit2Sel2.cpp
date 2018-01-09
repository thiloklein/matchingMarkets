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

// [[Rcpp::export]]
List stabit2Sel2(Rcpp::List Yr, Rcpp::List Xmatchr, Rcpp::List Cr, 
  Rcpp::List Cmatchr, Rcpp::List Dr, Rcpp::List dr, Rcpp::List Mr, Rcpp::List Hr, 
  arma::colvec nCollegesr, arma::colvec nStudentsr, Rcpp::List XXmatchr,
  Rcpp::List CCr, Rcpp::List CCmatchr, 
  Rcpp::List Lr, Rcpp::List studentIdsr, Rcpp::List collegeIdr, int n, int N,
  bool binary, int niter, int T, int censored, int thin, bool display_progress = true) {
  
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
  arma::mat sigmabaralpha = 10*arma::eye(kX,kX);
  arma::mat sigmabaralphainverse = arma::inv(sigmabaralpha);
  arma::mat alphabaroversigma = sigmabaralphainverse*alphabar;  

  // beta
  arma::colvec betabar = arma::zeros(kC,1);
  arma::mat sigmabarbeta = 10*arma::eye(kC,kC);
  arma::mat sigmabarbetainverse = arma::inv(sigmabarbeta);    
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
  
  // eta
  arma::colvec eta(n);
  
  // ---------------------------------------------
  // Matrices for parameter draws.
  // ---------------------------------------------
  
  arma::mat alphadraws(kX,(niter-(niter % thin))/thin); //floor(niter/thin);
  arma::mat betadraws(kC,(niter-(niter % thin))/thin);
  arma::mat kappadraws(1,(niter-(niter % thin))/thin);
  arma::mat etadraws(n,(niter-(niter % thin))/thin);
  arma::mat sigmasquarenudraws(1,(niter-(niter % thin))/thin);
  
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
    // Simulate the latent outcomes conditional on the latent group valuations, data, and parameters.
    // ---------------------------------------------
    
    double Yhat; 
    
    if(binary == TRUE){
      
      for(int t=0; t < T; t++){ // two-group-market identifiers.        
        for(int ij=0; ij < nStudents(t); ij++){
          Yhat = arma::as_scalar( Xmatch(t).row(ij)*alpha + ( Vc(t)(ij) - Cmatch(t).row(ij)*beta )*kappa );
          if(R(t)(ij) == 1){
            Y(t)(ij) = truncn2(Yhat, 1.0, 0, arma::datum::inf ); // mu, sigma, lower, upper
          } else{
            Y(t)(ij) = truncn2(Yhat, 1.0, -arma::datum::inf, 0 ); // mu, sigma, lower, upper
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
        
        for(unsigned int j=0; j < nStudents(t); j++){
          
          // --- Calculation of equilibrium bounds ---
          
          Vclowerbar = -arma::datum::inf;

          if(H(t)(i,j) == 0){  // non-equilibrium matches
            
            // --- Draw valuation (truncated by upper equilibrium bound) ---
                   
            Vhat = arma::as_scalar( C(t).row(M(t)(i,j))*beta );            
            // Eqn (5)
            Vcupperbar = std::max( Vc(t)(M(t)(collegeId(t)(j),j)), arma::min(Vc(t)(M(t)(iuvec,studentIds(L(t)(i))))) ); 
            Vc(t)(M(t)(i,j)) = truncn2(Vhat, 1.0, -arma::datum::inf, Vcupperbar ); // mu, sigma, lower, upper

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
            Vc(t)(M(t)(i,j)) = truncn2(Vhat, sqrt(sigmahatsquareV), Vclowerbar, arma::datum::inf ); // mu, sigma, lower, upper
            
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
      sum7 += trans(Vc(t)(d(t)) - Cmatch(t)*beta) * (Vc(t)(d(t)) - Cmatch(t)*beta);      
      sum8 += trans(Y(t) - Xmatch(t)*alpha) * (Vc(t)(d(t)) - Cmatch(t)*beta);
    }      
    sigmahatsquarekappa = 1/(1/sigmabarsquarekappa + (1/sigmasquarenu)*arma::as_scalar(sum7));
    kappahat = -sigmahatsquarekappa*(-kappabar/sigmabarsquarekappa - (1/sigmasquarenu)*arma::as_scalar(sum8));      
    if(censored == 1){ // from below, positive covariation of residuals
      kappa = truncn2(kappahat, sqrt(sigmahatsquarekappa), 0, arma::datum::inf ); // mu, sigma, lower, upper
    } else if(censored == 2){ // from above, negative covariation of residuals
      kappa = truncn2(kappahat, sqrt(sigmahatsquarekappa), -arma::datum::inf, 0 ); // mu, sigma, lower, upper
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

    if((iter % thin == 0) & (iter<(niter-thin+1))){
      alphadraws.col(iter/thin) = alpha;
      betadraws.col(iter/thin) = beta;
      kappadraws.col(iter/thin) = kappa;
      etadraws.col(iter/thin) = eta;
      
      if(binary == FALSE){
        sigmasquarenudraws.col(iter/thin) = sigmasquarenu;
      }
    }
  }
  
  // ---------------------------------------------  
  // Return the parameter draws.
  // ---------------------------------------------  
  
  if(binary == TRUE){
    return List::create(  
      // parameter draws
      Named("alphadraws") = alphadraws,
      Named("betadraws") = betadraws,
      Named("kappadraws") = kappadraws,
      Named("etadraws") = etadraws
    );  
  } else if(binary == FALSE){
    return List::create(  
      // parameter draws
      Named("alphadraws") = alphadraws,
      Named("betadraws") = betadraws,
      Named("kappadraws") = kappadraws,
      Named("etadraws") = etadraws,
      Named("sigmasquarenudraws") = sigmasquarenudraws
    );
  } else{
      return 0;
  }
}


// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================

double norm_rs(double a, double b)
{
   double  x;
   x = Rf_rnorm(0.0, 1.0);
   while( (x < a) || (x > b) ) x = norm_rand();
   return x;
}

// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// ======================================================================

double half_norm_rs(double a, double b)
{
   double   x;
   x = fabs(norm_rand());
   while( (x<a) || (x>b) ) x = fabs(norm_rand());
   return x;
}

// unif_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using uniform rejection sampling. 
// ======================================================================

double unif_rs(double a, double b)
{
   double xstar, logphixstar, x, logu;

   // Find the argmax (b is always >= 0)
   // This works because we want to sample from N(0,1)
   if(a <= 0.0) xstar = 0.0;
   else xstar = a;
   logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);

   x = R::runif(a, b);
   logu = log(R::runif(0.0, 1.0));
   while( logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
   {
      x = R::runif(a, b);
      logu = log(R::runif(0.0, 1.0));
   }
   return x;
}

// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// ======================================================================

double exp_rs(double a, double b)
{
  double  z, u, rate;

  rate = 1/a;

   // Generate a proposal on (0, b-a)
   z = R::rexp(rate);
   while(z > (b-a)) z = R::rexp(rate);
   u = R::runif(0.0, 1.0);

   while( log(u) > (-0.5*z*z))
   {
      z = R::rexp(rate);
      while(z > (b-a)) z = R::rexp(rate);
      u = R::runif(0.0,1.0);
   }
   return(z+a);
}

// truncn2( mu, sigma, lower, upper)
//
// generates one random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
//======================================================================

double truncn2(double mu, double sigma, double lower, double upper)
{
int change;
 double a, b;
 double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
 double z, tmp, lograt;

 change = 0;
 a = (lower - mu)/sigma;
 b = (upper - mu)/sigma;

 // First scenario
 if( (a == R_NegInf) || (b == R_PosInf))
   {
     if(a == R_NegInf)
       {
     change = 1;
     a = -b;
     b = R_PosInf;
       }

     // The two possibilities for this scenario
     if(a <= 0.45) z = norm_rs(a, b);
     else z = exp_rs(a, b);
     if(change) z = -z;
   }
 // Second scenario
 else if((a * b) <= 0.0)
   {
     // The two possibilities for this scenario
     if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
       {
     z = norm_rs(a, b);
       }
     else z = unif_rs(a,b);
   }
 // Third scenario
 else
   {
     if(b < 0)
       {
     tmp = b; b = -a; a = -tmp; change = 1;
       }

     lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
     if(lograt <= logt2) z = unif_rs(a,b);
     else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
     else z = exp_rs(a,b);
     if(change) z = -z;
   }
   double output;
   output = sigma*z + mu;
 return (output);
}



