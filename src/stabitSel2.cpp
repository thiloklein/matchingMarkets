#define ARMA_NO_DEBUG
// Uncomment the above line if you want to disable all run-time checks.
// This will result in faster code, but you first need to make sure that your code runs correctly!

// [[Rcpp::depends(RcppArmadillo,RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>

using namespace Rcpp;
arma::colvec mvrnormArma(arma::colvec mu, arma::mat sigma, int ncols);
double truncn(double bound, bool lb, double mu, double sigma);

// [[Rcpp::export]]
List stabitSel2(Rcpp::List Xr, Rcpp::List Rr, Rcpp::List Wr, 
  arma::colvec One, arma::colvec Two, int T,
  Rcpp::List offOutr, Rcpp::List offSelr,
  arma::mat sigmabarbetainverse, arma::mat sigmabaralphainverse,
  int niter, double n, arma::colvec l, Rcpp::List Pr, arma::colvec p,
  bool binary, bool selection, int censored, bool ntu, bool gPrior,
  bool display_progress=true) {

  // Market identifiers.
  int Two1=Two(0), TwoN=Two(1), One1=One(0), OneN=One(1);

  // ---------------------------------------------
  // Independent variables
  // ---------------------------------------------
    
  arma::field<arma::mat> X(T); // observed groups
  for(int i=0; i<T; i++){
    X(i) = Rcpp::as<arma::mat>(Xr[i]);  
  }    
  arma::field<arma::mat> W(TwoN); // all feasible groups
  for(int i=Two1; i<TwoN; i++){
    W(i) = Rcpp::as<arma::mat>(Wr[i]);  
  }  
  const int kX = X(0).n_cols; // number of parameters
  const int kW = W(0).n_cols; // number of parameters
  
  // ---------------------------------------------
  // Outcome variables.
  // ---------------------------------------------
  
  arma::field<arma::colvec> R(T), offOut(T), Y(T);
  for(int i=0; i<T; i++){
    R(i)      = Rcpp::as<arma::colvec>(Rr[i]);      // dependent variable
    offOut(i) = Rcpp::as<arma::colvec>(offOutr[i]); // offset
    if(binary == TRUE){
      Y(i)    = arma::zeros(R(i).n_rows,1); // latent outcome (set to zero)
    } else{
      Y(i)    = -offOut(i) + R(i);          // R = offOut + Xb + e
    }
  }
  
  arma::field<arma::colvec> V(TwoN), offSel(TwoN);
  for(int i=Two1; i<TwoN; i++){
    V(i) = arma::zeros(W(i).n_rows,1);              // latent match valuation
    offSel(i) = Rcpp::as<arma::colvec>(offSelr[i]); // offset
  }
  
  // ---------------------------------------------
  // Index of other group composed of residual borrowers in two-group markets.
  // ---------------------------------------------
  
  arma::field<arma::colvec> P(TwoN); // all feasible groups
  for(int i=Two1; i<TwoN; i++){
    P(i) = Rcpp::as<arma::colvec>(Pr[i]);  
  } 
  
  // ---------------------------------------------
  // Priors.
  // ---------------------------------------------
  
  // alpha
  arma::colvec alphabar = arma::zeros(kW,1);
  if(gPrior == FALSE){
    arma::mat sigmabaralpha = 10*arma::eye(kW,kW);
    arma::mat sigmabaralphainverse = arma::inv(sigmabaralpha);
  }
  arma::mat alphabaroversigma = sigmabaralphainverse*alphabar;  

  // beta
  arma::colvec betabar = arma::zeros(kX,1);
  if(gPrior == FALSE){
    arma::mat sigmabarbeta = 10*arma::eye(kX,kX);
    arma::mat sigmabarbetainverse = arma::inv(sigmabarbeta);    
  }
  arma::colvec betabaroversigma = sigmabarbetainverse*betabar;
    
  // delta
  double deltabar = 0;
  double sigmabarsquaredelta = 10;
  
  // ---------------------------------------------
  // Prior modes for parameters.
  // ---------------------------------------------

  // alpha
  arma::colvec alpha = alphabar;
  arma::colvec alphahat = arma::zeros(kW,1);
  arma::mat sigmahatalpha = arma::zeros(kW,kW);  

  // beta
  arma::colvec beta = betabar;
  arma::colvec betahat = arma::zeros(kX,1);
  arma::mat sigmahatbeta = arma::zeros(kX,kX);
  
  // delta
  double delta = deltabar;
  double deltahat = 0;
  double sigmahatsquaredelta = 0;
  
  // xi
  double a = 2, b = 1, ahat, bhat;
  double sigmasquarexiinverse = b*(a-1);
  double sigmasquarexi = 1/sigmasquarexiinverse;
  
  // sums
  arma::mat sum1    = arma::zeros(kX,kX); 
  arma::colvec sum2 = arma::zeros(kX,1);
  arma::mat sum3    = arma::zeros(kW,kW); 
  arma::colvec sum4 = arma::zeros(kW,1);
  arma::colvec sum5 = arma::zeros(1,1); 
  arma::colvec sum6 = arma::zeros(1,1);
  arma::colvec sum7 = arma::zeros(1,1);
  
  // eta
  arma::colvec eta(2*TwoN);
    
  // ---------------------------------------------
  // Matrices for parameter draws.
  // ---------------------------------------------
  
  arma::mat alphadraws(kW,niter);
  arma::mat betadraws(kX,niter);
  arma::mat deltadraws(1,niter);
  arma::mat etadraws(2*TwoN,niter);
  arma::mat sigmasquarexidraws(1,niter);
  
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
      for(int t=One1; t < OneN; t++){ // one-group-market identifiers.        
        Yhat = arma::as_scalar( X(t)*beta );
        if(R(t)(0) == 1){
          Y(t) = truncn(-offOut(t)(0), TRUE, Yhat, 1.0); // bound, lb, mu, sigma
        } else{
          Y(t) = truncn(-offOut(t)(0), FALSE, Yhat, 1.0); // bound, lb, mu, sigma
        }
      }     
      if(selection == TRUE){
        for(int t=Two1; t < TwoN; t++){ // two-group-market identifiers.        
          for(int G=0; G < 2; G++){
            Yhat = arma::as_scalar( X(t).row(G)*beta + (V(t)(G)-W(t).row(G)*alpha)*delta );
            if(R(t)(G) == 1){
              Y(t)(G) = truncn(-offOut(t)(G), TRUE, Yhat, 1.0); // bound, lb, mu, sigma
            } else{
              Y(t)(G) = truncn(-offOut(t)(G), FALSE, Yhat, 1.0); // bound, lb, mu, sigma
            }            
          }
        }                    
      } else{
        for(int t=Two1; t < TwoN; t++){ // two-group-market identifiers.        
          for(int G=0; G < 2; G++){
            Yhat = arma::as_scalar( X(t).row(G)*beta );
            if(R(t)(G) == 1){
              Y(t)(G) = truncn(-offOut(t)(G), TRUE, Yhat, 1.0); // bound, lb, mu, sigma
            } else{
              Y(t)(G) = truncn(-offOut(t)(G), FALSE, Yhat, 1.0); // bound, lb, mu, sigma
            } 
          }
        }            
      }
    }
        
    // ---------------------------------------------
    // Simulate the latent valuations conditional on the latent group outcomes, data, and parameters.
    // ---------------------------------------------
    
    double Vhat, Vupperbar, Vlowerbar, eqsum, neqmaxsum;
    
    if(selection == TRUE){ 

      if(ntu==TRUE){ // non-transferable utility
        for(int t=Two1; t < TwoN; t++){
          // considering _upper_bounds_ for all unobserved groups in market t.
          Vupperbar = arma::max( V(t).rows(0,1) ); // max valuation of observed groups as upper bound
          for(int G=2; G < l(t); G++){
            Vhat    = arma::as_scalar( W(t).row(G)*alpha );
            V(t)(G) = truncn(Vupperbar, FALSE, Vhat, 1.0); // bound, lb, mu, sigma
          }
        
          // considering _lower_bounds_ for observed groups in market t.
          for(int G=0; G < 2; G++){ // observed groups A and B.
            Vlowerbar = arma::max( V(t).rows(2,l(t)-1) ); // max valuation of unobserved groups as lower bound
            Vhat = arma::as_scalar( W(t).row(G)*alpha + ((Y(t)(G) - X(t).row(G)*beta)*delta)/(sigmasquarexi + pow(delta,2)) );
            if(V(t)(1-G) > Vlowerbar){ // no bounds if valuation of other observed group exceeds that of all unobserved groups.
              V(t)(G) = ::Rf_rnorm(Vhat, sqrt(sigmasquarexi/(sigmasquarexi+pow(delta,2))) );
            } else{
              V(t)(G) = truncn(Vlowerbar, TRUE, Vhat, sqrt(sigmasquarexi/(sigmasquarexi+pow(delta,2))) ); // bound, lb, mu, sigma
            }
          }
        }
      } else{ // transferable utility
        for(int t=Two1; t < TwoN; t++){
          // considering _upper_bounds_ for all unobserved groups in market t.
          eqsum = 0; // reset to zero
          for(int G=0; G<2; G++){
            eqsum += V(t)(G) + offSel(t)(G);
          }
          // first draws for "ego" groups...
          for(int G=2; G<p(t); G++){ // p(t)-1 gives index of last "ego" group
            Vupperbar = eqsum - ( V(t)(P(t)(G)-1) + offSel(t)(P(t)(G)-1) ) - offSel(t)(G);
            Vhat = arma::as_scalar( W(t).row(G)*alpha );
            V(t)(G) = truncn(Vupperbar, FALSE, Vhat, 1.0); // bound, lb, mu, sigma
          }
          // ...then for "other" groups because of interdependence between their valuations
          for(int G=p(t); G<l(t); G++){ // p(t) gives index of first "other" group
            Vupperbar = eqsum - ( V(t)(P(t)(G)-1) + offSel(t)(P(t)(G)-1) ) - offSel(t)(G);
            Vhat = arma::as_scalar( W(t).row(G)*alpha );
            V(t)(G) = truncn(Vupperbar, FALSE, Vhat, 1.0); // bound, lb, mu, sigma
          }
          
          // considering _lower_bounds_ for observed groups in market t.
          neqmaxsum = arma::max( V(t).rows(2,p(t)-1) + V(t).rows(p(t),l(t)-1) + offSel(t).rows(2,p(t)-1) + offSel(t).rows(p(t),l(t)-1) );
          for(int G=0; G<2; G++){
            Vlowerbar = neqmaxsum - ( V(t)(1-G) + offSel(t)(1-G) ) - offSel(t)(G);            
            Vhat = arma::as_scalar( W(t).row(G)*alpha + ( (Y(t)(G) - X(t).row(G)*beta )*delta ) / (sigmasquarexi + pow(delta,2)) );
            V(t)(G) = truncn(Vlowerbar, TRUE, Vhat, sqrt(sigmasquarexi/(sigmasquarexi+pow(delta,2))) ); // bound, lb, mu, sigma
          }  
        }
      } 
    }

    // ---------------------------------------------
    // Simulate each group of parameters conditional on latents, data, and all other parameters.
    // ---------------------------------------------
  
    // ---------------------------------------------
    // beta.
    // ---------------------------------------------
    
    sum1.zeros(kX,kX); // reset to zero
    sum2.zeros(kX,1);  // reset to zero
    for(int t = One1; t < OneN; t++){
      sum1 += trans(X(t))*X(t);
      sum2 += trans(X(t))*Y(t);
    }
    if(selection == TRUE){
      for(int t = Two1; t < TwoN; t++){
        sum1 += trans(X(t))*X(t);
        sum2 += trans(X(t))*( Y(t) - delta*(V(t).rows(0,1) - W(t).rows(0,1)*alpha) );
      }      
    } else{
      for(int t = Two1; t < TwoN; t++){
        sum1 += trans(X(t))*X(t);
        sum2 += trans(X(t))*Y(t);
      }
    }
    sigmahatbeta = arma::inv(sigmabarbetainverse + (1/sigmasquarexi)*sum1);
    betahat = -sigmahatbeta * (-betabaroversigma - (1/sigmasquarexi)*sum2);
    beta = mvrnormArma(betahat,sigmahatbeta,kX); 
    
    if(selection == TRUE){
      
      // ---------------------------------------------
      // alpha.
      // ---------------------------------------------
      
      sum3.zeros(kW,kW); // reset to zero
      sum4.zeros(kW,1);  // reset to zero
      for(int t=Two1; t<TwoN; t++){ // two-group-market identifiers.
        sum3 += trans(W(t))*W(t) + (pow(delta,2)/sigmasquarexi)*trans(W(t).rows(0,1))*W(t).rows(0,1);
        sum4 += -trans(W(t))*V(t) + (delta/sigmasquarexi)*(trans(W(t).rows(0,1))*(Y(t) - X(t)*beta - V(t).rows(0,1)*delta));
      } 
      sigmahatalpha = arma::inv(sigmabaralphainverse + sum3); 
      alphahat = -sigmahatalpha * (-alphabaroversigma + sum4); 
      alpha = mvrnormArma(alphahat,sigmahatalpha,kW);
      
      // ---------------------------------------------
      // delta.
      // ---------------------------------------------
      
      sum5.zeros(1,1); // reset to zero
      sum6.zeros(1,1); // reset to zero
      for(int t=Two1; t<TwoN; t++){ // two-group-market identifiers.
        for(int G=0; G < 2; G++){
          sum5 += pow(V(t)(G) - W(t).row(G)*alpha, 2);
          sum6 += (Y(t)(G) - X(t).row(G)*beta) * (V(t).row(G) - W(t).row(G)*alpha);
        }
      }      
      sigmahatsquaredelta = 1/(1/sigmabarsquaredelta + (1/sigmasquarexi)*arma::as_scalar(sum5));
      deltahat = -sigmahatsquaredelta*(-deltabar/sigmabarsquaredelta - (1/sigmasquarexi)*arma::as_scalar(sum6));      
      if(censored == 1){ // from below, positive covariation of residuals
        delta = truncn(0, TRUE, deltahat, sqrt(sigmahatsquaredelta) ); // bound, lb, mu, sigma
      } else if(censored == 2){ // from above, negative covariation of residuals
        delta = truncn(0, FALSE, deltahat, sqrt(sigmahatsquaredelta) ); // bound, lb, mu, sigma
      } else{ // not censored
        delta = ::Rf_rnorm(deltahat, sqrt(sigmahatsquaredelta));
      }
      
      // ---------------------------------------------
      // eta.
      // ---------------------------------------------
      for(int t=Two1; t<TwoN; t++){
        eta.rows(2*t,2*t+1) = V(t).rows(0,1) - W(t).rows(0,1)*alpha;
      }
    }
    
    // ---------------------------------------------
    // sigma xi
    // ---------------------------------------------
    
    if(binary == FALSE){
      sum7.zeros(1,1); // reset to zero
      for(int t = One1; t < OneN; t++){
        for(int G=0; G<2; G++){
          sum7 += pow(Y(t)(G) - X(t).row(G)*beta, 2);
        }
      }
      if(selection == TRUE){
        for(int t = Two1; t < TwoN; t++){
          for(int G=0; G<2; G++){
            sum7 += pow(Y(t)(G) - X(t).row(G)*beta - delta*(V(t)(G) - W(t).row(G)*alpha), 2);
          }
        }          
      } else{
        for(int t = Two1; t < TwoN; t++){
          for(int G=0; G<2; G++){
            sum7 += pow(Y(t)(G) - X(t).row(G)*beta, 2);
          }
        }  
      }
      ahat = a + n/2;          // shape
      bhat = 1/(1/b + arma::as_scalar(sum7)/2); // scale = 1/rate
      sigmasquarexiinverse = ::Rf_rgamma(ahat, bhat); // shape, scale
      sigmasquarexi = 1/sigmasquarexiinverse;
    }
    
    // ---------------------------------------------  
    // save this iterations draws.
    // ---------------------------------------------
    
    betadraws.col(iter) = beta;
    if(selection == TRUE){
      alphadraws.col(iter) = alpha;
      deltadraws.col(iter) = delta;
      etadraws.col(iter) = eta;
    }
    sigmasquarexidraws.col(iter) = sigmasquarexi;
  }
  
  // ---------------------------------------------  
  // Return the parameter draws.
  // ---------------------------------------------  

  if((binary==TRUE) & (selection==TRUE)){
    return List::create(  
      // parameter draws
      Named("betadraws") = alphadraws,
      Named("alphadraws") = betadraws,
      Named("deltadraws") = deltadraws,
      Named("etadraws") = etadraws
    );  
  } else if((binary==TRUE) & (selection==FALSE)){
    return List::create(  
      Named("alphadraws") = betadraws
    );  
  } else if((binary==FALSE) & (selection==TRUE)){
        return List::create(  
      // parameter draws
      Named("betadraws") = alphadraws,
      Named("alphadraws") = betadraws,
      Named("deltadraws") = deltadraws,
      Named("etadraws") = etadraws,
      Named("sigmasquarexidraws") = sigmasquarexidraws
    );
  } else if((binary==FALSE) & (selection==FALSE)){
    return List::create(  
      // parameter draws
      Named("alphadraws") = betadraws,
      Named("sigmasquarexidraws") = sigmasquarexidraws
    );      
  } else{
      return 0;
  }
}

// ---------------------------------------------  
// random multivariate normal sample generator using RcppArmadillo
// from http://gallery.rcpp.org/articles/simulate-multivariate-normal/
// ---------------------------------------------  

arma::colvec mvrnormArma(arma::colvec mu, arma::mat sigma, int ncols) {
  arma::rowvec y = as<arma::rowvec>(rnorm(ncols)); //arma::randn(1,ncols)
  return arma::trans( arma::trans(mu) + y*arma::chol(sigma));
}

// ---------------------------------------------  
// Returns one draw from truncated normal distribution (mu,sigma^2) with range
// (bound,+inf) if lb=TRUE and 
// (-inf,bound) if lb=FALSE
// from http://athens.src.uchicago.edu/jenni/econ319_2003/lecture.html
// ---------------------------------------------  

double truncn(double bound, bool lb, double mu, double sigma){
  
  double c, z, w;
  
  // 1. standardised cut-off c for truncation from below or above 
  if(lb == TRUE){
    c = (bound-mu)/sigma;
  } else{
    c = -(bound-mu)/sigma;
  }
  
  // 2. standardised draw using Geweke's (1991)
  if(c < 0.45){ // normal rejection sampling
    z = ::Rf_rnorm(0.0,1.0);
    while(z < c){
      z = ::Rf_rnorm(0.0,1.0);
    } 
  } else{ // exponential rejection sampling
    z = -log(1-::Rf_runif(0.0,1.0))/c;
    w = ::Rf_runif(0.0,1.0);
    while(w > exp(-0.5*pow(z,2))){
      z = -log(1-::Rf_runif(0.0,1.0))/c;
      w = ::Rf_runif(0.0,1.0);
    }
    z = z+c;
  }

  // 3. reverse standardisation
  if(lb == TRUE){
    return mu + sigma*z;
  } else{
    return mu - sigma*z;
  }
}
