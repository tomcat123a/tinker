#include <RcppArmadillo.h>
#include<ctime>
           // [[Rcpp::depends(RcppArmadillo)]]
           using namespace std;
/*		   
		   // [[Rcpp::export]]
		  
          extern "C" SEXP fa(SEXP ys, SEXP Xs) {
          
          Rcpp::NumericVector yr(ys);                 // creates Rcpp vector from SEXP
          Rcpp::NumericMatrix Xr(Xs);                 // creates Rcpp matrix from SEXP
          int n = Xr.nrow(), k = Xr.ncol();
          
          arma::mat X(Xr.begin(), n, k, false);       // reuses memory and avoids extra copy
          arma::colvec y(yr.begin(), yr.size(), false);
          
          arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
          arma::colvec resid = y - X*coef;            // residuals
          
          double sig2 = arma::as_scalar( arma::trans(resid)*resid/(n-k) );
          // std.error of estimate
          arma::colvec stderrest = arma::sqrt( sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
          
          return Rcpp::List::create(
            Rcpp::Named("coefficients") = coef,
            Rcpp::Named("stderr")       = stderrest,
            Rcpp::Named("1")       = 1
          ) ;
		  }
*/		  
inline	double S(double x,double y){
       if(y > abs(x)) {
	      return (0);
       }else{  
          if(x>0){
	          return (x-y);
          }else{
              return (x+y);
          }
       }	
	}	


 
		  
		  
		  
		// [[Rcpp::export]]


arma::colvec PALM(SEXP w0_r,SEXP X1_r,SEXP X2_r,SEXP t1_r,double lambda1,double lambda2=0,double alpha=1,double te=0.0001,
               int maxtimes=0,bool palmquiet=true){
// clock_t cl1 = clock();
 //clock_t cl2;
 //clock_t cl3 = 0;
  Rcpp::NumericVector wr(w0_r);                 // creates Rcpp vector from SEXP
  Rcpp::NumericMatrix X1r(X1_r);                 // creates Rcpp matrix from SEXP
  Rcpp::NumericMatrix X2r(X2_r);                 // creates Rcpp matrix from SEXP
  Rcpp::NumericVector t1r(t1_r);                 // creates Rcpp vector from SEXP
  const arma::mat X1(X1r.begin(), X1r.nrow(), X1r.ncol(), false);  
  const arma::mat X2(X2r.begin(), X2r.nrow(), X2r.ncol(), false);  
  const arma::colvec t1(t1r.begin(), t1r.size(), false);  
//  cl2 = clock() - cl1;
  arma::colvec w = wr;  
  long int p = w.size();
  const arma::mat X1_t1=X1.each_col()%t1;
     int count = -1;
  int ind = p + 1;     
  int prevind = p;     
  double wdif = 0;
  double l1err = 0;
  double parderi;
  double wtemp;
  double objv;
  double objvprev;
  arma::colvec tempv1_t1=X1_t1 * w;
  arma::colvec f1=1 - 1 / ( 1 + exp(-tempv1_t1));
  if(R_IsNA(X2(0,0))==false){
	arma::rowvec a = 0.999 / (sum( X1%X1 ) / 4 + lambda2 * sum( X2%X2 ) / 4 + lambda1 * (1 - alpha) );
	//std::cout<<"a(0)==" << std::endl << a(0) << std::endl;
    arma::colvec tempv2=X2 * w;
	arma::colvec f2=tempv2 % ( 1 / ( 2 + exp(-tempv2) + exp(tempv2) ) );
    objvprev =  - sum(log( 1 - f1 )) + lambda2 * sum(1 / (1 + exp(-tempv2)) % log(1 + exp(-tempv2) ) + 1 / (1 + exp(tempv2)) % log(1 + exp(tempv2) )       )
	            + lambda1 * ( alpha * sum( abs(w) ) + (1 - alpha) / 2 * sum(w % w) ) ;
			   
	do{
      if(ind==p+1){
        count++;
        ind=1;
        l1err=0;
		
      }
      if(wdif!=0){
        l1err=l1err+abs(wdif);
        tempv1_t1=tempv1_t1+X1_t1.col(prevind-1) * wdif;
        tempv2=tempv2+X2.col(prevind-1) * wdif;
        //cl2=clock();
        f1=1 - 1 / ( 1 + exp(-tempv1_t1));
        f2=tempv2 % ( 1 / ( 2 + exp(-tempv2) + exp(tempv2) ) );
		parderi=-sum(X1_t1.col(ind-1) % f1) - lambda2 * sum(X2.col(ind-1) % f2) + lambda1 * (1-alpha) * w(ind-1);
		//f2=tempv2 % ( 1 / ( 4 + tempv2%tempv2 + tempv2%tempv2%tempv2%tempv2/12 ) );
		//cl2=clock()-cl2;cl3=cl3+cl2;
      }else{
         parderi=-sum(X1_t1.col(ind-1) % f1) - lambda2 * sum(X2.col(ind-1) % f2) + lambda1 * (1-alpha) * w(ind-1);
      }
      wtemp=w(ind-1);
	  if(ind==1){
          w(ind-1)= wtemp-a(ind-1) * parderi;
	  }else{
		  w(ind-1)= S( wtemp - a(ind-1) * parderi, a(ind-1) * lambda1);
	  }
      wdif=w(ind-1)-wtemp;
      prevind=ind;
      ind++;
	  if(ind == p+1 ){
		 if( maxtimes==0){
			 objv=( - sum(log( 1 - f1 )) + lambda2 * sum(1 / (1 + exp(-tempv2)) % log(1 + exp(-tempv2) ) + 1 / (1 + exp(tempv2)) % log(1 + exp(tempv2) )       )
	            + lambda1 * ( alpha * sum( abs(w) )+ (1 - alpha) / 2 * sum(w % w) ) );
            if(  abs( objv - objvprev) < te * objvprev){
		     //  cout<<"objv" << objv;
		       break;
	        }else{
		       objvprev=objv;
	        }
	     }else{
	           if((maxtimes > 0 && count == maxtimes)){
		           break;
	          } 
	     }
	  }
	  //std::cout << ind -1 << " wdif= " << wdif << " wtemp - a(ind-1)*pa= " << wtemp - a(ind-2) * parderi << " a*lbd1= " <<a(ind-2) * lambda1 << " S= " << w(ind-2)<< std::endl;
    }while( true );
  }else{
	  arma::rowvec a = 0.999 / (sum( X1%X1 ) / 4 + lambda1 * (1 - alpha) );
      objvprev = ( - sum(log( 1 - f1 )) + lambda1 * ( alpha * sum( abs(w) ) + (1 - alpha) / 2 * sum(w % w) ) );
	  do{
      if(ind == p+1){
        count++;
        ind=1;
        l1err=0;
      }
      if(wdif!=0){
        l1err=l1err+abs(wdif);
        tempv1_t1=tempv1_t1 + X1_t1.col(prevind-1) * wdif;
        parderi=-sum( X1_t1.col(ind - 1) % f1) + lambda1 * (1 - alpha) * w(ind-1);
        f1=1 - 1 / ( 1 + exp(-tempv1_t1));
      }else{
        parderi=-sum( X1_t1.col(ind - 1) % f1) + lambda1 * (1 - alpha) * w(ind-1);
      }
      wtemp=w(ind-1);
      if(ind==1){
          w(ind-1)= wtemp-a(ind-1) * parderi;
	  }else{
		  w(ind-1)= S( wtemp - a(ind-1) * parderi, a(ind-1) * lambda1);
	  }
      wdif=w(ind-1)-wtemp;
      prevind=ind;
      ind++;
	  if(ind == p+1 ){
		 if( maxtimes==0){
			 objv=( - sum(log( 1 - f1 )) + lambda1 * ( alpha * sum( abs(w) )+ (1 - alpha) / 2 * sum(w % w) ) );
            if(  abs( objv - objvprev) < te * objvprev){
		       
		       break;
	        }else{
		       objvprev=objv;
	        }
	     }else{
	           if((maxtimes > 0 && count == maxtimes)){
		           break;
	          } 
	     }
	  }
    }while( true  );
  }
if(!palmquiet){
    cout<<"objv" << objv << endl;
   std::cout << "number of iterations: " << std::endl << count << std::endl;
   std::cout << "number of effective features detected:" << sum(sign(abs(w))) << std::endl;
}
//std::cout << "l1err== " << l1err << std::endl;
//cl2=clock()-cl1;
//cout << cl3/(double)CLOCKS_PER_SEC<< ":"<<cl2/(double)CLOCKS_PER_SEC <<endl;;
  return (w);
}


		  
		  
		// [[Rcpp::export]]


arma::colvec iPALM(SEXP w0_r,SEXP X1_r,SEXP X2_r,SEXP t1_r,double lambda1,double lambda2=0,double alpha=1,double te=0.0001,
               int maxtimes=500,double alpha_1 = 0.2, double beta_1 = 0.8, bool palmquiet=true){
 //clock_t cl1 = clock();
 //clock_t cl2;
 //clock_t cl3 = 0;
  Rcpp::NumericVector wr(w0_r);                 // creates Rcpp vector from SEXP
  Rcpp::NumericMatrix X1r(X1_r);                 // creates Rcpp matrix from SEXP
  Rcpp::NumericMatrix X2r(X2_r);                 // creates Rcpp matrix from SEXP
  Rcpp::NumericVector t1r(t1_r);                 // creates Rcpp vector from SEXP
  const arma::mat X1(X1r.begin(), X1r.nrow(), X1r.ncol(), false);  
  const arma::mat X2(X2r.begin(), X2r.nrow(), X2r.ncol(), false);  
  const arma::colvec t1(t1r.begin(), t1r.size(), false);  
 //  cl2 = clock() - cl1;
  arma::colvec w = wr;  
  arma::colvec wprev = w;
  
  long int p = w.size();
  const arma::mat X1_t1=X1.each_col()%t1;
  int count = -1;      
  int ind = p + 1;     
  int prevind = p;     
  double wdif = 0;
  double l1err = 0;
  double parderi;
  double wtemp;
  double wprevtemp;
  
  double alpha1 = alpha_1;
  double beta1 = beta_1;
  double wy;
  
  
  arma::colvec tempv1_t1=X1_t1 * w;
  arma::colvec f1=1 - 1 / ( 1 + exp(-tempv1_t1));
  if(R_IsNA(X2(0,0))==false){
	arma::rowvec a = 0.999 / (sum( X1%X1 ) / 4 + lambda2 * sum( X2%X2 ) / 4 + lambda1 * (1 - alpha) );
	//std::cout<<"a(0)==" << std::endl << a(0) << std::endl;
    arma::colvec tempv2=X2 * w;
	arma::colvec z_tempv2 = tempv2;
	arma::colvec z_tempv1_t1 = tempv1_t1;
	arma::colvec f2=tempv2 % ( 1 / ( 2 + exp(-tempv2) + exp(tempv2) ) );
    do{
      if(ind==p+1){
        count++;
        ind=1;
        l1err=0;
		
      }
      if(wdif!=0){
        l1err=l1err+abs(wdif);
        tempv1_t1=tempv1_t1+X1_t1.col(prevind-1) * wdif;
		z_tempv1_t1 = tempv1_t1 + beta1 * (w(ind-1)  - wprev(ind - 1)) * X1_t1.col(ind-1);
        tempv2=tempv2+X2.col(prevind-1) * wdif;
		z_tempv2=tempv2 + beta1 * (w(ind-1)  - wprev(ind - 1)) * X2.col(ind-1);
        //cl2=clock();
        f1=1 - 1 / ( 1 + exp(-z_tempv1_t1));
        f2=z_tempv2 % ( 1 / ( 2 + exp(-z_tempv2) + exp(z_tempv2) ) );
		parderi=-sum(X1_t1.col(ind-1) % f1) - lambda2 * sum(X2.col(ind-1) % f2) + lambda1 * (1-alpha) * w(ind-1);
		//f2=tempv2 % ( 1 / ( 4 + tempv2%tempv2 + tempv2%tempv2%tempv2%tempv2/12 ) );
		//cl2=clock()-cl2;cl3=cl3+cl2;
      }else{
         parderi=-sum(X1_t1.col(ind-1) % f1) - lambda2 * sum(X2.col(ind-1) % f2) + lambda1 * (1-alpha) * w(ind-1);
      }
      wtemp=w(ind-1);
	  wprevtemp=wprev(ind-1);
	  wprev(ind-1)=w(ind-1);
	  wy = wtemp + alpha1 * (wtemp - wprevtemp);
	  
	  if(ind==1){
		  
          w(ind-1)= wy-(1 - 2 * alpha1) / (1 + 2 * beta1) * a(ind-1) * parderi;
	  }else{
		  w(ind-1)= S( wy - (1 - 2 * alpha1) / (1 + 2 * beta1) * a(ind-1) * parderi, a(ind-1) * lambda1);
	  }
      wdif=w(ind-1)-wtemp;
      prevind=ind;
      ind++;
	  //std::cout << ind -1 << " wdif= " << wdif << " wtemp - a(ind-1)*pa= " << wtemp - a(ind-2) * parderi << " a*lbd1= " <<a(ind-2) * lambda1 << " S= " << w(ind-2)<< std::endl;
    }while( ! ( (ind == p+1 && count >= 1 && l1err < te&& 
            maxtimes == 0)||(maxtimes > 0 && count == maxtimes) ) );
  }
if(!palmquiet){
   std::cout << "number of iterations: " << std::endl << count << std::endl;
   std::cout << "number of effective features detected:" << sum(sign(abs(w))) << std::endl;
}
//std::cout << "l1err== " << l1err << std::endl;
//cl2=clock()-cl1;
//cout << cl3/(double)CLOCKS_PER_SEC<< ":"<<cl2/(double)CLOCKS_PER_SEC <<endl;;
  return (w);
}