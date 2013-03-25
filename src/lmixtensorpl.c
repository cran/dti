//  L-BFGS-B  version
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <stdlib.h>

#include <sys/time.h>
#include <time.h>

int dimx = 0, ngradc = 0;
int ii = 0;
double* siq_init, *bv, *grad;
double alpha, lambda;
extern void F77_NAME(rskmixl2)(double* param, int* npar, double* siq, 
             double* grad, double* bv, int* ng, double* result);

extern void F77_NAME(rskmixl1)(double* param, int* npar, double* siq, 
             double* grad, double* bv, int* ng, double* alpha, double* result);

extern void F77_NAME(rskmixl0)(double* param, int* npar, double* siq, 
             double* grad, double* bv, int* ng, double* lambda, 
             double* alpha, double* result);

extern void F77_NAME(drskml2)(double* param, int* npar, double* siq, 
             double* grad, double* bv, int* ng, double* result);

extern void F77_NAME(drskml1)(double* param, int* npar, double* siq, 
             double* grad, double* bv, int* ng, double* alpha, double* result);

extern void F77_NAME(drskml0)(double* param, int* npar, double* siq, 
             double* grad, double* bv, int* ng, double* lambda, 
             double* alpha, double* result);

extern void paroforient(double *dir, double *angles);
             
typedef struct
{
  int ngrad;
  double* siq;
  double* grad;
  double* bv;
  int fnscale;
} optimex2;

typedef struct
{
  int ngrad;
  double* siq;
  double* grad;
  double* bv;
  double alpha;
  int fnscale;
} optimex1;

typedef struct
{
  int ngrad;
  double* siq;
  double* grad;
  double* bv;
  double lambda;
  double alpha;
  int fnscale;
} optimex0;

typedef struct{
  int order;
  double lambda;
  double alpha;
  double* mix;
  double* orient;
  double* param;
  double value;
} mfunrskml2_ret;

typedef struct{
  int order;
  double lambda;
  double* mix;
  double* orient;
  double* param;
  double value;
} mfunrskml1_ret;

typedef struct{
  int order;
  double* mix;
  double* orient;
  double* param;
  double value;
} mfunrskml0_ret;





double rskmixl2(int param_length, double *param, void* ex){

  double result = 0;

  double* siq = Calloc(ngradc,double);
  int i;
  for(i=0;i<ngradc;i++){
    siq[i] = siq_init[i+ii*ngradc];
  }

  F77_CALL(rskmixl2)(param, &param_length, siq, grad, bv, &ngradc, &result);

  //check for infinity
//  if(result == R_PosInf || result == R_NegInf){
//    return 0;
//  }

  Free(siq);

  return result;
}

double rskmixl1(int param_length, double *param, void* ex){

  double result = 0;

  double* siq = Calloc(ngradc,double);
  int i;
  for(i=0;i<ngradc;i++){
    siq[i] = siq_init[i+ii*ngradc];
  }

  F77_CALL(rskmixl1)(param, &param_length, siq, grad, bv, &ngradc, &alpha, &result);

  //check for infinity
//  if(result == R_PosInf || result == R_NegInf){
//    return 0;
//  }

  Free(siq);

  return result;
}

double rskmixl0(int param_length, double *param, void* ex){

  double result = 0;

  double* siq = Calloc(ngradc,double);
  int i;
  for(i=0;i<ngradc;i++){
    siq[i] = siq_init[i+ii*ngradc];
  }

  F77_CALL(rskmixl0)(param, &param_length, siq, grad, bv, &ngradc, &lambda, &alpha, &result);

  //check for infinity
//  if(result == R_PosInf || result == R_NegInf){
//    return 0;
//  }

  Free(siq);

  return result;
}

void drskml2(int param_length, double* param, double* result, void* ex){

  int i;
  double* siq = Calloc(ngradc, double);
  for(i=0;i<ngradc;i++){
   siq[i] = siq_init[i+ii*ngradc];
  }

  F77_CALL(drskml2)(param, &param_length, siq, grad, bv, &ngradc, result);
  
  Free(siq);
}

void drskml1(int param_length, double* param, double* result, void* ex){

  int i;
  double* siq = Calloc(ngradc, double);
  for(i=0;i<ngradc;i++){
   siq[i] = siq_init[i+ii*ngradc];
  }

  F77_CALL(drskml1)(param, &param_length, siq, grad, bv, &ngradc, &alpha, result);
  
  Free(siq);
}

void drskml0(int param_length, double* param, double* result, void* ex){

  int i;
  double* siq = Calloc(ngradc, double);
  for(i=0;i<ngradc;i++){
   siq[i] = siq_init[i+ii*ngradc];
  }

  F77_CALL(drskml0)(param, &param_length, siq, grad, bv, &ngradc, &lambda, &alpha, result);
  
  Free(siq);
}

mfunrskml2_ret getparam2(int param_length, double* param, double fmin){
  
  mfunrskml2_ret ret_val;
  int i,j;

  int c_ord = (param_length-2)/3;
  
  double sw=1, lambda, alpha;
  double* w_tmp = Calloc(param_length, double);
  double* param_work = Calloc(param_length, double);
  int* o = Calloc(c_ord, int);
  
   for( i = 0; i < param_length; i++){
      param_work[i] = param[i];
   }
   for( i = 0; i < c_ord; i++){
      sw = sw + param[3*i];
      o[i] = i;
   }
//  calculate weights
   double* mix = (double*) R_alloc(c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      mix[i] = param[3*i]/sw;
   }
   
   revsort(mix, o, c_ord);
//  this sorts in decreasing order, now rearrange parameters
//  sorted index in o
   double* orient = (double*) R_alloc(2*c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      w_tmp[i] = param[3*o[i]];
      orient[2*i] = param[3*o[i]+1];
      orient[2*i+1] = param[3*o[i]+2];      
   }
   for(j = 0; j < c_ord; j++){
   while (orient[2*j] < 0) {
      orient[2*j] = orient[2*j] + M_PI;
   }
   while (orient[2*j] > M_PI) {
      orient[2*j] = orient[2*j] - M_PI;
   }
   while (orient[1+2*j] < 0) {
      orient[1+2*j] = orient[1+2*j] + 2 * M_PI;
   }
   while (orient[1+2*j] > 2 * M_PI) {
      orient[1+2*j] = orient[1+2*j] - 2 * M_PI;
   }
  }
   for( i = 0; i < c_ord; i++){
      param[3*i] = w_tmp[i];
      param[3*i+1] = orient[2*i];
      param[3*i+2] = orient[1+2*i];
   }
  lambda = param[3*c_ord];
  alpha = param[3*c_ord+1];
  ret_val.order = c_ord;
  ret_val.lambda = lambda;  
  ret_val.alpha = alpha;  
  ret_val.mix = mix;
  ret_val.orient = orient;  
  ret_val.param = param;
  ret_val.value = fmin;
//  Rprintf("getparam2: npar %i, fmin= %f ", param_length, fmin);
//  Rprintf("\n");                  
 
  Free(o);
  Free(param_work);
  Free(w_tmp);
  
  return ret_val; 
}

mfunrskml1_ret getparam1(int param_length, double* param, double fmin){
  
  mfunrskml1_ret ret_val;
  int i,j;

  int c_ord = (param_length-1)/3;
  
  double sw=1;
  double* w_tmp = Calloc(param_length, double);
  double* param_work = Calloc(param_length, double);
  int* o = Calloc(c_ord, int);
  
   for( i = 0; i < param_length; i++){
      param_work[i] = param[i];
   }
   for( i = 0; i < c_ord; i++){
      sw = sw + param[3*i];
      o[i] = i;
   }
//  calculate weights
   double* mix = (double*) R_alloc(c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      mix[i] = param[3*i]/sw;
   }
   
   revsort(mix, o, c_ord);
//  this sorts in decreasing order, now rearrange parameters
//  sorted index in o
   double* orient = (double*) R_alloc(2*c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      w_tmp[i] = param[3*o[i]];
      orient[2*i] = param[3*o[i]+1];
      orient[2*i+1] = param[3*o[i]+2];      
   }
   for(j = 0; j < c_ord; j++){
   while (orient[2*j] < 0) {
      orient[2*j] = orient[2*j] + M_PI;
   }
   while (orient[2*j] > M_PI) {
      orient[2*j] = orient[2*j] - M_PI;
   }
   while (orient[1+2*j] < 0) {
      orient[1+2*j] = orient[1+2*j] + 2 * M_PI;
   }
   while (orient[1+2*j] > 2 * M_PI) {
      orient[1+2*j] = orient[1+2*j] - 2 * M_PI;
   }
  }
   for( i = 0; i < c_ord; i++){
      param[3*i] = w_tmp[i];
      param[3*i+1] = orient[2*i];
      param[3*i+2] = orient[2*i+1];
   }
  ret_val.lambda = param[3*c_ord];  
  ret_val.order = c_ord;
  ret_val.param = param;
  ret_val.mix = mix;
  ret_val.orient = orient;  
  ret_val.value = fmin;
 
  Free(o);
  Free(param_work);
  Free(w_tmp);
  
  return ret_val; 
}

mfunrskml0_ret getparam0(int param_length, double* param, double fmin){
  
  mfunrskml0_ret ret_val;
  int i,j;

  int c_ord = param_length/3;
  
  double sw=1;
  double* w_tmp = Calloc(param_length, double);
  double* param_work = Calloc(param_length, double);
  int* o = Calloc(c_ord, int);
  
   for( i = 0; i < param_length; i++){
      param_work[i] = param[i];
   }
   for( i = 0; i < c_ord; i++){
      sw = sw + param[3*i];
      o[i] = i;
   }
//  calculate weights
   double* mix = (double*) R_alloc(c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      mix[i] = param[3*i]/sw;
   }
   
   revsort(mix, o, c_ord);
//  this sorts in decreasing order, now rearrange parameters
   double* orient = (double*) R_alloc(2*c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      w_tmp[i] = param[3*o[i]];
      orient[2*i] = param[3*o[i]+1];
      orient[2*i+1] = param[3*o[i]+2];      
   }
   for(j = 0; j < c_ord; j++){
   while (orient[2*j] < 0) {
      orient[2*j] = orient[2*j] + M_PI;
   }
   while (orient[2*j] > M_PI) {
      orient[2*j] = orient[2*j] - M_PI;
   }
   while (orient[1+2*j] < 0) {
      orient[1+2*j] = orient[1+2*j] + 2 * M_PI;
   }
   while (orient[1+2*j] > 2 * M_PI) {
      orient[1+2*j] = orient[1+2*j] - 2 * M_PI;
   }
  }
   for( i = 0; i < c_ord; i++){
      param[3*i] = w_tmp[i];
      param[3*i+1] = orient[2*i];
      param[3*i+2] = orient[2*i+1];
   }
  ret_val.order = c_ord;
  ret_val.mix = mix;
  ret_val.orient = orient;  
  ret_val.param = param;
  ret_val.value = fmin;
 
  Free(o);
  Free(param_work);
  Free(w_tmp);
  
  return ret_val; 
}

void mixtrl2( int* n1, int* siind, int* ngrad0, int* maxcomp, int* maxit, 
          double* grad_in, double* bv_in, double* lambda_in, double* alpha_in,
          double* factr, double* penIC, double* sigma2, double* vert, 
          double* siq_in, double* sigma2_ret, double* orient_ret, int* order_ret,
          double* alpha_ret, double* lambda_ret, double* mix_ret){
// mixtensor: prolate tensors with isotropic compartment, general EV
// optmethod: L-BFGS-B   
// n1 - number of voxel
// siind - initial parameters for orientations
// ngrad0 - number of gradient directions
// maxcomp - maximum number of mixture components
// maxit - maxium number of iterations in L-BFGS-B
// grad_in - array of gradients (3,ngrad)
// bv_in   - vector of b-values
// factr  - control-parameter in  L-BFGS-B (default 1e7)
// penIC   - penalty for model complexity
// sigma2  - vector of  residual variances for best model with initial estimates
// vert    - vertices in polyeder used for initioal orientations
// siq_in  - array of si/s0 (ngrad,n1)
// sigma2_ret - estimated residual variances
// orient_ret - estimated orientations (2,maxcomp,n1)
// order_ret - estimated model order (n1)
// alpha_ret  - estimated alpha values (n1)
// lambda_ret - estimated lambda values (n1)
// mix_ret - estimated mixture coefficients (maxcomp,n1)
  
   int maxcompc = *maxcomp;
   int trace = 0, nREPORT = 1, lmm = 5;
   int fncount = 5, grcount=2;    // number of calls to obj fct in optim
   int ord, iv, i, j, k, l, param_length, param_length_init; 
   double krit, sigma_init, si2new = 0, value = 0;
   double ttt = 0, pgtol = 0;
   double angles[2], dir[3];
   double *param = 0, *param_work = 0, *param_last = 0, 
          *param_tmp = 0, *siq = 0, *lower = 0, *upper = 0;
   char msg[60];
   
   int fail;              // failure code for optim: zero is OK
   double Fmin = 0.;          // minimal value of obj fct in optim 
   
   //Setting global variables
   dimx = *n1;
   siq_init = siq_in; grad = grad_in, ngradc = *ngrad0; bv = bv_in;
   alpha = *alpha_in;
   lambda = *lambda_in;
   
   //Gradient vectors corresponding to minima in spherical coordinates
   
   param_length_init=3*maxcompc+2;
   
   param = (double *) R_alloc(param_length_init, sizeof(double));
   param_work = (double *) R_alloc(param_length_init, sizeof(double));
   param_last = (double *) R_alloc(param_length_init, sizeof(double));
   param_tmp = (double *) R_alloc(param_length_init, sizeof(double));
   lower = (double *) R_alloc(param_length_init, sizeof(double));
   upper = (double *) R_alloc(param_length_init, sizeof(double));
   int *nbd = (int *) R_alloc(param_length_init, sizeof(int));
   siq = (double *) R_alloc(ngradc, sizeof(double));
   
   //initialize param
   for(i=0; i<param_length_init; i++){
       param[i] = 0;
   }
   
   for(ii = 0; ii < dimx; ii++){ 
       for(i=0; i<param_length_init; i++){
          lower[i] = R_NegInf;
          upper[i] = R_PosInf;
	  nbd[i] = 0;
       }
      param_length = param_length_init;
      sigma2_ret[ii] = sigma2[ii];
      
         ord=maxcompc+1;
         
         for (j = 0; j < maxcompc; j++){
            iv = siind[j+maxcompc*ii]-1;
            
            if(iv==-1) iv = j; // this should never happen
            
            dir[0]=vert[iv*3];  //vert dim(max(siind),3)
            dir[1]=vert[1+iv*3];
            dir[2]=vert[2+iv*3];
            
            paroforient(dir, angles);
            orient_ret[0+2*j+2*maxcompc*ii] = angles[0];
            orient_ret[1+2*j+2*maxcompc*ii] = angles[1];
            
            param[3*j] = 1;  
//initial value for w(j) corresponds to volume w(j)^2/(1+sum_k w(j)^2)=1/(1+maxcompc)
            param[3*j+1] = angles[0];
            param[3*j+2] = angles[1];
// set values for lower corresponding to weights to 0
            lower[3*j] = 1e-3;
	    nbd[3*j] = 1;
         }
// lower values for lambda and alpha to get lambda_1>=lambda_2>=0
         lower[3*maxcompc+1] = 0.6;
         upper[3*maxcompc+1] = 10;
	 nbd[3*maxcompc+1] = 2;
         lower[3*maxcompc] = 0.01;
	 nbd[3*maxcompc] = 1;
         // initialize EV-parameter
         param[3*maxcompc] = lambda;//lambda_2
         param[3*maxcompc+1] = alpha;
//alpha; lambda_1=(1+alpha)*lambda; FA=alpha/sqrt((1+alpha)^2+2)
//alpha=1.7 corresponds to a FA of 0.7
         
         sigma_init = sigma2[ii];
         krit = log(sigma_init) + penIC[0];
         
         
         for(l = 0; l < param_length; l++){
            param_work[l] = param[l];
            param_last[l] = param[l];
         }
         
         
         // use AIC/ngrad0, BIC/ngrad0 or AICC/ngrad0 respectively
         for(k = maxcompc; k > 0; k--){
            if(k < ord){
               if(k != maxcompc){
                  for (l = 0; l < param_length-2; l++){
                     param_tmp[l] = param_work[l];
                  }
                  param_tmp[3*k] = param_work[3*k+3];
                  param_tmp[3*k+1] = param_work[3*k+4];
                  param_length=3*k+2;
                  lower[3*k] = 0.01;
                  lower[3*k+1] = 0.6;
                  upper[3*k+1] = 10;
 	          nbd[3*k] = 1;
	          nbd[3*k+1] = 2;
                 
                  for(l= 0; l < param_length; l++){
                     param_work[l] = param_tmp[l];
                     param_last[l] = param_tmp[l];
                  }
               }
               
               optimex2 myoptimpar;
               
               //calculation and estimation of the best fitting model
               
               // initialize ret_val to avoid warnings
               mfunrskml2_ret ret_val;
               ret_val.order = 0;
               ret_val.lambda = 0;
               ret_val.alpha = 0;
               ret_val.mix = NULL;
               ret_val.orient = NULL;
               ret_val.param = NULL;
               ret_val.value = 0;
               
               for(l=0;l<ngradc;l++){
                  siq[l] = siq_init[l+ii*ngradc];
               }
               
               lbfgsb(param_length, lmm, param_work, lower, upper, nbd, 
                      &Fmin, rskmixl2, drskml2, &fail, &myoptimpar, *factr, pgtol,  
                      &fncount, &grcount, *maxit, msg, trace,  nREPORT);
 
               ret_val = getparam2(param_length, param_work, Fmin);
	       
               value = ret_val.value;
               ord = ret_val.order;
               
               if(ord < k)
               {
                  ttt = krit;
                  
                  for(l =0; l < param_length; l++){
                     param_work[l] = param_last[l];
                  }
               }
               else
               {
                  si2new = value/ngradc;
                  
                  ttt = log(si2new) + penIC[ord];
                  
                  for(l =0; l < param_length; l++){
                     param_work[l] = ret_val.param[l];
                  }
               }
               
               
               //use directions corresponding to largest weights as initial directions
               if(ttt < krit){    
                  krit = ttt;
                  
                  order_ret[ii] = ret_val.order;
                  
                  lambda_ret[ii] = ret_val.lambda;
                  alpha_ret[ii] = ret_val.alpha;
                  
                  for(l = 0; l < ord; l++){
                     mix_ret[l+maxcompc*ii] = ret_val.mix[l];
                     orient_ret[0+2*l+2*maxcompc*ii] = ret_val.orient[l*2];
                     orient_ret[1+2*l+2*maxcompc*ii] = ret_val.orient[1+l*2];
                  }
                  for(l = ord; l < maxcompc; l++){
                     mix_ret[l+maxcompc*ii] = 0;
                  }
                  
                  sigma2_ret[ii] = si2new;
               } // end if(ttt < krit)
            } // end if(k < ord)
         } // end model order
      R_CheckUserInterrupt();
   } // end ii
}
void mixtrl1( int* n1, int* siind, int* ngrad0, int* maxcomp, int* maxit, 
          double* grad_in, double* bv_in, double* lambda_in, double* alpha_in,  
          double* factr, double* penIC, double* sigma2, double* vert, 
          double* siq_in, double* sigma2_ret, double* orient_ret, 
          int* order_ret, double* lambda_ret, double* mix_ret){
// mixtensor: prolate tensors with isotropic compartment, fixed FA
// optmethod: L-BFGS-B   
// n1 - number of voxel
// siind - initial parameters for orientations
// ngrad0 - number of gradient directions
// maxcomp - maximum number of mixture components
// maxit - maxium number of iterations in L-BFGS-B
// grad_in - array of gradients (3,ngrad)
// bv_in   - vector of b-values
// alpha_in - (lambda1-lambda_2)/lambda_2
// factr  - control-parameter in  L-BFGS-B (default 1e7)
// penIC   - penalty for model complexity
// sigma2  - vector of  residual variances for best model with initial estimates
// vert    - vertices in polyeder used for initioal orientations
// siq_in  - array of si/s0 (ngrad,n1)
// sigma2_ret - estimated residual variances
// orient_ret - estimated orientations (2,maxcomp,n1)
// order_ret - estimated model order (n1)
// lambda_ret - estimated lambda values (n1)
// mix_ret - estimated mixture coefficients (maxcomp,n1)
   
   int maxcompc = *maxcomp;
   int trace = 0, nREPORT = 1, lmm = 5;
   int fncount = 5, grcount=2;    // number of calls to obj fct in optim
   int ord, iv, i, j, k, l, param_length, param_length_init;
   double krit, sigma_init, si2new = 0, value = 0;
   double ttt = 0, pgtol = 0;
   double angles[2], dir[3];
   double *param = 0, *param_work = 0, *param_last = 0, 
          *param_tmp = 0, *siq = 0, *lower = 0, *upper = 0;
   char msg[60];
   
   int fail;              // failure code for optim: zero is OK
   double Fmin = 0.;          // minimal value of obj fct in optim 
   
   //Setting global variables
   dimx = *n1;
   siq_init = siq_in; grad = grad_in, ngradc = *ngrad0; bv = bv_in;
   alpha = *alpha_in;
   lambda = *lambda_in;
   
   
   //Gradient vectors corresponding to minima in spherical coordinates
   
   param_length_init=3*maxcompc+1;
   
   param = (double *) R_alloc(param_length_init, sizeof(double));
   param_work = (double *) R_alloc(param_length_init, sizeof(double));
   param_last = (double *) R_alloc(param_length_init, sizeof(double));
   param_tmp = (double *) R_alloc(param_length_init, sizeof(double));
   lower = (double *) R_alloc(param_length_init, sizeof(double));
   upper = (double *) R_alloc(param_length_init, sizeof(double));
   int *nbd = (int *) R_alloc(param_length_init, sizeof(int));
   siq = (double *) R_alloc(ngradc, sizeof(double));
   
   //initialize param
   for(i=0; i<param_length_init; i++){
       param[i] = 0;
   }
   
   for(ii = 0; ii < dimx; ii++){ 
      for(i=0; i<param_length_init; i++){
         lower[i] = R_NegInf;
         upper[i] = R_PosInf;
	 nbd[i] = 0;
      }
      param_length = param_length_init;
      sigma2_ret[ii] = sigma2[ii];
      
         ord=maxcompc+1;
         
         for (j = 0; j < maxcompc; j++){
            iv = siind[j+maxcompc*ii]-1;
            
            if(iv==-1) iv = j; // this should never happen
            
            dir[0]=vert[iv*3];  //vert dim(max(siind),3)
            dir[1]=vert[1+iv*3];
            dir[2]=vert[2+iv*3];
            
            paroforient(dir, angles);
            orient_ret[0+2*j+2*maxcompc*ii] = angles[0];
            orient_ret[1+2*j+2*maxcompc*ii] = angles[1];
            
            param[3*j] = 1;  
//initial value for w(j) corresponds to volume w(j)^2/(1+sum_k w(j)^2)=1/(1+maxcompc)
            param[3*j+1] = angles[0];
            param[3*j+2] = angles[1];
// set values for lower corresponding to weights to 0
            lower[3*j] = 1e-3;
	    nbd[3*j] = 1;
         }
         lower[3*maxcompc] = 0.01;
	 nbd[3*maxcompc] = 1;
// lower value for alpha to get lambda_1>=lambda_2
         // initialize EV-parameter
         param[3*maxcompc] = lambda;//lambda_2
//alpha; lambda_1=(1+alpha)*lambda; FA=alpha/sqrt((1+alpha)^2+2)
//alpha=1.7 corresponds to a FA of 0.7
         
         sigma_init = sigma2[ii];
         krit = log(sigma_init) + penIC[0];
         
         
         for(l = 0; l < param_length; l++){
            param_work[l] = param[l];
            param_last[l] = param[l];
         }
         
         
         // use AIC/ngrad0, BIC/ngrad0 or AICC/ngrad0 respectively
         for(k = maxcompc; k > 0; k--){
            if(k < ord){
               if(k != maxcompc){
                  for (l = 0; l < param_length-1; l++){
                     param_tmp[l] = param_work[l];
                  }
                  param_tmp[3*k] = param_work[3*k+3];
                  param_length=3*k+1;
                  lower[3*k] = 0.01;
 	          nbd[3*k] = 1;
                  
                  for(l= 0; l < param_length; l++){
                     param_work[l] = param_tmp[l];
                     param_last[l] = param_tmp[l];
                  }
               }
               
               optimex1 myoptimpar;
               
               //calculation and estimation of the best fitting model
               
               // initialize ret_val to avoid warnings
               mfunrskml1_ret ret_val;
               ret_val.order = 0;
               ret_val.lambda = 0;
               ret_val.mix = NULL;
               ret_val.orient = NULL;
               ret_val.param = NULL;
               ret_val.value = 0;
               
               for(l=0;l<ngradc;l++){
                  siq[l] = siq_init[l+ii*ngradc];
               }
               
               lbfgsb(param_length, lmm, param_work, lower, upper, nbd, 
                      &Fmin, rskmixl1, drskml1, &fail, &myoptimpar, *factr, pgtol,  
                      &fncount, &grcount, *maxit, msg, trace,  nREPORT);

               ret_val = getparam1(param_length, param_work, Fmin);
	       
               value = ret_val.value;
               ord = ret_val.order;
               
               if(ord < k)
               {
                  ttt = krit;
                  
                  for(l =0; l < param_length; l++){
                     param_work[l] = param_last[l];
                  }
               }
               else
               {
                  si2new = value/ngradc;
                  
                  ttt = log(si2new) + penIC[ord];
                  
                  for(l =0; l < param_length; l++){
                     param_work[l] = ret_val.param[l];
                  }
               }
               
               
               //use directions corresponding to largest weights as initial directions
               if(ttt < krit){    
                  krit = ttt;
                  
                  order_ret[ii] = ret_val.order;
                  
                  lambda_ret[ii] = ret_val.lambda;
                  
                  for(l = 0; l < ord; l++){
                     mix_ret[l+maxcompc*ii] = ret_val.mix[l];
                     orient_ret[0+2*l+2*maxcompc*ii] = ret_val.orient[l*2];
                     orient_ret[1+2*l+2*maxcompc*ii] = ret_val.orient[1+l*2];
                  }
                  for(l = ord; l < maxcompc; l++){
                     mix_ret[l+maxcompc*ii] = 0;
                  }
                  
                  sigma2_ret[ii] = si2new;
               } // end if(ttt < krit)
            } // end if(k < ord)
         } // end model order
      R_CheckUserInterrupt();
   } // end ii
}
void mixtrl0( int* n1, int* siind, int* ngrad0, int* maxcomp, int* maxit, 
          double* grad_in, double* bv_in, double* lambda_in, double* alpha_in,  
          double* factr, double* penIC, double* sigma2, double* vert, 
          double* siq_in, double* sigma2_ret, double* orient_ret, 
          int* order_ret, double* mix_ret){
// mixtensor: prolate tensors with isotropic compartment, fixed eigenvalues
// optmethod: L-BFGS-B   
// n1 - number of voxel
// siind - initial parameters for orientations
// ngrad0 - number of gradient directions
// maxcomp - maximum number of mixture components
// maxit - maxium number of iterations in L-BFGS-B
// grad_in - array of gradients (3,ngrad)
// bv_in   - vector of b-values
// lambda_in - lambda_2
// alpha_in - (lambda1-lambda_2)/lambda_2
// factr  - control-parameter in  L-BFGS-B (default 1e7)
// penIC   - penalty for model complexity
// sigma2  - vector of  residual variances for best model with initial estimates
// vert    - vertices in polyeder used for initioal orientations
// siq_in  - array of si/s0 (ngrad,n1)
// sigma2_ret - estimated residual variances
// orient_ret - estimated orientations (2,maxcomp,n1)
// order_ret - estimated model order (n1)
// mix_ret - estimated mixture coefficients (maxcomp,n1)
   
   int maxcompc = *maxcomp;
   int trace = 0, nREPORT = 1, lmm = 5;
   int fncount = 5, grcount=2;    // number of calls to obj fct in optim
   int ord, iv, i, j, k, l, param_length, param_length_init; 
   double krit, sigma_init, si2new = 0, value = 0;
   double ttt = 0, pgtol = 0;
   double angles[2], dir[3];
   double *param = 0, *param_work = 0, *param_last = 0, 
          *param_tmp = 0, *siq = 0, *lower = 0, *upper = 0;
   char msg[60];
   
   int fail;              // failure code for optim: zero is OK
   double Fmin = 0.;          // minimal value of obj fct in optim 
   
   //Setting global variables
   alpha = *alpha_in; lambda = *lambda_in;;
   dimx = *n1;
   siq_init = siq_in; grad = grad_in, ngradc = *ngrad0; bv = bv_in;
   
   
   //Gradient vectors corresponding to minima in spherical coordinates
   
   param_length_init=3*maxcompc;
   
   param = (double *) R_alloc(param_length_init, sizeof(double));
   param_work = (double *) R_alloc(param_length_init, sizeof(double));
   param_last = (double *) R_alloc(param_length_init, sizeof(double));
   param_tmp = (double *) R_alloc(param_length_init, sizeof(double));
   lower = (double *) R_alloc(param_length_init, sizeof(double));
   upper = (double *) R_alloc(param_length_init, sizeof(double));
   int *nbd = (int *) R_alloc(param_length_init, sizeof(int));
   siq = (double *) R_alloc(ngradc, sizeof(double));
   
   //initialize param
   for(i=0; i<param_length_init; i++){
       param[i] = 0;
   }
   
   for(ii = 0; ii < dimx; ii++){ 
       for(i=0; i<param_length_init; i++){
          lower[i] = R_NegInf;
          upper[i] = R_PosInf;
	  nbd[i] = 0;
       }
      param_length = param_length_init;
      sigma2_ret[ii] = sigma2[ii];
      
         ord=maxcompc+1;
         
         for (j = 0; j < maxcompc; j++){
            iv = siind[j+maxcompc*ii]-1;
            
            if(iv==-1) iv = j; // this should never happen
            
            dir[0]=vert[iv*3];  //vert dim(max(siind),3)
            dir[1]=vert[1+iv*3];
            dir[2]=vert[2+iv*3];
            
            paroforient(dir, angles);
            orient_ret[0+2*j+2*maxcompc*ii] = angles[0];
            orient_ret[1+2*j+2*maxcompc*ii] = angles[1];
            
            param[3*j] = 1;  
//initial value for w(j) corresponds to volume w(j)^2/(1+sum_k w(j)^2)=1/(1+maxcompc)
            param[3*j+1] = angles[0];
            param[3*j+2] = angles[1];
// set values for lower corresponding to weights to 0
            lower[3*j] = 1e-3;
	    nbd[3*j] = 1;
         }
// lower value for alpha to get lambda_1>=lambda_2
         // initialize EV-parameter
         
         sigma_init = sigma2[ii];
         krit = log(sigma_init) + penIC[0];
         
//          Rprintf("param from orient:\n");
         
         for(l = 0; l < param_length; l++){
            param_work[l] = param[l];
            param_last[l] = param[l];
//             Rprintf(" %f ", param[l]);
         }
         
//          Rprintf("\n");
         
         // use AIC/ngrad0, BIC/ngrad0 or AICC/ngrad0 respectively
         for(k = maxcompc; k > 0; k--){
            if(k < ord){
               if(k != maxcompc){
                  for (l = 0; l < param_length-2; l++){
                     param_tmp[l] = param_work[l];
                  }
                  param_length=3*k;
                  
                  for(l= 0; l < param_length; l++){
                     param_work[l] = param_tmp[l];
                     param_last[l] = param_tmp[l];
                  }
               }
               
               
               //neuer check wenn elses drin sind
               
               optimex0 myoptimpar;
               
               //calculation and estimation of the best fitting model
               
               // initialize ret_val to avoid warnings
               mfunrskml0_ret ret_val;
               ret_val.order = 0;
               ret_val.mix = NULL;
               ret_val.orient = NULL;
               ret_val.param = NULL;
               ret_val.value = 0;
               
               for(l=0;l<ngradc;l++){
                  siq[l] = siq_init[l+ii*ngradc];
               }
               
               lbfgsb(param_length, lmm, param_work, lower, upper, nbd, 
                      &Fmin, rskmixl0, drskml0, &fail, &myoptimpar, *factr, pgtol,  
                      &fncount, &grcount, *maxit, msg, trace,  nREPORT);

               ret_val = getparam0(param_length, param_work, Fmin);
	       
               value = ret_val.value;
               ord = ret_val.order;
               
               if(ord < k)
               {
                  ttt = krit;
                  
                  for(l =0; l < param_length; l++){
                     param_work[l] = param_last[l];
                  }
               }
               else
               {
                  si2new = value/ngradc;
                  
                  ttt = log(si2new) + penIC[ord];
                  
                  for(l =0; l < param_length; l++){
                     param_work[l] = ret_val.param[l];
                  }
               }
               
               
               //use directions corresponding to largest weights as initial directions
               if(ttt < krit){    
                  krit = ttt;
                  
                  order_ret[ii] = ret_val.order;
                                    
                  for(l = 0; l < ord; l++){
                     mix_ret[l+maxcompc*ii] = ret_val.mix[l];
                     orient_ret[0+2*l+2*maxcompc*ii] = ret_val.orient[l*2];
                     orient_ret[1+2*l+2*maxcompc*ii] = ret_val.orient[1+l*2];
                  }
                  for(l = ord; l < maxcompc; l++){
                     mix_ret[l+maxcompc*ii] = 0;
                  }
                  
                  sigma2_ret[ii] = si2new;
               } // end if(ttt < krit)
            } // end if(k < ord)
         } // end model order
      R_CheckUserInterrupt();
   } // end ii
}

