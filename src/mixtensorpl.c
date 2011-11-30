#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <stdlib.h>

#include <sys/time.h>
#include <omp.h>
#include <time.h>

int dim_x = 0, dim_y = 0, dim_z = 0, ngrad0c = 0;
int i1 = 0, i2 = 0, i3 = 0;
double* siq_init, *penc, *grad, pen = 100;

extern void F77_NAME(mfunpl0)(double* param, double* siq, double* grad, int* m, 
		       int* lpar, int* ngrad, double* pen, double* tmp,
		       double* tmp2, double* result);

extern void F77_NAME(mfunpl0g)(double* param, double* siq, double* grad,
				int* m, int* lpar, int* ngrad,
				double* tmp1,     //(ngrad*m),# z(n,m)
				double* tmp2,//(m*m),# v(m,m)
				double* tmp3,//(ngrad),# w(n) need w((m+1):n) for solver dgelsy
				double* tmp4,//(ngrad*m),# dkgj(n,m)
				double* tmp5,//(ngrad*m),# dkgj2(n,m)
				double* tmp6,//(ngrad*m),# ddkdphig(n,m)
				double* tmp7,//(ngrad*m),# ddkdetag(n,m)
				double* tmp8,//(m*m),# dvdth(m,m)
				double* tmp9,//(m*m*m),# dvdphi(m,m,m)
				double* tmp10,//(m*m*m),# dvdeta(m,m,m)
				double* tmp11,//(ngrad*m*3),# dzdpars(n,m,3)
				double* tmp12,//(m*lpar),# dwdpars(m,lpar)
				double* tmp13,//(m*lpar),# dwdpars2(m,lpar)
				double* tmp14,//(ngrad*m),# zs(n,m)
				double* tmp15,//(ngrad*m),# work1(n,m)
				double* tmp16,//(ngrad*m),# work2(n,m)
				double* tmp17,//(ngrad),# scopy(n)
				double* pen, 
				double* result //result dfdpar double(lpar)
			);

extern void F77_NAME(mfunpl0h)(double* param, double* siq, double* grad, int* m,
			       int* lpar, int* ngrad, double* tmp1, 
			       double* tmp2, double* tmp3, double* tmp4,
			       double* result);

extern void F77_NAME(mfunpl0w)(double* param, double* w, double* siq, 
			       double* grad, int* m, int* param_length, 
			       int* ngrad, double* tmp1, double* result);

extern void F77_NAME(mfunpl0h)(double* param, double* siq, double* grad,
			       int* m, int *lpar, int* ngrad, double* tmp1,
			       double* w, double* tmp2, double* tmp3,
			       double* result);

typedef struct
{
  int ngrad;
  double* siq;
  double* grad;
  double* pen;
  int fnscale;
} optimex;

typedef struct
{
  int ngrad;
  double *siq;
  double *grad;
  double *w;
  double *z;
  double *qv;
  double *dqv;
  double *fv;
  double *dfv;
  double *work1;
} optimexpl;

typedef struct{
  int order;
  double lev_par;
  double lev_sum_w;
  double* mix;
  double* orient;
  double* param;
  double value;
  double* w;
} mfunplwghts0_ret;

/* C++ code
public static T[] SubArrayDeepClone<T>(this T[] data, int index, int length){
  // ggf. fuer schoenes array kopieren (http://stackoverflow.com/questions/943635/c-arrays-getting-a-sub-array-from-an-existing-array)
    T[] arrCopy = new T[length];
    Array.Copy(data, index, arrCopy, 0, length);
    using (MemoryStream ms = new MemoryStream())
    {
        var bf = new BinaryFormatter();
        bf.Serialize(ms, arrCopy);
        ms.Position = 0;
        return (T[])bf.Deserialize(ms);
    }
    }*/

static R_INLINE int get_ind2d_img(int x, int y){
  return x + y*dim_x;
}

static R_INLINE int get_ind3d_img(int x, int y, int z){
  return x + y*dim_x + z*dim_x*dim_y;
}

static R_INLINE int get_ind4d_img(int x, int y, int z, int i){
  return x + y*dim_x + z*dim_x*dim_y + i*dim_x*dim_y*dim_z;
}

static R_INLINE int get_ind5d_img(int x, int y, int z, int i, int k, int dim){
  return x + y*dim_x + z*dim_x*dim_y + i*dim_x*dim_y*dim_z + k*dim_x*dim_y*dim_z*dim;
}

static R_INLINE int get_ind2d(int x, int y, int dim){
  return x + y*dim;
}


int compare_doubles (const void *a, const void *b){
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  
  return (*da < *db) - (*da > *db);
}


int contains(double* array, int length, double value){
  int i;
  for(i = 0; i < length; i++){
    if(value == array[i]){
      return 1;
    }
  }
  return 0;
}

int contains_int(int* array, int length, int value){
  int i;
  for(i = 0; i < length; i++){
    if(value == array[i]){
      return 1;
    }
  }
  return 0;
}


double* transpose(double* matrix, int dim1, int dim2){
  int i, j;
  double* t_matrix = Calloc(dim1*dim2, double);
  for(i = 0; i < dim1; i++){
    for(j = 0; j< dim2; j++){
      t_matrix[get_ind2d(i,j,dim1)] = matrix[get_ind2d(j,i,dim2)];
    }
  }
  return t_matrix;
}
  

void paroforient(double *dir, double *angles){
  angles[0] = acos(dir[2]);        // angles[0]=theta; angles[1]=phi
  double sth = sin(angles[0]);    
  angles[1] = 0;                  

  if (sth < 1e-8){
    angles[0] = 0;
  }else{
    double z = dir[0]/sth;
    if(abs(z)>=1){
      angles[1] = z < 0 ? 0 : M_PI;
    }else{
      angles[1] = acos(z) * sign(dir[1]);
    }
    if(angles[1]<0) angles[1]=angles[1]+2*M_PI;
  }

  return;
}


double mfunpl0(int param_length, double *param, void* ex){

  int m = (param_length-1)/2;
  double* z = Calloc(ngrad0c*m, double);
  double* w = Calloc(ngrad0c, double);
  double result = 0;

  double* siq = Calloc(ngrad0c,double);
  int i;
  for(i=0;i<ngrad0c;i++){
    siq[i] = siq_init[get_ind4d_img(i1,i2,i3,i)];
  }

  F77_CALL(mfunpl0)(param, siq, grad, &m, &param_length, &ngrad0c, 
		    &pen, z, w, &result);


  //check for infinity
  if(result == R_PosInf || result == R_NegInf){
    return 0;
  }

  Free(z);
  Free(w);
  Free(siq);
  return result;
}


double mfunpl0h(int param_length, double* param, void* ex){

  int m = (param_length-1)/2;
  double result = 0;
  int i;

  double* siq = Calloc(ngrad0c, double);
  for(i=0;i<ngrad0c;i++){
    siq[i] = siq_init[get_ind4d_img(i1,i2,i3,i)];
  }
  
  double* z =  Calloc(ngrad0c*m, double);
  double* w =  Calloc(ngrad0c, double);
  double* b =  Calloc(ngrad0c, double);
  double* work1 = Calloc(ngrad0c, double);
  
  F77_CALL(mfunpl0h)(param, siq, grad, &m, &param_length, &ngrad0c, 
		     z, w, b, work1, &result);

  Free(siq);
  Free(z);
  Free(w);
  Free(b);
  Free(work1);

  return result;
}


void gmfunpl0(int param_length, double* param, double* result, void* ex){

  int m = (param_length-1)/2;

  double* z = Calloc(ngrad0c*m, double);
  double* v = Calloc(m*m, double);
  double* w = Calloc(ngrad0c, double);
  double* dkgj = Calloc(ngrad0c*m, double);
  double* dkgj2 = Calloc(ngrad0c*m, double);
  double* ddkdphig = Calloc(ngrad0c*m, double);
  double* ddkgetag = Calloc(ngrad0c*m, double);
  double* dvdth = Calloc(m*m, double);
  double* dvdphi = Calloc(m*m*m, double);
  double* dvdeta = Calloc(m*m*m, double);
  double* dzdpars = Calloc(ngrad0c*m*3, double);
  double* dwdpars = Calloc(m*param_length, double);
  double* dwdpars2 = Calloc(m*param_length, double);
  double* zs = Calloc(ngrad0c*m, double);
  double* work1 = Calloc(ngrad0c*m, double);
  double* work2 = Calloc(ngrad0c*m, double);
  double* scopy = Calloc(ngrad0c, double);

  int i;
  double* siq = Calloc(ngrad0c, double);
  for(i=0;i<ngrad0c;i++){
	siq[i] = siq_init[get_ind4d_img(i1,i2,i3,i)];
  }

  F77_CALL(mfunpl0g)(param, siq, grad, &m, &param_length, &ngrad0c,
		     z, v, w, dkgj, dkgj2, ddkdphig, ddkgetag,
		     dvdth, dvdphi, dvdeta, dzdpars, dwdpars, dwdpars2, zs, work1, work2, scopy,
		     &pen, result);
  Free(siq);
  Free(z);
  Free(v);
  Free(w);
  Free(dkgj);
  Free(dkgj2);
  Free(ddkdphig);
  Free(ddkgetag);
  Free(dvdth);
  Free(dvdphi);
  Free(dvdeta);
  Free(dzdpars);
  Free(dwdpars);
  Free(dwdpars2);
  Free(zs);
  Free(work1);
  Free(work2);
  Free(scopy);
}




/**
 *  Fall mixtensor, BFGS/CG
 *
 */
mfunplwghts0_ret mfunplwghts0(int param_length, double* param, double* siq){
  
  mfunplwghts0_ret ret_val;
  int i,j;

  int m = (param_length-1)/2, flag = 0, ord = 0;

  if (param[0] < 0){
    param[0] = 0;
  }
  
  double pen = 100;
  double tmp2 = 0;
  double* param_work = Calloc(param_length, double);
  double* tmp = Calloc(ngrad0c*m, double);
  double* w = Calloc(ngrad0c, double);
  double* w_red = Calloc(m, double);
  double* w_red2 = Calloc(m, double);
  
  //copy param
   for( i = 0; i < param_length; i++){
      param_work[i] = param[i];
   }

  //initialize w
  for(i=0;i<ngrad0c;i++){
      w[i] = 0;
  }
  
  F77_CALL(mfunpl0)(param_work, siq, grad, &m, &param_length, &ngrad0c, 
		    &pen, tmp, w, &tmp2);
			
  tmp2 = 0;
  Free(tmp);
  tmp = Calloc(ngrad0c*m, double);

  for(i = 0; i < m; i++){
    if(w[i] > 0){
      w_red[i] = w[i];
    }else{
      w_red[i] = 0;
    }
  }  

  for(i = 0; i < m; i++){
    if(w_red[i] == R_PosInf || w_red[i] == R_NegInf){
      flag = 1;
    }
  }
  
  
  if(flag == 0){
    F77_CALL(mfunpl0w)(param_work, w_red, siq, grad, &m, &param_length, &ngrad0c, tmp, &tmp2);
  }else{
    tmp2 = pow(10,20);
  }


  for(i = 0; i < m; i++){
	w_red2[i] = w_red[i];
  }
  
  qsort(w_red, m, sizeof(double), compare_doubles);
  int* o = Calloc(m, int);
  for(i = 0; i<m; i++){
	o[i] = -1;
  }

  for(i=0; i<m; i++){
    for(j=0; j<m; j++){
      if(w_red2[j] == w_red[i] && contains_int(o,m,j) == 0){
		o[i] = j;
		j = m;
      }
    }
  }
  
  double sum_w = 0;
  for (i = 0; i < m; i++){
    if(w_red[i] > 0){
      ord++;
      sum_w += w_red[i];	
    }
  }


  /** correct order:
    * if order == 0, C for-loops wouldn't start
    * and arrays wouldn't be allocated
    * R solves this with o <- o[1:ord] = o[1:0] = o[1]
    * so c_ord will be 1 iff ord == 0, else c_ord = ord
    **/
  int c_ord = ord;
  if (c_ord == 0) c_ord = 1;

  // here ord will not be replaced, because zero-check is important
  double* mix = (double*) R_alloc(ord, sizeof(double));
  
  if(ord > 0){
	for(i = 0; i < ord; i++){
		mix[i] = w_red[i]/sum_w;
	}
  }else{
	mix = NULL;
  }
  
  double* or = (double*) R_alloc(2*c_ord, sizeof(double));
  for(i = 0; i < c_ord; i++){
	or[get_ind2d(0,i,2)] = param[2*o[i]+1];
	or[get_ind2d(1,i,2)] = param[2*o[i]+2];
  }
  
  for(j = 0; j < c_ord; j++){
	while (or[get_ind2d(0,j,2)] < 0) {
		or[get_ind2d(0,j,2)] = or[get_ind2d(0,j,2)] + M_PI;
	}
	while (or[get_ind2d(0,j,2)] > M_PI) {
		or[get_ind2d(0,j,2)] = or[get_ind2d(0,j,2)] - M_PI;
	}
	while (or[get_ind2d(1,j,2)] < 0) {
		or[get_ind2d(1,j,2)] = or[get_ind2d(1,j,2)] + 2 * M_PI;
	}
	while (or[get_ind2d(1,j,2)] > 2 * M_PI) {
		or[get_ind2d(1,j,2)] = or[get_ind2d(1,j,2)] - 2 * M_PI;
	}
  }
  
  for(i = 0; i < c_ord; i++){       //par[0] stays
	param[2*i+1] = or[2*i];
	param[2*i+2] = or[2*i+1];
	
  }
  
  ret_val.lev_par = param[0];  //lev[1]
  ret_val.lev_sum_w = -log(sum_w);   //lev[2]
  ret_val.order = ord;
  ret_val.param = param;
  ret_val.value = tmp2;
  ret_val.orient = or;  
  ret_val.mix = mix;
 
  Free(tmp);
  Free(param_work);
  Free(w);
  Free(w_red);
  Free(w_red2);
  Free(o);
  
  return ret_val; 
}


mfunplwghts0_ret mfunplwghts0h(int param_length, double* param, double* siq){
  
  if(param[0]<0){
    param[0]=0;
  }

  int i,j;
  
  mfunplwghts0_ret ret_val;
  int m = (param_length-1)/2;
  int ord = 0;

  double* param_work = Calloc(param_length, double);
  double* z = Calloc(ngrad0c*m, double);
  double* w = Calloc(ngrad0c, double);
  double* w_red = Calloc(m, double);
  double* w_red2 = Calloc(m, double);
  double* b = Calloc(ngrad0c, double);
  double* work1 = Calloc(ngrad0c, double);
  double result = 0;

  F77_CALL(mfunpl0h)(param, siq, grad, &m, &param_length, &ngrad0c, z,
		     w, b, work1, &result);

  result *= result;

    //copy param
  for( i = 0; i < param_length; i++){
     param_work[i] = param[i];
  }
  
  for(i = 0; i < m; i++){
    w_red[i] = w[i];
    w_red2[i] = w[i];
  }
 
  qsort(w_red, m, sizeof(double), compare_doubles);
  int* o = Calloc(m, int);
  for(i = 0; i<m; i++){
	o[i] = -1;
  }

  for(i=0; i<m; i++){
    for(j=0; j<m; j++){
      if(w_red2[j] == w_red[i] && contains_int(o,m,j) == 0){
		o[i] = j;
		j = m;
      }
    }
  }
  
  double sum_w = 0;
  for (i = 0; i < m; i++){
    if(w_red[i] > 0){
      ord++;
      sum_w += w_red[i];	
    }
  }


  /** correct order:
    * if order == 0, C for-loops wouldn't start
    * and arrays wouldn't be allocated
    * R solves this with o <- o[1:ord] = o[1:0] = o[1]
    * so c_ord will be 1 iff ord == 0, else c_ord = ord
    **/
  int c_ord = ord;
  if (c_ord == 0) c_ord = 1;

  // here ord will not be replaced, because zero-check is important
  double* mix = (double*) R_alloc(ord, sizeof(double));
  
  if(ord > 0){
	for(i = 0; i < ord; i++){
		mix[i] = w_red[i]/sum_w;
	}
  }else{
	mix = NULL;
  }
  
  double* or = (double*) R_alloc(2*c_ord, sizeof(double));
  for(i = 0; i < c_ord; i++){
	or[get_ind2d(0,i,2)] = param[2*o[i]+1];
	or[get_ind2d(1,i,2)] = param[2*o[i]+2];
  }

  for(j = 0; j < c_ord; j++){
	while (or[get_ind2d(0,j,2)] < 0) {
		or[get_ind2d(0,j,2)] = or[get_ind2d(0,j,2)] + M_PI;
	}
	while (or[get_ind2d(0,j,2)] > M_PI) {
		or[get_ind2d(0,j,2)] = or[get_ind2d(0,j,2)] - M_PI;
	}
	while (or[get_ind2d(1,j,2)] < 0) {
		or[get_ind2d(1,j,2)] = or[get_ind2d(1,j,2)] + 2 * M_PI;
	}
	while (or[get_ind2d(1,j,2)] > 2 * M_PI) {
		or[get_ind2d(1,j,2)] = or[get_ind2d(1,j,2)] - 2 * M_PI;
	}
  }
  
  for(i = 0; i < c_ord; i++){       //par[0] stays
	param[2*i+1] = or[2*i];
	param[2*i+2] = or[2*i+1];
	
  }  
  
  ret_val.lev_par = param[0];  //lev[1]
  ret_val.lev_sum_w = -log(sum_w);   //lev[2]
  ret_val.order = ord;;
  ret_val.param = param;
  ret_val.value = result;
  ret_val.orient = or;  
  ret_val.mix = mix;
  
  return ret_val; 
}

mfunplwghts0_ret mfunwghts(int param_length, double* param, double siq){
  
  int m = (param_length-1)/3;

  mfunplwghts0_ret ret_val;
  int i, j, ord = 0;
  int lpar = 2*m+1;
  int w_length = param_length - lpar;
  double*  w = (double*)Calloc(w_length,double);
  double* w2 = (double*)Calloc(w_length,double);

  for(i=lpar-1; i<param_length; i++){
	 w[i-(lpar-1)] = param[i];
	w2[i-(lpar-1)] = param[i];
  }


  qsort(w, w_length, sizeof(double), compare_doubles);
  int* o = Calloc(w_length, int);
  for(i = 0; i<w_length; i++){
	o[i] = -1;
  }

  for(i=0; i<w_length; i++){
	for(j=0; j<w_length; j++){
		if(w2[j] == w[i] && contains_int(o,w_length,j) == 0){
			o[i] = j;
			j = w_length;
		}
	}
  }
  
  double sum_w = 0;
  for (i = 0; i < w_length; i++){
	if(w[i] > 0){
		ord++;
		sum_w += w[i];
	}
  }

  /** correct order:
    * if order == 0, C for-loops wouldn't start
    * and arrays wouldn't be allocated
    * R solves this with o <- o[1:ord] = o[1:0] = o[1]
    * so c_ord will be 1 iff ord == 0, else c_ord = ord
    **/
  int c_ord = ord;
  if (c_ord == 0) c_ord = 1;

  Free(w2);

  w2 = (double*) R_alloc(c_ord, sizeof(double));
  for(i = 0; i < c_ord; i++)
	w2[i] = w[i];

  // here ord will not be replaced, because zero-check is important
  double* mix = (double*) R_alloc(ord, sizeof(double));
  if(ord > 0){
	for(i = 0; i < ord; i++){
		mix[i] = w[i]/sum_w;
	}
  }else{
	mix = NULL;
  }

  double* or = (double*) R_alloc(2*c_ord, sizeof(double));
  for(i = 0; i < c_ord; i++){
	or[get_ind2d(0,i,2)] = param[2*o[i]+1];
	or[get_ind2d(1,i,2)] = param[2*o[i]+2];
  }
  
  for(j = 0; j < c_ord; j++){
	while (or[get_ind2d(0,j,2)] < 0) {
		or[get_ind2d(0,j,2)] = or[get_ind2d(0,j,2)] + M_PI;
	}
	while (or[get_ind2d(0,j,2)] > M_PI) {
		or[get_ind2d(0,j,2)] = or[get_ind2d(0,j,2)] - M_PI;
	}
	while (or[get_ind2d(1,j,2)] < 0) {
		or[get_ind2d(1,j,2)] = or[get_ind2d(1,j,2)] + 2 * M_PI;
	}
	while (or[get_ind2d(1,j,2)] > 2 * M_PI) {
		or[get_ind2d(1,j,2)] = or[get_ind2d(1,j,2)] - 2 * M_PI;
	}
  }

  
  
  for(i = 0; i < c_ord; i++){       //par[0] stays
	param[2*i+1] = or[2*i];
	param[2*i+2] = or[2*i+1];
  }

  ret_val.lev_par = param[0];  //lev[1]
  ret_val.lev_sum_w = -log(sum_w);   //lev[2]
  ret_val.order = ord;
  ret_val.param = param;
  ret_val.orient = or;  
  ret_val.mix = mix;
  ret_val.value = 0.;
  ret_val.w = w2;
 
  Free(o);
  Free(w);

  return ret_val; 
}

mfunplwghts0_ret mfunwghtsi(int param_length, double* param, double siq){
  
  int m = (param_length-2)/3;
  double w0 = 1;
  mfunplwghts0_ret ret_val;

  if (m > 0) {
	int i, j, ord = 0;
	int lpar = 2*m+1;
	int w_length = param_length - lpar -1;
	double*  w = (double*)Calloc(w_length,double);
	double* w2 = (double*)Calloc(w_length,double);
	
	for(i=lpar; i<param_length; i++){
		 w[i-lpar] = param[i];
		w2[i-lpar] = param[i];
	}
	
	w0 = param[lpar-1];
	
	qsort(w, w_length, sizeof(double), compare_doubles);
	int* o = Calloc(w_length, int);
	for(i = 0; i<w_length; i++){
		o[i] = -1;
	}
	
	for(i=0; i<w_length; i++){
		for(j=0; j<w_length; j++){
			if(w2[j] == w[i] && contains_int(o,w_length,j) == 0){
				o[i] = j;
				j = w_length;
			}
		}
	}
	
	double sum_w = 0.;
	for (i = 0; i < w_length; i++){
		if(w[i] > 0){
			ord++;
			sum_w += w[i];	
		}
	}
	
	sum_w += fmax2(w0,0.);
	
	/** correct order:
	* if order == 0, C for-loops wouldn't start
	* and arrays wouldn't be allocated
	* R solves this with o <- o[1:ord] = o[1:0] = o[1]
	* so c_ord will be 1 iff ord == 0, else c_ord = ord
	**/
	int c_ord = ord;
	if (c_ord == 0) c_ord = 1;
	
	// here ord will not be replaced, because zero-check is important
	double* mix = (double*) R_alloc(ord, sizeof(double));
	if(ord > 0){
		for(i = 0; i < ord; i++){
			mix[i] = w[i]/sum_w;
		}
	}else{
		mix = NULL;
	}
	
	double* or = (double*) R_alloc(2*c_ord, sizeof(double));
	for(i = 0; i < c_ord; i++){
		or[get_ind2d(0,i,2)] = param[2*o[i]+1];
		or[get_ind2d(1,i,2)] = param[2*o[i]+2];
	}
	
	for(j = 0; j < m; j++){
	
		while (or[get_ind2d(0,j,2)] < 0) {
			or[get_ind2d(0,j,2)] = or[get_ind2d(0,j,2)] + M_PI;
		}
		while (or[get_ind2d(0,j,2)] > M_PI) {
			or[get_ind2d(0,j,2)] = or[get_ind2d(0,j,2)] - M_PI;
		}
		while (or[get_ind2d(1,j,2)] < 0) {
			or[get_ind2d(1,j,2)] = or[get_ind2d(1,j,2)] + 2 * M_PI;
		}
		while (or[get_ind2d(1,j,2)] > 2 * M_PI) {
			or[get_ind2d(1,j,2)] = or[get_ind2d(1,j,2)] - 2 * M_PI;
		}
	}
	
	for(i = 0; i < c_ord; i++){       //par[0] stays
		param[2*i+1] = or[2*i];
		param[2*i+2] = or[2*i+1];
	}
	
	w2 = (double*) R_alloc(c_ord+1, sizeof(double));

	for(i = 0; i < c_ord; i++) w2[i+1] = w[i];

	w2[0] = w0;

	ret_val.lev_sum_w = -log(sum_w);   //lev[2]
	ret_val.order = ord;
	ret_val.param = param;
	ret_val.orient = or;  
	ret_val.mix = mix;
	ret_val.w = w2;
	
	Free(o);
	Free(w);
	
  } else {
	ret_val.lev_sum_w = 0;   //lev[2]
	ret_val.order = 0 ;
	ret_val.param = &param[0];
	ret_val.orient = NULL;
	ret_val.mix = NULL;
	ret_val.w = &w0;
  }

  ret_val.lev_par = param[0];  //lev[1]
  ret_val.value = 0.;

  return ret_val; 
}

void mixture2( int* method, int* optmethod, int* n1, int* n2, int* n3, 
	       int* mask, int* siind, int* ngrad, int* ngrad0, int* maxcomp,
	       int* maxit, 
	       double* pen, double* grad_in, double* reltol,
	       double* th, double* penIC, double* sigma2, double* vert, 
	       double* orient, double* siq_in,
	       double* sigma2_ret, double* orient_ret, 
	       int* order_ret, double* lev_ret, double* mix_ret){
	
	
	// method=1 == mixtensor, method=2 == mixtensoriso
	// optmethod:
	//   1 : BFGS
	//   2 : CG
	//   3 : Nelder-Mead
	//   4 : L-BFGS-B
	//
	// i1,i2,i3 and n1,n2,n3: dim_x, dim_y, dim_z current/total respectively
	// ngrad : # gradients
	// maxcomp : maximal order of mix ???
	
	int maxcompc = *maxcomp, methodc = *method, optmethodc = *optmethod;
	int trace = 0, nREPORT = 1;
	int fncount = 5, grcount=2;    // number of calls to obj fct in optim
	int type = 1; // for CG, takes 1, 2 or 3 see R help(optim)
	int ord, iv, i, j, k, l, param_length, param_length_init, tmp_int; 
	double krit, sigma_init, si2new = 0, value = 0;
	double ttt = 0, abstol = R_NegInf, intol = 1;
	double alpha = 1, beta = 0.5, gamma = 2; //Optim parameters for Nelder Mead
	double angles[2], dir[3];
	double *param = 0, *param_red = 0, *param_ret = 0, *param_work = 0, *param_last = 0;
	
	int fail;              // failure code for optim: zero is OK
	double Fmin = 0.;          // minimal value of obj fct in optim 
	
	//Setting global variables
	dim_x = *n1; dim_y = *n2; dim_z = *n3;
	siq_init = siq_in; penc = pen; grad = grad_in, ngrad0c = *ngrad0;
	
	// not yet initialized
	double* w = Calloc(maxcompc, double);
	
	//Gradient vectors corresponding to minima in spherical coordinates
	
	if(optmethodc != 4){
		param_length_init=2*maxcompc+1;
//		if(optmethodc != 1){
//			 param_ret = (double *) R_alloc(param_length, sizeof(double));
//		 }
	}else{
		if(methodc == 1){
			param_length_init=3*maxcompc+1;
		}else{
			param_length_init=3*maxcompc+2;
		}
	}
	
	param =  (double *) Calloc(param_length_init, double);
	param_ret = (double *) Calloc(param_length_init, double);
	
	//initialize param
	for(i=0; i<param_length_init; i++){
		 param[i] = 0;
	}
	
	for(i1 = 0; i1 < dim_x; i1++){ for(i2 = 0; i2 < dim_y; i2++){ for(i3 = 0; i3 < dim_z; i3++){
		param_length = param_length_init;
		sigma2_ret[get_ind3d_img(i1,i2,i3)] = sigma2[get_ind3d_img(i1,i2,i3)];
		
		//only analyze voxel within mask
		if (mask[get_ind3d_img(i1,i2,i3)]==1){
			ord=maxcompc+1;
			
			for (j = 0; j < maxcompc; j++){
				iv = siind[(j+2)+(maxcompc+2)*get_ind3d_img(i1,i2,i3)]-1;
				
				if(iv==-1) iv = j; // this should never happen
				
				dir[0]=vert[get_ind2d(0,iv,3)];  //vert dim(max(siind),3)
				dir[1]=vert[get_ind2d(1,iv,3)];
				dir[2]=vert[get_ind2d(2,iv,3)];
				
				paroforient(dir, angles);
				orient_ret[0+2*j+2*maxcompc*get_ind3d_img(i1,i2,i3)] = angles[0];
				orient_ret[1+2*j+2*maxcompc*get_ind3d_img(i1,i2,i3)] = angles[1];
				
				param[2*j+1] = angles[0];
				param[2*j+2] = angles[1];
			}
			
			// initialize EV-parameter
			tmp_int = siind[1+(maxcompc+2)*get_ind3d_img(i1,i2,i3)];
			if(tmp_int > 0){
				param[0] = th[tmp_int-1];
			}else{
				param[0] = 0.001;
			}
			
			sigma_init = sigma2[get_ind3d_img(i1,i2,i3)];
			krit = log(sigma_init) + penIC[0];
			
// 			Rprintf("param from orient:\n");
			
			// working param array
			param_work = Calloc(param_length, double);
			param_last = Calloc(param_length, double);
			for(l = 0; l < param_length; l++){
				param_work[l] = param[l];
				param_last[l] = param[l];
// 				Rprintf(" %f ", param[l]);
			}
			
// 			Rprintf("\n");
			
			// use AIC/ngrad0, BIC/ngrad0 or AICC/ngrad0 respectively
			for(k = maxcompc; k > 0; k--){
				if(k < ord){
					if(k != maxcompc){
						double* param_tmp = Calloc(param_length, double);
						for (l = 0; l < param_length; l++){
							param_tmp[l] = param_work[l];
						}
						if(methodc == 1){
							if(optmethodc != 4){
								param_length=2*k+1;
								/*if(optmethodc != 1){
									param_ret = (double *) R_alloc(param_length, sizeof(double));
								}*/
							}else{
								param_length=3*k+1;
							}
						}else{
							param_length=3*k+2;
						}
						
						param_work = Calloc(param_length,double);
						param_last = Calloc(param_length,double);
						for(l= 0; l < param_length; l++){
							param_work[l] = param_tmp[l];
							param_last[l] = param_tmp[l];
						}
						
						Free(param_tmp);
					}
					
					if(optmethodc == 4){
						int param_red_length;
						if(methodc == 1){
							param_red =  (double *) Calloc(3*k+1, double);
							param_red_length = 3*k+1;
							for(l = 0; l < param_red_length; l++){
								param_red[l] = param_work[l];
							}
							
							if(k == maxcompc){
								for(l=2*k;l<3*k;l++){
									param_red[l]=1/k;
								}
							}else{
								for(l=2*k;l<3*k;l++){
									param_red[l] = w[l-2*k];
								}
							}
						}else{
							param_red =  (double *) Calloc(3*k+2, double);
							param_red_length = 3*k+2;
							for(l = 0; l < param_red_length; l++){
								param_red[l] = param_work[l];
							}
							
							if(k == maxcompc){
								for(l=2*k;l<3*k+1;l++){
									param_red[l]=1/k;
								}
							}else{
								for(l=2*k;l<3*k+1;l++){
									param_red[l] = w[l-2*k];
								}
							}
						}
					}
					
					//neuer check wenn elses drin sind
					
					optimex myoptimpar;
					
					//calculation and estimation of the best fitting model
					
					// initialize ret_val to avoid warnings
					mfunplwghts0_ret ret_val;
					ret_val.order = 0;
					ret_val.lev_par = 0;
					ret_val.lev_sum_w = 0;
					ret_val.mix = NULL;
					ret_val.orient = NULL;
					ret_val.param = NULL;
					ret_val.value = 0;
					ret_val.w = NULL;
					
					double* siq = Calloc(ngrad0c, double);
					for(l=0;l<ngrad0c;l++){
						siq[l] = siq_init[get_ind4d_img(i1,i2,i3,l)];
					}
					
					if(methodc == 1){
						switch (optmethodc){ // R code guarantees method is 1, 2, 3, 4
							case 1: 
								vmmin(param_length, param_work, &Fmin, mfunpl0, gmfunpl0,
								*maxit, trace, mask, abstol, *reltol, nREPORT,
								&myoptimpar, &fncount, &grcount, &fail);
								
								ret_val = mfunplwghts0(param_length, param_work, siq);
								break;
							case 2: 
								cgmin(param_length, param_work, param_ret, &Fmin, mfunpl0,
								gmfunpl0, &fail, abstol, intol, &myoptimpar, 
								type, trace, &fncount, &grcount, *maxit);
								
								ret_val = mfunplwghts0(param_length, param_ret, siq);
								break;
							case 3:
								nmmin(param_length, param_work, param_ret, &Fmin, mfunpl0h,
								&fail, abstol, intol, &myoptimpar, alpha, beta,
								gamma, trace, &fncount, *maxit);
								
								ret_val = mfunplwghts0h(param_length, param_ret, siq);
								break;
							case 4:
// 								lbfgsb(param_length, 
								value = Fmin;
								break; //L-BFGS-B
						}
					
					}else{
						//call optim mixtensoriso method L-BFGS-B
						//other methods seem less numerically stable in this situation
					} 
					
					Free(siq);
					
					if(optmethodc != 4){
						value = ret_val.value;
					}
					
					ord = ret_val.order;
					
					if( (ret_val.lev_par < 0) || (ret_val.lev_sum_w < 0) || (ord < k) )
					{
						ttt = krit;
						
						for(l =0; l < param_length; l++){
							param_work[l] = param_last[l];
						}
					}
					else
					{
						si2new = value/ngrad0c;
						
						ttt = log(si2new) + penIC[ord];
						
						for(l =0; l < param_length; l++){
							param_work[l] = ret_val.param[l];
						}
					}
					
					if(optmethodc == 4){
						w = ret_val.w;
					}
					
					//use directions corresponding to largest weights as initial directions
					if(ttt < krit){    
						krit = ttt;
						
						order_ret[get_ind3d_img(i1,i2,i3)] = ret_val.order;
						
						lev_ret[0+2*get_ind3d_img(i1,i2,i3)] = ret_val.lev_par;
						lev_ret[1+2*get_ind3d_img(i1,i2,i3)] = ret_val.lev_sum_w;
						
						for(l = 0; l < ord; l++){
							mix_ret[l+maxcompc*get_ind3d_img(i1,i2,i3)] = ret_val.mix[l];
							orient_ret[0+2*l+2*maxcompc*get_ind3d_img(i1,i2,i3)] = ret_val.orient[get_ind2d(0,l,2)];
							orient_ret[1+2*l+2*maxcompc*get_ind3d_img(i1,i2,i3)] = ret_val.orient[get_ind2d(1,l,2)];
						}
						for(l = ord; l < maxcompc; l++){
							mix_ret[l+maxcompc*get_ind3d_img(i1,i2,i3)] = 0;
						}
						
						sigma2_ret[get_ind3d_img(i1,i2,i3)] = si2new;
					} // end if(ttt < krit)
				} // end if(k < ord)
			} // end model order
		} // end mask
		R_CheckUserInterrupt();
	      } // end i3
	      R_CheckUserInterrupt();
	    } // end i2
	} // end i1
	
	Free(w);
	Free(param);
	Free(param_ret);
	Free(param_work);
	Free(param_red);
}

void mixturetest(double* test2d, double *test3d, double *test4d, double* test5d,
		 double* testimgback, int* x, int* y, int* z){
  int i;
  dim_x = *x; dim_y = *y; dim_z = *z;

  Rprintf("Testimgback\n");
 Rprintf("test4d[%i]=%g\n",get_ind4d_img(1,2,3,4),test4d[3+4*get_ind3d_img(4,4,4)]);
  Rprintf("test4d[%i]=%g\n",get_ind4d_img(3,2,4,1),test4d[1+4*get_ind3d_img(2,5,6)]);
  Rprintf("test4d[1,5,3,1]=%g\n",test4d[3+4*get_ind3d_img(4,5,6)]);

  
  Rprintf("Test2d:\n");
  Rprintf("test2d[1,2]=%g\n",test2d[get_ind2d(3,4,5)]);
  Rprintf("test2d[3,2]=%g\n",test2d[get_ind2d(2,0,5)]);
  Rprintf("test2d[1,5]=%g\n",test2d[get_ind2d(4,4,5)]);

  Rprintf("Test3d:\n");
  Rprintf("test3d[1,2,3]=%g\n",test3d[get_ind3d_img(3,1,2)]);
  Rprintf("test3d[3,2,1]=%g\n",test3d[get_ind3d_img(2,2,4)]);
  Rprintf("test3d[1,5,3]=%g\n",test3d[get_ind3d_img(2,0,2)]);
  for(i=0;i<10;i++){
    Rprintf("%g, ", test4d[i]);
    
  }
  
  Rprintf("Test4d:\n");

  Rprintf("test4d[%i]=%g\n",get_ind4d_img(1,2,3,4),test4d[get_ind4d_img(4,4,4,4)]);
  Rprintf("test4d[%i]=%g\n",get_ind4d_img(3,2,4,1),test4d[get_ind4d_img(1,2,6,3)]);
  Rprintf("test4d[1,5,3,1]=%g\n",test4d[get_ind4d_img(4,5,6,3)]);

  Rprintf("Test5d:\n");
  Rprintf("test5d[%i]=%g\n",(4+3*get_ind4d_img(1,2,3,1)),test5d[get_ind5d_img(2,4,5,2,4,3)]);

  Rprintf("test5d[3,2,1,2,3]=%g\n",test5d[get_ind5d_img(3,5,6,2,2,3)]);
  Rprintf("test5d[1,5,3,3,2]=%g (%i)\n",(get_ind5d_img(2,1,4,1,0,3)),test5d[get_ind5d_img(2,3,6,1,2,3)]);
  //  int count =0;

  // PROTECT(&i);
  // PROTECT(&count);
  /*
  for (i=0;i<20;i++){
    Rprintf("sigma2 step %i: %g\n",i, sigma2[i]);
        Rprintf("mask step %i: %g\n",i, mask[i]);
    Rprintf("siind step %i: %g\n",i, siind[i]);
  }

  Rprintf("\n");
  */
  //  for (i=20*20*10;i<30*30*20;i++){
  //  Rprintf("%i", siind[i]);
  // }
  
  /*
  
  Rprintf("printing mix C before op\n\n");
	
  for(i=0; i<1000;i++){
    if (i-count==200){
      Rprintf("step %i of %i\n", i, *l_m);
      count=i;
    }
    Rprintf("%i", mix[i]);
    }*/
  /*
  Rprintf("adding mix\n\n");
  i=0;
  count=0;
  for(i=0; i<1000;i++){
    if (i-count==200){
      Rprintf("step %i of %i\n", i, *l_m);
      count=i;
    }

    mix[i]=mask[i]+siind[i];
     R_CheckUserInterrupt();
  }

  for(i=0; i<1000;i++){
    if (i-count==200){
      Rprintf("step %i of %i\n", i, *l_m);
      count=i;
    }
    Rprintf("%i", mix[i]);
   }
   R_CheckUserInterrupt();
  
    
   //   UNPROTECT(2);

   */
	

}
/*
void convolveCall(SEXP a, SEXP b)

{
  int i;
  double* test=REAL(a);
  double ra,rb;

  PROTECT(test);
  PROTECT(ra);
  PROTECT(rb);
  
  for(i=0;i<4;i++){
    ra=REAL(a)[i];
    Rprintf("a in %i: %d", i, ra);
     Rprintf("a in %i: %d", i, test[i]);
  }
  for(i=0;i<4;i++) Rprintf("b in %i: %d", i, REAL(b)[i]);

  UNPROTECT(3);

  } */

