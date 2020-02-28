#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#include <Rinternals.h> // for SEXP
#include <R_ext/RS.h>

SEXP interface_tracking(SEXP data_dir_coords, SEXP data_FA_values, SEXP data_mask,
			        SEXP dim_x, SEXP dim_y, SEXP dim_z, SEXP roi_x_s, SEXP roi_x_e,
			        SEXP roi_y_s, SEXP roi_y_e, SEXP roi_z_s, SEXP roi_z_e,
			        SEXP voxelext_x, SEXP voxelext_y, SEXP voxelext_z,
			        SEXP min_anisotropy, SEXP max_angle);

SEXP interface_tracking_mixtensor(
             	SEXP data_dir_coords, SEXP data_order,
             	SEXP data_FA_values, SEXP data_mask,
             	SEXP data_mix, SEXP maxorder,
             	SEXP dim_x, SEXP dim_y, SEXP dim_z,
             	SEXP roi_x_s, SEXP roi_x_e,
             	SEXP roi_y_s, SEXP roi_y_e,
             	SEXP roi_z_s, SEXP roi_z_e,
             	SEXP voxelext_x, SEXP voxelext_y, SEXP voxelext_z,
             	SEXP min_anisotropy,
             	SEXP max_angle);

static const R_CallMethodDef callMethods[]  = {
  {"interface_tracking", (DL_FUNC) &interface_tracking, 17},
  {"interface_tracking_mixtensor", (DL_FUNC) &interface_tracking_mixtensor, 20},
  {NULL, NULL, 0}
};

extern void dtens(int* n1, double* param, double* sig_in, int* ngrad,
  double* btb_in, double* sdcoef, double* sig_tmp, double* vinv_tmp,
  int* maxit, double* reltol);
extern void mixture( int* n1, int* siind, int* ngrad0, int* maxcomp,
  int* maxit, double* pen, double* grad_in, double* reltol, double* th,
  double* penIC, double* sigma2, double* vert, double* siq_in,
  double* sigma2_ret, double* orient_ret, int* order_ret, double* lev_ret,
  double* mix_ret );
extern void mixtrl0b( int* n1, int* siind, double* wi, int* ngrad, int* maxcomp,
  int* maxit, double* grad_in, double* bv_in, double* lambda_in,
  double* alpha_in, double* factr, double* penIC, double* sigma2,
  double* vert, double* si_in, double* sigma2_ret, double* orient_ret,
  int* order_ret, double* mix_ret);
extern void mixtrl1b( int* n1, int* siind, double* wi, int* ngrad, int* maxcomp,
  int* maxit, double* grad_in, double* bv_in, double* lambda_in,
  double* alpha_in, double* factr, double* penIC, double* sigma2,
  double* vert, double* si_in, double* sigma2_ret, double* orient_ret,
  int* order_ret, double* lambda_ret, double* mix_ret);
extern void mixtrl2b( int* n1, int* siind, double* wi, int* ngrad, int* maxcomp,
  int* maxit, double* grad_in, double* bv_in, double* lambda_in,
  double* alpha_in, double* factr, double* penIC, double* sigma2,
  double* vert, double* si_in, double* sigma2_ret, double* orient_ret,
  int* order_ret, double* alpha_ret, double* lambda_ret, double* mix_ret);

static R_NativePrimitiveArgType dtens_t[]={INTSXP, REALSXP, REALSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType mixture_t[]={INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType mixtrl0b_t[]={INTSXP, INTSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType mixtrl1b_t[]={INTSXP, INTSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType mixtrl2b_t[]={INTSXP, INTSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP,
  REALSXP};
static const R_CMethodDef CMethods[] = {
            {"dtens", (DL_FUNC) &dtens, 10, dtens_t},
            {"mixture", (DL_FUNC) &mixture, 18, mixture_t},
            {"mixtrl0b", (DL_FUNC) &mixtrl0b, 19, mixtrl0b_t},
            {"mixtrl1b", (DL_FUNC) &mixtrl1b, 20, mixtrl1b_t},
            {"mixtrl2b", (DL_FUNC) &mixtrl2b, 21, mixtrl2b_t},
            {NULL, NULL, 0,NULL}
};

void F77_NAME(adcradii)( double* vert, int* nv, double* tens, int* ntens,
  double* radii);
void F77_NAME(awsrgdti)( double* , double* , double* , int* , int* ,
  int* , int* , int* , int* , double* , double* , double* , double* ,
  double* , double* , double* , double* , double* , double* , double* , double* ,
  int* , double* , double* , double* , double* , double* , double* , double* ,
  double* , int* , int* , int* , double* , int* , double* , double* );
void F77_NAME(awssidti)( double* s0, double* si, int* mask, double* th,
  double* bi, double* ani, double* andir, double* det, double* bcov,
  double* solvebtb, double* sigma2, double* sigma2h, int* n1, int* n2,
  int* n3, int* ngrad, double* h, double* vext, double* rho, double* lambda,
  double* thnew, double* sigma2n, double* swsi, double* eps);
void F77_NAME(cfibers)( double* fibers, int* sind, int* nf, int* nsi,
  double* delta, int* nnf);
void F77_NAME(d2rall)( double* D, double* rho, int* nvox);
void F77_NAME(dti2dga)( double* D, int* n, int* mask, double* ga, double* md,
  double* adir);
void F77_NAME(dti2dfa)( double* D, int* n, int* mask, double* fa, double* md,
  double* adir);
void F77_NAME(dti3dall)( double* D, int* n, double* fa, double* ga,
  double* md, double* adir, double* ev);
void F77_NAME(dti3dand)( double* D, int* n, double* andir);
void F77_NAME(dti3devall)( double* D, int* n, double* andir, double* evalues);
void F77_NAME(dti3dev)( double* D, int* n, double* ev);
void F77_NAME(dti3dreg)( double* D, int* n);
void F77_NAME(dtieigen)( double* D, int* n, double* fa, double* ev,
  double* adir);
void F77_NAME(dtiind3d)( double* D, int* n, double* fa, double* ga, double* md,
  double* adir, double* bary);
void F77_NAME(ellradii)( double* vert, int* nv, double* tens, int* ntens,
  double* radii);
void F77_NAME(fibersta)( double* fibers, int* nfibers, int* start,
  int* nstart);
void F77_NAME(gethani)( double* x, double* y, double* value, double* a,
  double* vext, double* eps, double* bw);
void F77_NAME(getmask)( double* s0, int* n1, int* n2, int* n3, int* ns,
  double* level, int* msize, double* prob, double* s0m, int* mask);
void F77_NAME(getsii30)( double* si, double* vsi, int* ngrad,
  int* nvox, int* m, double* dgrad, int* nv, double* th, int* nth, int* indth,
  double* egrad, int* isample, int* ntry, double* sms, double* z, int* siind,
  double* mval, int* ns, int* mask);
void F77_NAME(getsii31)( double* si, double* vsi, int* ngrad, int* nvox,
  int* m, double* dgrad, int* nv, int* iandir, double* th, int* nth,
  int* indth, double* egrad, int* isample, int* ntry, double* sms, double* z,
  int* siind, double* mval, int* ns, int* mask, double* dgradv, double* maxc);
void F77_NAME(getsiibv)( double* si, int* ngrad, int* nvox,
  int* m, double* dgrad, double* bv, int* nv, double* alpha, double* lambda,
  double* egrad, int* isample, int* ntry, double* sms, double* z0, double* z,
  int* siind, double* wi, double* mval, int* ns);
void F77_NAME(getvofh)( double* a, double* bw, double* vext, double* vol);
void F77_NAME(hg1f1)( double* a, double* b, double* z, int* n, double* fz);
void F77_NAME(iandir)( double* vico, int* nvico, double* andir, int* nvox,
  int* landir, int* iandi);
void F77_NAME(initdata)( double* si, int* n1, int* n2, int* n3, int* nb,
  double* maxvalue);
void F77_NAME(lconnect)( int* segm, int* n1, int* n2, int* n3, int* i1,
  int* i2, int* i3, int* ind1, int* ind2, int* ind3, int* mask);
void F77_NAME(mcorr)( double* res, int* mask, int* n1, int* n2, int* n3,
  int* nv, double* sigma, double* mean, double* scorr, int* l1, int* l2,
  int* l3);
void F77_NAME(means0)( double* s0, int* n, int* ng0, int* level, double* ms0,
  int* mask);
void F77_NAME(mixandir)( double* ori, double* mix, int* ord, int* mo,
  int* nobj, double* andir);
void F77_NAME(mixtradi)( double* vert, int* nv, double* ev, double* ori,
  double* mix, int* ord, int* mo, int* nobj, double* radii);
void F77_NAME(odfradii)( double* vert, int* nv, double* tens, int* ntens,
  double* radii);
void F77_NAME(outlier)( double* si, int* n, int* nb, int* s0ind, int* siind,
  int* ls0, double* sinew, int* ind);
void F77_NAME(outlierp)( double* si, int* n, int* nb, int* s0ind, int* ls0,
  int* siind, int* lsi, double* sinew, int* nb1);
void F77_NAME(pgtsii30)( double* si, double* vsi, int* ngrad,
  int* nvox, int* m, double* dgrad, int* nv, double* th, int* nth, int* indth,
  double* egrad, int* isample, int* ntry, double* sms, double* z, int* siind,
  double* mval, int* ns);
void F77_NAME(pgtsii31)( double* si, double* vsi, int* ngrad, int* nvox,
  int* m, double* dgrad, int* nv, int* iandir, double* th, int* nth,
  int* indth, double* egrad, int* isample, int* ntry, double* sms, double* z,
  int* siind, double* mval, int* ns, double* dgradv, double* maxc);
void F77_NAME(projdt2)( double* th, int* n1, int* n2, int* n3, double* thnew,
   double* ani, double* dir, double* det, double* eps);
void F77_NAME(r2dall)( double* rho, double* D, int* nvox);
void F77_NAME(rho2d0)( double* rho, double* D);
void F77_NAME(reducefe)( double* fibers, int* nsegm, int* startf, int* endf,
  int* nfibers, int* keep, double* maxd);
void F77_NAME(reducefi)( double* fibers, int* nsegm, int* startf, int* endf,
  int* nfibers, int* keep, double* maxd);
void F77_NAME(roifiber)( double* fiber, double* newfiber, int* sizef,
  int* ifiber, int* mlf, int* startf, int* lengthf, int* nfibers, int* roi,
  int* nx, int* ny, int* nz, double* vext, int* sizenf, int* nnfiber);
void F77_NAME(selisamp)( int* isample, int* nguess, int* maxcomp,
  double* dgrad, int* ndg, int* ind, double* maxc);
void F77_NAME(smsigma)( double* sigma2, int* n1, int* n2, int* n3, double* h,
  double* vext, double* sigma2h);
void F77_NAME(sweeps0)( double* si, double* s0, int* n, int* ng0, int* ng1,
  int* level, double* siq, double* ms0, double* vsi, int* mask);
void F77_NAME(sweeps0p)( double* si, double* s0, int* n, int* ng0, int* ng1,
  int* level, double* siq, int* ng2);
void F77_NAME(tensres)( double* th0, double* D, double* s, int* nvox, int* nb,
  double* b, double* res, double* rss);
void F77_NAME(thcorr)( double* w, int* n1, int* n2, int* n3, double* scorr,
  int* l1, int* l2, int* l3);
void F77_NAME(touchfi)( double* fibers1, int* nsegm1, int* startf1, int* endf1,
  int* nfibers1, int* keep, double* fibers2, int* nsegm2, double* maxdist);

static R_NativePrimitiveArgType adcradii_t[]={REALSXP, INTSXP, REALSXP,
  INTSXP, REALSXP};
static R_NativePrimitiveArgType awsrgdti_t[]={REALSXP, REALSXP, REALSXP,
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType awssidti_t[]={REALSXP, REALSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cfibers_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType d2rall_t[]={REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType dti2dga_t[]={REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType dti2dfa_t[]={REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType dti3dall_t[]={REALSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType dti3dand_t[]={REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType dti3devall_t[]={REALSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType dti3dev_t[]={REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType dti3dreg_t[]={REALSXP, INTSXP};
static R_NativePrimitiveArgType dtieigen_t[]={REALSXP, INTSXP, REALSXP,
  REALSXP, REALSXP};
static R_NativePrimitiveArgType dtiind3d_t[]={REALSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType ellradii_t[]={REALSXP, INTSXP, REALSXP,
  INTSXP, REALSXP};
static R_NativePrimitiveArgType fibersta_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP};
static R_NativePrimitiveArgType gethani_t[]={REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType getmask_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType getsii30_t[]={REALSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType getsii31_t[]={REALSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType getsiibv_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType getvofh_t[]={REALSXP, REALSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType hg1f1_t[]={REALSXP, REALSXP, REALSXP,
  INTSXP, REALSXP};
static R_NativePrimitiveArgType iandir_t[]={REALSXP, INTSXP, REALSXP,
  INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType initdata_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType lconnect_t[]={INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType mcorr_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType means0_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, INTSXP};
static R_NativePrimitiveArgType mixandir_t[]={REALSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType mixtradi_t[]={REALSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType odfradii_t[]={REALSXP, INTSXP, REALSXP,
  INTSXP, REALSXP};
static R_NativePrimitiveArgType outlier_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType outlierp_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType pgtsii30_t[]={REALSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType pgtsii31_t[]={REALSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType projdt2_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType r2dall_t[]={REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType rho2d0_t[]={REALSXP, REALSXP};
static R_NativePrimitiveArgType reducefe_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType reducefi_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType roifiber_t[]={REALSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType selisamp_t[]={INTSXP, INTSXP,
  INTSXP, REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType smsigma_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType sweeps0_t[]={REALSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType sweeps0p_t[]={REALSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType tensres_t[]={REALSXP, REALSXP, REALSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType thcorr_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType touchfi_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP};

static const R_FortranMethodDef FMethods[] = {
            {"adcradii", (DL_FUNC) &adcradii_ , 5,adcradii_t},
            {"awsrgdti", (DL_FUNC) &awsrgdti_ , 37, awsrgdti_t},
            {"awssidti", (DL_FUNC) &awssidti_ , 24, awssidti_t},
            {"cfibers", (DL_FUNC) &cfibers_ , 6, cfibers_t},
            {"d2rall", (DL_FUNC) &d2rall_ , 3, d2rall_t},
            {"dti2dga", (DL_FUNC) &dti2dga_ , 6, dti2dga_t},
            {"dti2dfa", (DL_FUNC) &dti2dfa_ , 6, dti2dfa_t},
            {"dti3dall", (DL_FUNC) &dti3dall_ , 7, dti3dall_t},
            {"dti3dand", (DL_FUNC) &dti3dand_ , 3, dti3dand_t},
            {"dti3devall", (DL_FUNC) &dti3devall_ , 4, dti3devall_t},
            {"dti3dev", (DL_FUNC) &dti3dev_ , 3, dti3dev_t},
            {"dti3dreg", (DL_FUNC) &dti3dreg_ , 2, dti3dreg_t},
            {"dtieigen", (DL_FUNC) &dtieigen_ , 5, dtieigen_t},
            {"dtiind3d", (DL_FUNC) &dtiind3d_ , 7, dtiind3d_t},
            {"ellradii", (DL_FUNC) &ellradii_ , 5, ellradii_t},
            {"fibersta", (DL_FUNC) &fibersta_ , 4, fibersta_t},
            {"gethani", (DL_FUNC) &gethani_ , 7, gethani_t},
            {"getmask", (DL_FUNC) &getmask_ , 10, getmask_t},
            {"getsii30", (DL_FUNC) &getsii30_ , 19, getsii30_t},
            {"getsii31", (DL_FUNC) &getsii31_ , 22, getsii31_t},
            {"getsiibv", (DL_FUNC) &getsiibv_ , 19, getsiibv_t},
            {"getvofh", (DL_FUNC) &getvofh_ , 4, getvofh_t},
            {"hg1f1", (DL_FUNC) &hg1f1_ , 5, hg1f1_t},
            {"iandir", (DL_FUNC) &iandir_ , 6, iandir_t},
            {"initdata", (DL_FUNC) &initdata_ , 6, initdata_t},
            {"lconnect", (DL_FUNC) &lconnect_ , 11, lconnect_t},
            {"mcorr", (DL_FUNC) &mcorr_ , 12, mcorr_t},
            {"means0", (DL_FUNC) &means0_ , 6, means0_t},
            {"mixandir", (DL_FUNC) &mixandir_ , 6, mixandir_t},
            {"mixtradi", (DL_FUNC) &mixtradi_ , 9, mixtradi_t},
            {"odfradii", (DL_FUNC) &odfradii_ , 5, odfradii_t},
            {"outlier", (DL_FUNC) &outlier_ , 8, outlier_t},
            {"outlierp", (DL_FUNC) &outlierp_ , 9, outlierp_t},
            {"pgtsii30", (DL_FUNC) &pgtsii30_ , 18, pgtsii30_t},
            {"pgtsii31", (DL_FUNC) &pgtsii31_ , 21, pgtsii31_t},
            {"projdt2", (DL_FUNC) &projdt2_ , 9, projdt2_t},
            {"r2dall", (DL_FUNC) &r2dall_ , 3, r2dall_t},
            {"rho2d0", (DL_FUNC) &rho2d0_ , 2, rho2d0_t},
            {"reducefe", (DL_FUNC) &reducefe_ , 7, reducefe_t},
            {"reducefi", (DL_FUNC) &reducefi_ , 7, reducefi_t},
            {"roifiber", (DL_FUNC) &roifiber_ , 15, roifiber_t},
            {"selisamp", (DL_FUNC) &selisamp_ , 7, selisamp_t},
            {"smsigma", (DL_FUNC) &smsigma_ , 7, smsigma_t},
            {"sweeps0", (DL_FUNC) &sweeps0_ , 10, sweeps0_t},
            {"sweeps0p", (DL_FUNC) &sweeps0p_ , 8, sweeps0p_t},
            {"tensres", (DL_FUNC) &tensres_ , 8, tensres_t},
            {"thcorr", (DL_FUNC) &thcorr_ , 8, thcorr_t},
            {"touchfi", (DL_FUNC) &touchfi_ , 9, touchfi_t},
            {NULL, NULL, 0,NULL}
};

void R_init_dti(DllInfo *dll)
         {
             R_registerRoutines(dll, CMethods, callMethods, FMethods , NULL);
             R_useDynamicSymbols(dll,FALSE);
             R_forceSymbols(dll,TRUE);
         }
