#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>


#include <sys/time.h>
#include <omp.h>

typedef struct
{
  int ngrad;
  double *siq;
  double *grad;
  int p;
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

// double F77_NAME(dotprod3)(double *, double *);
// double z1 = F77_CALL(dotprod3)(dir, test);
//void F77_NAME(nnls)(double *, int, int, int, double *, double *, double *, double *, double *, int * , int *);
void F77_NAME(nnls)(double *a, int *mda, int *m, int *n, double *b, double *X, double *rnorm, double *w, double *zz, int *index, int *mode);

void paroforient(double *dir, double *angles)
{
  angles[0] = acos(dir[2]);
  double sth = sin(angles[0]);
  angles[1] = 0.;

  if (sth < 1e-8)
  {
    angles[0] = 0.;
  }
  else
  {
    double z = dir[0]/sth;
    angles[1] = (abs(z) >= 1.) ? ((z<0.) ? 0. : M_PI) : acos(z)*sign(dir[1]);
    if (angles[1] < 0.) angles[1] += M_2PI;
  }
  return;
}

static R_INLINE int get_ind2d(int i, int j, int n1)
{
  return i + j*n1;
}

static R_INLINE int get_ind3d(int i, int j, int k, int n1, int n2)
{
  return i + j*n1 + k*n1*n2;
}

void getweights(int m, double *cw, double *z)
{
  double ew;
  double sw = 0.0;

  for (int i = 0; i < m; i++)
  {
    // use log and exp to have free optim problem
    ew = exp(cw[i]);
    z[i] = (1. - sw) * ew/(1.+ew);
    sw += z[i];
  }
  z[m] = 1. - sw;
  return;
}

double fn1(int n, double *par, void *ex)
{
  optimex ext = *((optimex*)ex);

  int m = (n-1)/3;     // order of mix tensor model
  double *z = (double *) R_alloc(ext.ngrad, sizeof(double));
  for (int k = 0; k < ext.ngrad; k++)
    z[k] = 0;
  double dir[3];       // used for direction vector calculated from angles in par
  double c1 = exp(par[0]), c2 = exp(par[1]), sw = 0, w, ew, sth, z1 = 0;
  int i, j;

  double erg = 0; // result

  for (i = 0; i < m; i++)
  {
    if (i == m-1 )
    {
      w = 1.0 - sw;
    }
    else
    {
      // use log and exp to have free optim problem
      ew = exp(par[3*i+4]);
      w = (1.0 - sw)*ew/(1.0+ew);
      sw += w;
    }
    sth = sin(par[i*3+2]);
    dir[0] = sth*cos(par[i*3+3]);
    dir[1] = sth*sin(par[i*3+3]);
    dir[2] = cos(par[i*3+2]);
    for (j = 0; j < ext.ngrad; j++)
    {
      z1 = dir[0]*ext.grad[j] + dir[1]*ext.grad[ext.ngrad+j] + dir[2]*ext.grad[2*ext.ngrad+j]; 
      z[j] += (ext.p == 0) ? w*exp(-c2 - c1*z1*z1) : w*exp(-ext.p*log(1 + ( c2 + c1*z1*z1 )/ext.p));
    }
  }

  // calculate residual variance
  for (j = 0; j < ext.ngrad; j++)
  {
    z1 = ext.siq[j] - z[j];
    erg += z1*z1;
  }

  // finished
  return erg;
}

double fnpl(int n, double *par, void *ex)
{
  optimexpl ext = *((optimexpl*)ex);

  int m = (n-1)/2;     // order of mix tensor model
  for (int k = 0; k < ext.ngrad; k++)
  {
    for (int kk = 0; kk < m; kk++)
      ext.z[kk*ext.ngrad + k] = 0;
    ext.fv[k] = ext.siq[k];
  }
  double dir[3];       // used for direction vector calculated from angles in par
//  double c1 = exp(par[0]), sth, z1;
  double c1 = par[0], sth, z1;
  int i, i2, j;
  int ind[10], mode = 0;
  double work2[10];

  double erg = 0; // result

  for (i = 0; i < m; i++)
  {
    i2 = 2*i;
    sth = sin(par[i2+1]);
    dir[0] = sth*cos(par[i2+2]);
    dir[1] = sth*sin(par[i2+2]);
    dir[2] = cos(par[i2+1]);
    for (j = 0; j < ext.ngrad; j++)
    {
      z1 = dir[0]*ext.grad[j] + dir[1]*ext.grad[ext.ngrad+j] + dir[2]*ext.grad[2*ext.ngrad+j]; 
      ext.z[j + i*ext.ngrad] += exp(-c1*z1*z1);
    }
  }
 // siq will be replaced, need to copy it if C-version of optim is used
  F77_CALL(nnls)(ext.z, &ext.ngrad, &ext.ngrad, &m, ext.fv, ext.w, &erg, work2, ext.work1, ind, &mode);

  // finished
  return erg;
}

void grpl(int n, double *par, double *gr, void *ex)
{
  optimexpl ext = *((optimexpl*)ex);

  int m = (n-1)/2;     // order of mix tensor model
  for (int k = 0; k < ext.ngrad; k++)
  {
    for (int kk = 0; kk < m; kk++)
      ext.z[kk*ext.ngrad + k] = 0;
    ext.fv[k] = ext.siq[k]; // zwischenspeicher fuer siq in nnls, spaeter neu initialisiert
  }
  int i, i2, j;
  int ind[10], mode = 0;
  double dir[3];       // used for direction vector calculated from angles in par
  double c1 = par[0], sth, cpsi, spsi, z1;
  double work2[10];
  double erg = 0; // result

  for (i = 0; i < m; i++)
  {
    i2 = 2*i;
    sth = sin(par[i2+1]);
    cpsi = cos(par[i2+2]);
    spsi = sin(par[i2+2]);
    dir[0] = sth*cos(par[i2+2]);
    dir[1] = sth*sin(par[i2+2]);
    dir[2] = cos(par[i2+1]);
    for (j = 0; j < ext.ngrad; j++)
    {
      z1 = dir[0]*ext.grad[get_ind2d(j, 0, ext.ngrad)] + dir[1]*ext.grad[get_ind2d(j, 1, ext.ngrad)] + dir[2]*ext.grad[get_ind2d(j, 2, ext.ngrad)]; 
      ext.qv[get_ind2d(i, j, m)] = z1*z1;
      ext.z[get_ind2d(j, i, ext.ngrad)] += exp(-c1*ext.qv[get_ind2d(i, j, m)]);
      ext.dqv[get_ind3d(0, i, j, 2, m)] = 2.0*(dir[2]*(cpsi*ext.grad[get_ind2d(j, 0, ext.ngrad)] + 
                                      spsi*ext.grad[get_ind2d(j, 1, ext.ngrad)])-sth*ext.grad[get_ind2d(j, 2, ext.ngrad)])*ext.qv[get_ind2d(i, j, m)];
      ext.dqv[get_ind3d(1, i, j, 2, m)] = 2.0*sth*(cpsi*ext.grad[get_ind2d(j, 1, ext.ngrad)]-spsi*ext.grad[get_ind2d(j, 0, ext.ngrad)])*ext.qv[get_ind2d(i, j, m)];
    }
  }

  F77_CALL(nnls)(ext.z, &ext.ngrad, &ext.ngrad, &m, ext.fv, ext.w, &erg, work2, ext.work1, ind, &mode);

  for (j = 0; j < ext.ngrad; j++)
  {
    ext.fv[j] = 0.0;
    ext.dfv[get_ind2d(0,j,n)] = 0.0;
    for (i = 0; i < m; i++)
    {
      z1 = exp(-c1*ext.qv[get_ind2d(i, j, m)]);
      ext.fv[j] += ext.w[i]*z1;
      ext.dfv[get_ind2d(0, j, n)] -= ext.w[i]*ext.qv[get_ind2d(i, j, m)]*z1;
      ext.dfv[get_ind2d(2*i+1, j, n)] = -ext.w[i]*c1*ext.dqv[get_ind3d(0, i, j, 2, m)]*z1;
      ext.dfv[get_ind2d(2*i+2, j, n)] = -ext.w[i]*c1*ext.dqv[get_ind3d(1, i, j, 2, m)]*z1;
    }
  }

  for (i = 0; i < n; i++)
  {
    gr[i] = 0.0;
    for (j = 0; j < ext.ngrad; j++)
    {
       gr[i] = gr[i] - ext.dfv[get_ind2d(i, j, n)]*(ext.siq[j]-ext.fv[j]);
    }
    gr[i] = 2.0*gr[i];
//    Rprintf("%f ", gr[i]);
  }
//  Rprintf("\n");
}


void mixture(int *method, int *r, int *mask, double *siq, int *siind, int *n, double *grad, int *maxcomp, int *ep, int *maxit, double *reltol,
             double *order, double *lev, double *mix, double *orient, double *p, double *sigma2){

  int nv = *r, ngrad = *n, mc = *maxcomp;
  int iv, mc0, lpar;
  double rss, krit, ttt;
  double *siiq = (double *) R_alloc(ngrad, sizeof(double));
  int fail;              // failure code for optim: zero is OK
  int fncount;           // number of calls to obj fct in optim
  double Fmin = 0.;          // minimal value of obj fct in optim
  int i, k;
  int gradind;
  double angles[2];
  double dir[3];
  double par[11];
  double cpar[11];
  double x[11];
  double tmp;
  double theta, phi;
//  struct timeval tp1, tp2;
//  struct timezone tzp;
//  double timm = 0;

// #pragma omp parallel for private(i, siiq, mc0, gradind, tmp, par, dir, angles, rss, krit, k, lpar, cpar, x, fail, Fmin, ttt, fncount, reltol, maxit)
  for (iv = 0; iv < nv; iv++)
  {
    if (mask[iv] != 1)
    {
      order[iv] = 0;
      lev[get_ind2d(0, iv, 2)] = 0.;
      lev[get_ind2d(1, iv, 2)] = 0.;
      for (i = 0; i < mc; i++)
      {
        mix[get_ind2d(i, iv, mc)] = 0;
        orient[get_ind3d(0, i, iv, 2, mc)] = 0.;
        orient[get_ind3d(1, i, iv, 2, mc)] = 0.;
      }
//          if (*method == 3) p[get_ind3d(i1, i2, i3, n1, n2)] = 0;
      sigma2[iv] = 0.;
    }
    else
    {
      // prepare signal
      for (i = 0; i < ngrad; i++)
      {
        siiq[i] = siq[get_ind2d(iv, i, nv)];
      }
      // determine useful prime estimates
      mc0 = siind[get_ind2d(0, iv, mc+1)]; // order of modell 
      gradind = siind[get_ind2d(1, iv, mc+1)];
      tmp = -log(siq[get_ind2d(iv, gradind-1, nv)]);
      par[0] = log(tmp*.8);
      par[1] = log(tmp*.2);
      for (i = 0; i < mc0; i++)
      {
        gradind = siind[get_ind2d(i+1, iv, mc+1)];
        dir[0] = grad[get_ind2d(gradind-1, 0, ngrad)];
        dir[1] = grad[get_ind2d(gradind-1, 1, ngrad)];
        dir[2] = grad[get_ind2d(gradind-1, 2, ngrad)];
        paroforient(dir, angles);    // INDEX????

        par[3*i + 2] = angles[0];
        par[3*i + 3] = angles[1]; 
//            Rprintf("par[] %i, %f %f\n", 3*i+2, par[3*i+2], par[3*i+3]); 
      }
//         if (*method == 3) par[lpar] = 0;

      // estimate models and select
      rss = R_PosInf;
      krit = R_PosInf;
      for (k = mc0; k > 0; k--)
      {
        for (i = 0; i < 3*k+1; i++) cpar[i] = par[i];
       //if (*method == 3) cpar[3*k+1] = ??;
        lpar = (*method == 3) ? 3*k+2 : 3*k+1;
        // use log and exp to have free optim problem
        if (k>1) for (i = 0; i < k-1; i++) cpar[3*i+4] = -log(k-i-1.); // par[3*(2:k)-1] <- -log((k-1):1)
        optimex myoptimpar;
        myoptimpar.ngrad = ngrad;
        myoptimpar.siq = siiq;
        myoptimpar.grad = grad;
        switch (*method)
        { // R code guarantees method is 1, 2, 3
          case 1:
            myoptimpar.p = 0; // unused here
            nmmin(lpar, cpar, x, &Fmin, fn1,
                  &fail, R_NegInf, *reltol, &myoptimpar,
                  1.0, 0.5, 2.0, 0,
                  &fncount, *maxit);
            break;
          case 2:
            myoptimpar.p = *ep; // exp for Jian model
            nmmin(lpar, cpar, x, &Fmin, fn1,
                  &fail, R_NegInf, *reltol, &myoptimpar,
                  1.0, 0.5, 2.0, 0,
                  &fncount, *maxit);
            break;
//          case 3:
//            myoptimpar.p = 0; // unused here
//            nmmin(lpar, par, x, Fmin, fn3,
//                  &fail, R_NegInf, *reltol, &myoptimpar,
//                  1.0, 0.5, 2.0, 0,
//                  &fncount, *maxit);
        }

        if (Fmin < rss) rss = Fmin;
        ttt = Fmin + (6.*k+2.)/(ngrad - 3.*mc - 1.) * rss;
        if (ttt < krit)
        {
          krit = ttt;
          order[iv] = k;
          lev[get_ind2d(0, iv, 2)] = x[0];
          lev[get_ind2d(1, iv, 2)] = x[1];
          if (k == 1)
          {
            mix[get_ind2d(0, iv, mc)] = 1.;
          }
          else
          {
            double *zmm = (double *) R_alloc(k, sizeof(double));
            double *zm = (double *) R_alloc(k-1, sizeof(double));
            for (i = 0; i < k-1; i++) zm[i] = x[3*i+4];
            getweights(k-1, zm, zmm);
            for (i = 0; i < k; i++) mix[get_ind2d(i, iv, mc)] = zmm[i];
          }

          for (i = 0; i < k; i++)
          {
            theta = x[3*i+2];
            while (theta < 0.) 
            {
            //  Rprintf("theta %f %f", theta, M_PI); 
              theta += M_PI;
            //  Rprintf("theta %f\n", theta); 
            }
            while (theta > M_PI) theta -= M_PI;
            phi = x[3*i+3];
            while (phi < 0.) phi += M_2PI;
            while (phi > M_2PI) phi -= M_2PI;
            orient[get_ind3d(0, i, iv, 2, mc)] = theta;
            orient[get_ind3d(1, i, iv, 2, mc)] = phi;
          }
//         if (*method == 3) p[i1*n2*n3 + i2*n3 + i3] = x[lpar-1];
        }
      }
      sigma2[iv] = rss/(ngrad-3.*mc0-1.);
    }
    R_CheckUserInterrupt();
  }
//  gettimeofday(&tp1, &tzp);
//  gettimeofday(&tp2, &tzp);
//  if (tp1.tv_usec > tp2.tv_usec) tp2.tv_usec += 1000000;
//  timm += tp2.tv_usec - tp1.tv_usec;
//  Rprintf("zeit %f\n", timm/1000000.);
  return;
}

void mixturepl(int *method,
               int *r,
               int *mask,
               double *siq,
               int *siind,
               int *n,
               double *grad,
               int *maxcomp,
               int *ep,
               int *maxit,
               double *reltol,
               int *order,
               double *lev,
               double *mix,
               double *orient,
               double *sigma2){

  int nv = *r, ngrad = *n, mc = *maxcomp;
  int iv, mc0, lpar, ord, maxc;
  double rss, krit, ttt, sw = 0.;
  double *siiq = (double *) R_alloc(ngrad, sizeof(double));
  int fail;                  // failure code for optim: zero is OK
  int fncount, grcount=0;               // number of calls to obj fct in optim
  double Fmin = 0.;          // minimal value of obj fct in optim
  int i, k;
  int gradind;
  double angles[2];
  double dir[3];
  double par[9];
  double cpar[9];
  double x[9];
  double tmp;
  int *mmask = (int *) R_alloc(9, sizeof(int));
  for (i = 0; i < 9; i++) mmask[i] = 1;
  double *w = (double *) R_alloc(9, sizeof(double));
  double theta, phi;
  int *ind = (int *) R_alloc(4, sizeof(int));
  double *z = (double *) R_alloc(ngrad*4, sizeof(double));
  double *qv = (double *) R_alloc(4*ngrad, sizeof(double));
  double *dqv = (double *) R_alloc(2*4*ngrad, sizeof(double));
  double *fv = (double *) R_alloc(ngrad, sizeof(double));
  double *dfv = (double *) R_alloc(9*ngrad, sizeof(double));
  double *work1 = (double *) R_alloc(ngrad, sizeof(double));
  optimexpl myoptimpar;
  myoptimpar.ngrad = ngrad;
  myoptimpar.grad = grad;
  myoptimpar.w = w;
  myoptimpar.z = z;
  myoptimpar.qv = qv;
  myoptimpar.dqv = dqv;
  myoptimpar.fv = fv;
  myoptimpar.dfv = dfv;
  myoptimpar.work1 = work1;

  for (iv = 0; iv < nv; iv++)
  {
    if (mask[iv] != 1)
    {
      order[iv] = 0;
      lev[get_ind2d(0, iv, 2)] = 0.;
      lev[get_ind2d(1, iv, 2)] = 0.;
      for (i = 0; i < mc; i++)
      {
        mix[get_ind2d(i, iv, mc)] = 0.;
        orient[get_ind3d(0, i, iv, 2, mc)] = 0.;
        orient[get_ind3d(1, i, iv, 2, mc)] = 0.;
      }
      sigma2[iv] = 0.;
    }
    else
    {
      // prepare signal
      for (i=0; i < ngrad; i++) siiq[i] = siq[get_ind2d(iv, i, nv)];

      // determine useful prime estimates
      mc0 = siind[get_ind2d(0, iv, mc+1)]; // order of model: possibly constant == mc since intial estimates guarantee mc0 == mc!!
      ord = mc0 + 1;
      gradind = siind[get_ind2d(1, iv, mc+1)];
      par[0] = -log(siq[get_ind2d(iv, gradind-1, nv)])*.2;
      for (i = 0; i < mc0; i++)
      {
        if (i>0) gradind = siind[get_ind2d(i+1, iv, mc+1)];
        dir[0] = grad[get_ind2d(gradind-1, 0, ngrad)];
        dir[1] = grad[get_ind2d(gradind-1, 1, ngrad)];
        dir[2] = grad[get_ind2d(gradind-1, 2, ngrad)];
        paroforient(dir, angles);    // INDEX????

        par[2*i + 1] = angles[0];
        par[2*i + 2] = angles[1];
      }

      // estimate models and select
      rss = R_PosInf;
      krit = R_PosInf;
      maxc = mc;
      for (k = mc0; k > 0; k--)
      {
        if (k < ord)
        {
          lpar = 2*k+1;
//          Rprintf("par in ");
          for (i = 0; i < lpar; i++) {cpar[i] = par[i]; /*Rprintf("%f ", cpar[i]);*/}
//          Rprintf("\n");
          myoptimpar.siq = siiq;
          switch (*method)
          { // R code guarantees method is 1, 2
            case 1:
              nmmin(lpar, cpar, x, &Fmin, fnpl,
                    &fail, R_NegInf, *reltol, &myoptimpar,
                    1.0, 0.5, 2.0, 0,
                    &fncount, *maxit);
              break;
//            case 2:
//              myoptimpar.p = *ep; // exp for Jian model
//              nmmin(lpar, cpar, x, &Fmin, fnpl,
//                    &fail, R_NegInf, *reltol, &myoptimpar,
//                    1.0, 0.5, 2.0, 0,
//                    &fncount, *maxit);
//              break;
            case 3:
              vmmin(lpar, cpar, &Fmin,
                     fnpl, grpl, *maxit, 0,
                     mmask, R_NegInf, *reltol, 100,
                     &myoptimpar, &fncount, &grcount, &fail);
              for (i = 0; i < lpar; i++) x[i] = cpar[i];
              break;
          }
//          Rprintf("par out ");
//          for (i = 0; i < lpar; i++) {Rprintf("%f ", x[i]);}
//          Rprintf("\n");

          if (Fmin*Fmin < rss) rss = Fmin*Fmin;
          Fmin = fnpl(lpar, x, &myoptimpar); // should return the same value
          ord = 0;
          sw = 0.;
          for (i = 0; i < k; i++)
          {
            ind[i] = i;
            if (myoptimpar.w[i] > 0.)
            {
              sw += myoptimpar.w[i];
              ord++;
            }
          }
          tmp = -log(sw);
          if ((x[0] < 0.) | (tmp < 0.))
          {
            //   parameters not interpretable reduce order
            ttt = R_PosInf;
            rss = R_PosInf;
            maxc = k-1;
          }
          else
          {
            ttt = Fmin*Fmin + (6.*ord+2.)/(ngrad - 3.*maxc - 1.) * rss;
            revsort(myoptimpar.w, ind, k);
//            Rprintf("weights ");
//            for (i = 0; i < ord; i++) Rprintf("%f ", myoptimpar.w[i]/sw);
//            Rprintf("\n");
            par[0] = x[0];
            for (i = 0; i < ord; i++)
            {
              par[2*i+1] = x[2*ind[i]+1];
              par[2*i+2] = x[2*ind[i]+2];
            }
          }

          if (ttt < krit)
          {
            krit = ttt;
            order[iv] = ord;
            lev[get_ind2d(0, iv, 2)] = x[0];
            lev[get_ind2d(1, iv, 2)] = tmp;
            if (ord == 1)
            {
              mix[get_ind2d(0, iv, mc)] = 1.;
              for (i = 1; i < mc; i++) mix[get_ind2d(i, iv, mc)] = 0.;
            }
            else
            {
              for (i = 0; i < ord; i++) mix[get_ind2d(i, iv, mc)] = myoptimpar.w[i]/sw;
              for (i = ord; i < mc; i++) mix[get_ind2d(i, iv, mc)] = 0.;
            }

            for (i = 0; i < k; i++)
            {
              theta = x[2*ind[i]+1];
              while (theta < 0.) 
              {
                theta += M_PI;
              }
              while (theta > M_PI) theta -= M_PI;
              phi = x[2*ind[i]+2];
              while (phi < 0.) phi += M_2PI;
              while (phi > M_2PI) phi -= M_2PI;
              orient[get_ind3d(0, i, iv, 2, mc)] = theta;
              orient[get_ind3d(1, i, iv, 2, mc)] = phi;
            }
          }
        }
      }
      sigma2[iv] = rss/(ngrad-3.*maxc-1.);
    }
    R_CheckUserInterrupt();
  }
  return;
}
