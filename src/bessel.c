#include <R.h>
#include <Rmath.h>

double F77_SUB(besseli)(double *x, double *nu, double *expo) {
   return bessel_i(*x,*nu,*expo);
}

double F77_SUB(bessliex)(double *x, double *nu, double *expo, double *bi) {
   return bessel_i_ex(*x,*nu,*expo, bi);
}

double F77_SUB(gammaf)(double *x) {
   return gammafn(*x);
}
