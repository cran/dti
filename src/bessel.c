#include <R.h>
#include <Rmath.h>

double F77_SUB(besseli)(double *x, double *nu, double *expo) {
   return bessel_i(*x,*nu,*expo);
}

