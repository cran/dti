#include <R.h>
#include <Rmath.h>

double F77_SUB(besseli)(double *x, double *nu, double *expo) {
   return bessel_i(*x,*nu,*expo);
}

double F77_SUB(gammaf)(double *x) {
   return gammafn(*x);
}
double F77_SUB(digammaf)(double *x) {
   return digamma(*x);
}
double F77_SUB(lgammaf)(double *x) {
   return lgammafn(*x);
}
