#include "utils.ih"

void  print_su3(su3 const *in)
{
  printf("[ %c %14.12f %c %14.12f * 1i, %c %14.12f %c %14.12f * 1i, %c %14.12f %c %14.12f * 1i;\n",
          (creal(in->c00) > 0 ? '+' : '-'), fabs(creal(in->c00)), (cimag(in->c00) > 0 ? '+' : '-'), fabs(cimag(in->c00)), 
          (creal(in->c01) > 0 ? '+' : '-'), fabs(creal(in->c01)), (cimag(in->c01) > 0 ? '+' : '-'), fabs(cimag(in->c01)), 
          (creal(in->c02) > 0 ? '+' : '-'), fabs(creal(in->c02)), (cimag(in->c02) > 0 ? '+' : '-'), fabs(cimag(in->c02))); 
  printf("  %c %14.12f %c %14.12f * 1i, %c %14.12f %c %14.12f * 1i, %c %14.12f %c %14.12f * 1i;\n",
          (creal(in->c10) > 0 ? '+' : '-'), fabs(creal(in->c10)), (cimag(in->c10) > 0 ? '+' : '-'), fabs(cimag(in->c10)), 
          (creal(in->c11) > 0 ? '+' : '-'), fabs(creal(in->c11)), (cimag(in->c11) > 0 ? '+' : '-'), fabs(cimag(in->c11)), 
          (creal(in->c12) > 0 ? '+' : '-'), fabs(creal(in->c12)), (cimag(in->c12) > 0 ? '+' : '-'), fabs(cimag(in->c12))); 
  printf("  %c %14.12f %c %14.12f * 1i, %c %14.12f %c %14.12f * 1i, %c %14.12f %c %14.12f * 1i];\n",
          (creal(in->c20) > 0 ? '+' : '-'), fabs(creal(in->c20)), (cimag(in->c20) > 0 ? '+' : '-'), fabs(cimag(in->c20)), 
          (creal(in->c21) > 0 ? '+' : '-'), fabs(creal(in->c21)), (cimag(in->c21) > 0 ? '+' : '-'), fabs(cimag(in->c21)), 
          (creal(in->c22) > 0 ? '+' : '-'), fabs(creal(in->c22)), (cimag(in->c22) > 0 ? '+' : '-'), fabs(cimag(in->c22))); 
}