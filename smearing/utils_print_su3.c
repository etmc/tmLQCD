#include "utils.ih"

void  print_su3(su3 *in)
{
  printf("[ %12.14f + %12.14f * i, %12.14f + %12.14f * i, %12.14f + %12.14f * i; \n",
      creal(in->c00), cimag(in->c00), 
      creal(in->c01), cimag(in->c01), 
      creal(in->c02), cimag(in->c02)) ; 
  printf("  %12.14f + %12.14f * i, %12.14f + %12.14f * i, %12.14f + %12.14f * i; \n",
      creal(in->c10), cimag(in->c10), 
      creal(in->c11), cimag(in->c11), 
      creal(in->c12), cimag(in->c12)) ; 
  printf("  %12.14f + %12.14f * i, %12.14f + %12.14f * i, %12.14f + %12.14f * i  ] \n",
      creal(in->c20), cimag(in->c20), 
      creal(in->c21), cimag(in->c21), 
      creal(in->c22), cimag(in->c22)) ; 
}