#include "utils.ih"

void  print_su3(su3 *in)
{
  printf("[ %12.14f + %12.14f * i, %12.14f + %12.14f * i, %12.14f + %12.14f * i; \n",
      in->c00.re,in->c00.im, 
      in->c01.re,in->c01.im, 
      in->c02.re,in->c02.im ) ; 
  printf("  %12.14f + %12.14f * i, %12.14f + %12.14f * i, %12.14f + %12.14f * i; \n",
      in->c10.re,in->c10.im, 
      in->c11.re,in->c11.im, 
      in->c12.re,in->c12.im ) ; 
  printf("  %12.14f + %12.14f * i, %12.14f + %12.14f * i, %12.14f + %12.14f * i  ] \n",
      in->c20.re,in->c20.im, 
      in->c21.re,in->c21.im, 
      in->c22.re,in->c22.im ) ; 
}