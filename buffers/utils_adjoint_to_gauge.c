#include "utils.ih"

#include <complex.h>

void adjoint_to_gauge(gauge_field_t out, adjoint_field_t const in)
{
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
      out[x][mu].c00 =  (0.5773502691896258 * in[x][mu].d8 + in[x][mu].d3) * I; 
      out[x][mu].c01 =  in[x][mu].d2                       +  in[x][mu].d1 * I; 
      out[x][mu].c02 =  in[x][mu].d5                       +  in[x][mu].d4 * I; 
      out[x][mu].c10 = -in[x][mu].d2                       +  in[x][mu].d1 * I; 
      out[x][mu].c11 =  (0.5773502691896258 * in[x][mu].d8 - in[x][mu].d3) * I; 
      out[x][mu].c12 =  in[x][mu].d7                       +  in[x][mu].d6 * I; 
      out[x][mu].c20 = -in[x][mu].d5                       +  in[x][mu].d4 * I; 
      out[x][mu].c21 = -in[x][mu].d7                       +  in[x][mu].d6 * I; 
      out[x][mu].c22 =  (1.154700538379252 * in[x][mu].d8) * I; 
    }
}