#include "utils.ih"

#include <complex.h>

void adjoint_to_gauge(gauge_field_t out, adjoint_field_t const in)
{
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
      out.field[x][mu].c00 =  (0.5773502691896258 * in.field[x][mu].d8 + in.field[x][mu].d3) * I; 
      out.field[x][mu].c01 =  in.field[x][mu].d2 +  in.field[x][mu].d1 * I; 
      out.field[x][mu].c02 =  in.field[x][mu].d5 +  in.field[x][mu].d4 * I; 
      out.field[x][mu].c10 = -in.field[x][mu].d2 +  in.field[x][mu].d1 * I; 
      out.field[x][mu].c11 =  (0.5773502691896258 * in.field[x][mu].d8 - in.field[x][mu].d3) * I; 
      out.field[x][mu].c12 =  in.field[x][mu].d7 +  in.field[x][mu].d6 * I; 
      out.field[x][mu].c20 = -in.field[x][mu].d5 +  in.field[x][mu].d4 * I; 
      out.field[x][mu].c21 = -in.field[x][mu].d7 +  in.field[x][mu].d6 * I; 
      out.field[x][mu].c22 =  (1.154700538379252 * in.field[x][mu].d8) * I; 
    }
}