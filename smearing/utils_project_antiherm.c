#include "utils.ih"

void project_antiherm(su3 *omega)
{
  static const double fac_3 = 1.00 / 3.00;
  double tr_omega = (omega->c00.im + omega->c11.im + omega->c22.im) * fac_3; 

  omega->c00.re = 0.00;
  omega->c00.im -= tr_omega;

  omega->c11.re = 0.00;
  omega->c11.im -= tr_omega;

  omega->c22.re = 0.00;
  omega->c22.im -= tr_omega;

  omega->c01.re -= omega->c10.re;  omega->c01.re *= 0.50; 
  omega->c01.im += omega->c10.im;  omega->c01.im *= 0.50; 
  omega->c10.re  = -omega->c01.re; 
  omega->c10.im  =  omega->c01.im; 

  omega->c02.re -= omega->c20.re;  omega->c02.re *= 0.50; 
  omega->c02.im += omega->c20.im;  omega->c02.im *= 0.50;
  omega->c20.re  = -omega->c02.re; 
  omega->c20.im  =  omega->c02.im; 

  omega->c12.re -= omega->c21.re;  omega->c12.re *= 0.50; 
  omega->c12.im += omega->c21.im;  omega->c12.im *= 0.50; 
  omega->c21.re  = -omega->c12.re; 
  omega->c21.im  =  omega->c12.im; 
}
