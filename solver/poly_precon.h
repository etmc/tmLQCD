#ifndef _POLY_PRECON_H
#define _POLY_PRECON_H

void poly_precon(spinor * const, spinor * const, const double prec, const int n);
void poly_nonherm_precon(spinor * const R, spinor * const S, 
			 const double e, const int n);

#endif
