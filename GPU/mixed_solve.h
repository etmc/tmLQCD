#ifndef _MIXED_SOLVE_H_

void initnn();

extern "C" int mixed_solve  (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps,const int rel_prec, const int N);
extern "C" int mixed_solveD (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps,const int rel_prec, const int N);

extern "C" int mixed_solve_eo  (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec, const int N);
extern "C" int mixed_solve_eoD (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec, const int N);


extern "C" int bind_texture_spin(dev_spinor* s, int i);
extern "C" int unbind_texture_spin(int i);

extern "C" int bind_texture_nn(int* nn);
extern "C" int unbind_texture_nn();

#define _MIXED_SOLVE_H_
#endif
