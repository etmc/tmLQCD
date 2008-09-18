#ifndef _DFL_PROJECTOR_H
#define _DFL_PROJECTOR_H

#include "su3spinor.h"

void project(spinor * const out, spinor * const in);
void project_left(spinor * const out, spinor * const in);
void project_right(spinor * const out, spinor * const in);
void project_left_D(spinor * const out, spinor * const in);
void D_project_right(spinor * const out, spinor * const in);
int check_projectors();
void check_little_D_inversion();
void check_local_D();
void free_dfl_projector();

void little_project(complex * const out, complex * const in, const int  N);
void little_P_L_D(complex * const out, complex * const in);
void little_D_P_R(complex * const out, complex * const in);
void little_P_R(complex * const out, complex * const in);
void little_P_L(complex * const out, complex * const in);

extern double dfl_little_D_prec;
extern int dfl_sloppy_prec;


#endif
