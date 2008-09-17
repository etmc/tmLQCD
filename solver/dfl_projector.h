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

extern double dfl_little_D_prec;
extern int dfl_sloppy_prec;


#endif
