/***********************************************************************
 *
 * Copyright (C) 2017 Bartosz Kostrzewa
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************************************/

#ifndef STRUCT_ACCESSORS_H
#define STRUCT_ACCESSORS_H

#include "su3.h"
#include <stdlib.h>

static inline double su3_get_elem_linear(const su3* const matrix, int cc, int reim){
  switch(cc){
    case 0:
      if(reim==0) return( creal(matrix->c00) );
      else return( cimag(matrix->c00) );
      break;
    case 1:
      if(reim==0) return( creal(matrix->c01) );
      else return( cimag(matrix->c01) );
      break;
    case 2:
      if(reim==0) return( creal(matrix->c02) );
      else return( cimag(matrix->c02) );
      break;
    case 3:
      if(reim==0) return( creal(matrix->c10) );
      else return( cimag(matrix->c10) );
      break;
    case 4:
      if(reim==0) return( creal(matrix->c11) );
      else return( cimag(matrix->c11) );
      break;
    case 5:
      if(reim==0) return( creal(matrix->c12) );
      else return( cimag(matrix->c12) );
      break;
    case 6:
      if(reim==0) return( creal(matrix->c20) );
      else return( cimag(matrix->c20) );
      break;
    case 7:
      if(reim==0) return( creal(matrix->c21) );
      else return( cimag(matrix->c21) );
      break;
    case 8:
      if(reim==0) return( creal(matrix->c22) );
      else return( cimag(matrix->c22) );
      break;
    default:
      exit(-222);
  }
}

static inline double su3_get_elem(const su3* const matrix, int c0, int c1, int reim){
  return su3_get_elem_linear(matrix, 3*c0+c1, reim);
}

static inline double spinor_get_elem_linear(const spinor* const matrix, int sc, int reim){
  switch(sc){
    case 0:
      if(reim==0) return( creal(matrix->s0.c0) );
      else return( cimag(matrix->s0.c0) );
      break;
    case 1:
      if(reim==0) return( creal(matrix->s0.c1) );
      else return( cimag(matrix->s0.c1) );
      break;
    case 2:
      if(reim==0) return( creal(matrix->s0.c2) );
      else return( cimag(matrix->s0.c2) );
      break;
    case 3:
      if(reim==0) return( creal(matrix->s1.c0) );
      else return( cimag(matrix->s1.c0) );
      break;
    case 4:
      if(reim==0) return( creal(matrix->s1.c1) );
      else return( cimag(matrix->s1.c1) );
      break;
    case 5:
      if(reim==0) return( creal(matrix->s1.c2) );
      else return( cimag(matrix->s1.c2) );
      break;
    case 6:
      if(reim==0) return( creal(matrix->s2.c0) );
      else return( cimag(matrix->s2.c0) );
      break;
    case 7:
      if(reim==0) return( creal(matrix->s2.c1) );
      else return( cimag(matrix->s2.c1) );
      break;
    case 8:
      if(reim==0) return( creal(matrix->s2.c2) );
      else return( cimag(matrix->s2.c2) );
      break;
    case 9:
      if(reim==0) return( creal(matrix->s3.c0) );
      else return( cimag(matrix->s3.c0) );
      break;
    case 10:
      if(reim==0) return( creal(matrix->s3.c1) );
      else return( cimag(matrix->s3.c1) );
      break;
    case 11:
      if(reim==0) return( creal(matrix->s3.c2) );
      else return( cimag(matrix->s3.c2) );
      break;
    default:
      exit(-223);
  }
}

static inline double spinor_get_elem(const spinor* const matrix, int s, int c, int reim){
  return spinor_get_elem_linear(matrix, 3*s+c, reim);
}

static inline void spinor_set_elem_linear(spinor* const matrix, int sc, const double rein, const double imin){
  switch(sc){
    case 0:
      *(reinterpret_cast<double*const>(&(matrix->s0.c0))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s0.c0))+1) = imin;
      break;
    case 1:
      *(reinterpret_cast<double*const>(&(matrix->s0.c1))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s0.c1))+1) = imin;
      break;
    case 2:
      *(reinterpret_cast<double*const>(&(matrix->s0.c2))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s0.c2))+1) = imin;
      break;
    case 3:
      *(reinterpret_cast<double*const>(&(matrix->s1.c0))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s1.c0))+1) = imin;
      break;
    case 4:
      *(reinterpret_cast<double*const>(&(matrix->s1.c1))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s1.c1))+1) = imin;
      break;
    case 5:
      *(reinterpret_cast<double*const>(&(matrix->s1.c2))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s1.c2))+1) = imin;
      break;
    case 6:
      *(reinterpret_cast<double*const>(&(matrix->s2.c0))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s2.c0))+1) = imin;
      break;
    case 7:
      *(reinterpret_cast<double*const>(&(matrix->s2.c1))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s2.c1))+1) = imin;
      break;
    case 8:
      *(reinterpret_cast<double*const>(&(matrix->s2.c2))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s2.c2))+1) = imin;
      break;
    case 9:
      *(reinterpret_cast<double*const>(&(matrix->s3.c0))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s3.c0))+1) = imin;
      break;
    case 10:
      *(reinterpret_cast<double*const>(&(matrix->s3.c1))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s3.c1))+1) = imin;
      break;
    case 11:
      *(reinterpret_cast<double*const>(&(matrix->s3.c2))  ) = rein;
      *(reinterpret_cast<double*const>(&(matrix->s3.c2))+1) = imin;
      break;
    default:
      exit(-224);
  }
}

static inline void spinor_set_elem(spinor* const matrix, int s, int c, const double rein, const double imin){
  spinor_set_elem_linear(matrix, 3*s+c, rein, imin);
}

#endif