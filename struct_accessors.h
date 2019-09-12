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
      if(reim==0) return( ((double*)(&matrix->c00))[1] );
      else return( ((double*)(&matrix->c00))[2] );
      break;
    case 1:
      if(reim==0) return( ((double*)(&matrix->c01))[1] );
      else return( ((double*)(&matrix->c01))[2] );
      break;
    case 2:
      if(reim==0) return( ((double*)(&matrix->c02))[1] );
      else return( ((double*)(&matrix->c02))[2] );
      break;
    case 3:
      if(reim==0) return( ((double*)(&matrix->c10))[1] );
      else return( ((double*)(&matrix->c10))[2] );
      break;
    case 4:
      if(reim==0) return( ((double*)(&matrix->c11))[1] );
      else return( ((double*)(&matrix->c11))[2] );
      break;
    case 5:
      if(reim==0) return( ((double*)(&matrix->c12))[1] );
      else return( ((double*)(&matrix->c12))[2] );
      break;
    case 6:
      if(reim==0) return( ((double*)(&matrix->c20))[1] );
      else return( ((double*)(&matrix->c20))[2] );
      break;
    case 7:
      if(reim==0) return( ((double*)(&matrix->c21))[1] );
      else return( ((double*)(&matrix->c21))[2] );
      break;
    case 8:
      if(reim==0) return( ((double*)(&matrix->c22))[1] );
      else return( ((double*)(&matrix->c22))[2] );
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
      if(reim==0) return( ((double*)(&matrix->s0.c0))[1] );
      else return( ((double*)(&matrix->s0.c0))[2] );
      break;
    case 1:
      if(reim==0) return( ((double*)(&matrix->s0.c1))[1] );
      else return( ((double*)(&matrix->s0.c1))[2] );
      break;
    case 2:
      if(reim==0) return( ((double*)(&matrix->s0.c2))[1] );
      else return( ((double*)(&matrix->s0.c2))[2] );
      break;
    case 3:
      if(reim==0) return( ((double*)(&matrix->s1.c0))[1] );
      else return( ((double*)(&matrix->s1.c0))[2] );
      break;
    case 4:
      if(reim==0) return( ((double*)(&matrix->s1.c1))[1] );
      else return( ((double*)(&matrix->s1.c1))[2] );
      break;
    case 5:
      if(reim==0) return( ((double*)(&matrix->s1.c2))[1] );
      else return( ((double*)(&matrix->s1.c2))[2] );
      break;
    case 6:
      if(reim==0) return( ((double*)(&matrix->s2.c0))[1] );
      else return( ((double*)(&matrix->s2.c0))[2] );
      break;
    case 7:
      if(reim==0) return( ((double*)(&matrix->s2.c1))[1] );
      else return( ((double*)(&matrix->s2.c1))[2] );
      break;
    case 8:
      if(reim==0) return( ((double*)(&matrix->s2.c2))[1] );
      else return( ((double*)(&matrix->s2.c2))[2] );
      break;
    case 9:
      if(reim==0) return( ((double*)(&matrix->s3.c0))[1] );
      else return( ((double*)(&matrix->s3.c0))[2] );
      break;
    case 10:
      if(reim==0) return( ((double*)(&matrix->s3.c1))[1] );
      else return( ((double*)(&matrix->s3.c1))[2] );
      break;
    case 11:
      if(reim==0) return( ((double*)(&matrix->s3.c2))[1] );
      else return( ((double*)(&matrix->s3.c2))[2] );
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
      *(((double*const)(&(matrix->s0.c0)))  ) = rein;
      *(((double*const)(&(matrix->s0.c0)))+1) = imin;
      break;
    case 1:
      *(((double*const)(&(matrix->s0.c1)))  ) = rein;
      *(((double*const)(&(matrix->s0.c1)))+1) = imin;
      break;
    case 2:
      *(((double*const)(&(matrix->s0.c2)))  ) = rein;
      *(((double*const)(&(matrix->s0.c2)))+1) = imin;
      break;
    case 3:
      *(((double*const)(&(matrix->s1.c0)))  ) = rein;
      *(((double*const)(&(matrix->s1.c0)))+1) = imin;
      break;
    case 4:
      *(((double*const)(&(matrix->s1.c1)))  ) = rein;
      *(((double*const)(&(matrix->s1.c1)))+1) = imin;
      break;
    case 5:
      *(((double*const)(&(matrix->s1.c2)))  ) = rein;
      *(((double*const)(&(matrix->s1.c2)))+1) = imin;
      break;
    case 6:
      *(((double*const)(&(matrix->s2.c0)))  ) = rein;
      *(((double*const)(&(matrix->s2.c0)))+1) = imin;
      break;
    case 7:
      *(((double*const)(&(matrix->s2.c1)))  ) = rein;
      *(((double*const)(&(matrix->s2.c1)))+1) = imin;
      break;
    case 8:
      *(((double*const)(&(matrix->s2.c2)))  ) = rein;
      *(((double*const)(&(matrix->s2.c2)))+1) = imin;
      break;
    case 9:
      *(((double*const)(&(matrix->s3.c0)))  ) = rein;
      *(((double*const)(&(matrix->s3.c0)))+1) = imin;
      break;
    case 10:
      *(((double*const)(&(matrix->s3.c1)))  ) = rein;
      *(((double*const)(&(matrix->s3.c1)))+1) = imin;
      break;
    case 11:
      *(((double*const)(&(matrix->s3.c2)))  ) = rein;
      *(((double*const)(&(matrix->s3.c2)))+1) = imin;
      break;
    default:
      exit(-224);
  }
}

static inline void spinor_set_elem(spinor* const matrix, int s, int c, const double rein, const double imin){
  spinor_set_elem_linear(matrix, 3*s+c, rein, imin);
}

#endif
