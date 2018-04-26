/***********************************************************************
 *
 * Copyright (C) 2018 Bartosz Kostrzewa 
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
 ************************************************************************/

#ifndef TENSORS_H
#define TENSORS_H

typedef struct epsilon4_t {
  int N;
  double eps_val[24];
  int eps_idx[24][4];
} epsilon4_t; 

typedef struct epsilon3_t {
  int N;
  double eps_val[6];
  int eps_idx[6][3];
} epsilon3_t;

static inline epsilon3_t new_epsilon3(void) {
  epsilon3_t ret;

  ret.N = 6;

  int i = 0;
  int p = 0;
  for( int i1 = 1; i1 <= 3; i1++ ){
    for( int i2 = 1; i2 <= 3; i2++ ){
      for( int i3 = 1; i3 <= 3; i3++ ){
        // for eps_123 we have: (1 - 2)(1 - 3)(2 - 3) = -2 
        //                      -> minus sign
        p = -(i1 - i2)*(i1 - i3)*(i2 - i3);
        if( p != 0 ){
          ret.eps_val[i] = p > 0 ? 1 : -1;
          ret.eps_idx[i][0] = i1-1;
          ret.eps_idx[i][1] = i2-1;
          ret.eps_idx[i][2] = i3-1;
          i++;
        }
      } 
    }
  }
  return(ret);
}

// note that this is the Euclidean eps_ijkl, for which we have eps_1234 = 1,
// whereas in Minkowski space we have eps_0123 = -1
static inline epsilon4_t new_epsilon4(void) {
  epsilon4_t ret;

  ret.N = 24;

  int i = 0;
  int p = 0;
  for( int i1 = 1; i1 <= 4; i1++ ){
    for( int i2 = 1; i2 <= 4; i2++ ){
      for( int i3 = 1; i3 <= 4; i3++ ){
        for( int i4 = 1; i4 <= 4; i4++ ){
          // for eps_1234 we have: (1 - 2)(1 - 3)(1 - 4)(2 - 3)(2 - 4)(3 - 4) = 12 
          //                       -> NO minus sign
          p = (i1 - i2)*(i1 - i3)*(i1 - i4)*(i2 - i3)*(i2 - i4)*(i3 - i4);
          if( p != 0 ){
            ret.eps_val[i] = p > 0 ? 1 : -1;
            ret.eps_idx[i][0] = i1-1;
            ret.eps_idx[i][1] = i2-1;
            ret.eps_idx[i][2] = i3-1;
            ret.eps_idx[i][3] = i4-1;
            i++;
          }
        }
      } 
    }
  }
  return(ret);
}

#endif

