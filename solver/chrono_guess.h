/* $Id$ */

#ifndef _CHRONO_GUESS_H
#define _CHRONO_GUESS_H

#include "solver/matrix_mult_typedef.h"

int chrono_guess(spinor * const trial, spinor * const phi, spinor ** const v, int index_array[], 
		 const int _N, int * _n, const int V, matrix_mult f);

#endif
