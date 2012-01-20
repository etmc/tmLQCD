/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
/* 
 * QUICKSORT 
 *
 * Sorts a double array using a non-recursive quicksort algorithm in
 * ascending order carrying along an int array.
 *  
 */

void quicksort(int n, double arr[], int idx[]){
  double v, td;
  int i, j, l, r, ti, tos, stack[32];
  
  l = 0; r = n-1; tos = -1;
  for (;;){
      while (r > l){ 
	  v = arr[r]; i = l; j = r-1;
	  for (;;){ 
	    while (arr[i] < v) i ++;
	    /* j > l prevents underflow */
	    while (arr[j] >= v && j > l) j --;
	    if (i >= j) break;
	    td = arr[i]; arr[i] = arr[j]; arr[j] = td;
	    ti = idx[i]; idx[i] = idx[j]; idx[j] = ti;
	  }
	  td = arr[i]; arr[i] = arr[r]; arr[r] = td;
	  ti = idx[i]; idx[i] = idx[r]; idx[r] = ti;
	  if (i-l > r-i){ 
	    stack[++tos] = l; stack[++tos] = i-1; l = i+1; 
	  }
	  else{	    
	    stack[++tos] = i+1; stack[++tos] = r; r = i-1; 
	  }
	  if(tos > 31) {
	    fprintf(stderr,"Error in quicksort! Aborting...!");fflush(stderr);
	    exit(31);
	  }
	} 
      if (tos == -1) break;
      r = stack[tos--]; l = stack[tos--]; 
    }
}
