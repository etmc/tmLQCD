/* $Id$ */

#include<stdlib.h>
#include<stdio.h>
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
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
