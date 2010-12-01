
















/////////////////////////////////////////////
// geometry- and nearest-neighbour indices //
/////////////////////////////////////////////


#ifdef MPI

// builds an array  iseven[global position]  to check wether is even or odd

void init_iseven() {

  int x0, x1, x2, x3;
  int ix;
  
  for (x0 = -1; x0 < T+1; x0++) {
    for (x1 = 0; x1 < LX; x1++) {
      for (x2 = 0; x2 < LY; x2++) {
        for (x3 = 0; x3 < LZ; x3++) {
          
          ix = Index(x0, x1, x2, x3);
          
	  if ((x0 + x1 + x2 + x3 + g_proc_coords[0]*T  + g_proc_coords[1]*LX + 
		                   g_proc_coords[2]*LY + g_proc_coords[3]*LZ) % 2 == 0) {
	    iseven[ix] = 1;
	  } 
	  else {
	    iseven[ix] = 0; 
	  }
        
        }}}}
     
}


// initialize nearest-neighbour table for gpu with even-odd enabled

void init_nnspinor_eo_mpi() {
									
  int x, y, z, t, ind, nnpos, j;					// mixed_solve_eo(...) allocates 8 integers per even or odd lattice site: size_t nnsize = 8*VOLUME*sizeof(int);
  									
  for (t = 0; t < T; t++) {						// loop goes over all INTERN latice sites !!
    for (x = 0; x < LX; x++) {						// doesn't refer to any EXTERN BOUNDARIES !!  ->  CORRESPONDS TO THE WHOLE LATTICE (I.E. WHEN NO SUBLATTICES ARE ASSIGNED) !!
      for (y = 0; y < LY; y++) {					//						  because of the behaviour of  g_iup[][] in the non-parallel case
        for (z = 0; z < LZ; z++) {
        								// NOTICE: g_ipt, g_iup, g_idn, and g_lexic2eosub  refer to pos. of lin. proj. pos. of the lattice
          ind = g_ipt[t][x][y][z];					// g_ipt[t][x][y][z] 	returns the linearly projected position of (t,x,y,z) of the lattice
          								//	indexes computed in geometry_eo() from geometry_eo.c
          								//	memory for the index array allocated by init_geometry_indices() from init_geometry_indices.c
          if ((t+x+y+z)%2 == 0) { // EVEN
          //if ((t + x + y + z + g_proc_coords[0]*T  + g_proc_coords[1]*LX + 
	  //	               g_proc_coords[2]*LY + g_proc_coords[3]*LZ) % 2 == 0) {
          
            nnpos = g_lexic2eosub[ind];					// g_lexic2eosub[ind] 	returns the position of [ind] in the sub-eo-notation
            										      ////////////////
            for (j = 0; j < 4; j++) { // plus direction			// here are also the  // BOUNDARIES //  included and properly mapped:
            								//		      ////////////////
              nn_eo[8*nnpos+j] = g_lexic2eosub[ g_iup[ind][j] ];	// g_iup[ind][j] 	returns the position of the nearest neighbour of [ind] in direction +[j]
            }								//				-->  for the non-parallized code g_iup[][] maps INTERN !!
            for (j = 0; j < 4; j++) { // minus direction
              nn_eo[8*nnpos+4+j] = g_lexic2eosub[ g_idn[ind][j] ];	// g_idn[ind][j] 	returns the position of the nearest neighbour of [ind] in direction -[j]
            }
          }
          
          else {		  // ODD
          
            nnpos = g_lexic2eosub[ind];
            
            for (j = 0; j < 4; j++) { // plus direction
              nn_oe[8*nnpos+j] = g_lexic2eosub[ g_iup[ind][j] ];	// nn_oe	      will return the nearest neigbours
            }								// nn_eo  and  nn_oe  strictly refer to the 4d-spacetime lattice
            
            for (j = 0; j < 4; j++) { // minus direction
              nn_oe[8*nnpos+4+j] = g_lexic2eosub[ g_idn[ind][j] ];
            }
          }
  }}}} // for loops
}




// the following functions can all be used to properly initialize the fields  eoidx_even[]  and  eoidx_odd[]  for addressing the gauge fields:


void init_idxgauge_mpi() {		// works!

  int t, x, y, z;
  int pos_eo, pos_global;
  
  for (t = -1; t < T+1; t++) {
    for (x = 0; x < LX; x++) {
      for (y = 0; y < LY; y++) {
        for (z = 0; z < LZ; z++) {
        
        //pos_global = g_ipt[t][x][y][z];
        pos_global = Index(t,x,y,z);
        pos_eo     = g_lexic2eosub[pos_global];
        
        //if ((t+x+y+z)%2 == 0) { // EVEN
        if ((t + x + y + z + g_proc_coords[0]*T  + g_proc_coords[1]*LX + 
	  	             g_proc_coords[2]*LY + g_proc_coords[3]*LZ) % 2 == 0) {
          eoidx_even[pos_eo] = g_eo2lexic[pos_eo];
        }
        else  {			// ODD
          eoidx_odd[pos_eo] = g_eo2lexic[(VOLUME+RAND)/2+pos_eo];
        }
  }}}} // for loop over the INTERN lattice
  
  //printf("This was init_idxgauge_mpi().\n");
  
}



/*
void init_idxgauge_mpi() {		// works!

  int t, x, y, z;
  int pos_eo, pos_global;
  
  for (t = -1; t < T+1; t++) {
    for (x = 0; x < LX; x++) {
      for (y = 0; y < LY; y++) {
        for (z = 0; z < LZ; z++) {
        
        //pos_global = g_ipt[t][x][y][z];
        pos_global = Index(t,x,y,z);
        pos_eo     = g_lexic2eosub[pos_global];
        
        //if ((t+x+y+z)%2 == 0) { // EVEN
        if ((t + x + y + z + g_proc_coords[0]*T  + g_proc_coords[1]*LX + 
	  	             g_proc_coords[2]*LY + g_proc_coords[3]*LZ) % 2 == 0) {
          eoidx_even[pos_eo] = pos_global;
        }
        else  {			// ODD
          eoidx_odd[pos_eo] = pos_global;
        }
  }}}} // for loop over the INTERN lattice
  
  //printf("This was init_idxgauge_mpi().\n");
  
}
*/


/*
void init_idxgauge_mpi() {		// works!

  int pos_eo, pos_global_even, pos_global_odd;
  
  for (pos_eo = 0; pos_eo < (VOLUME+RAND)/2; pos_eo++) {
      // even
      pos_global_even = g_eo2lexic[pos_eo];
      eoidx_even[pos_eo] = pos_global_even;
      // odd
      pos_global_odd = g_eo2lexic[(VOLUME+RAND)/2 + pos_eo];
      eoidx_odd[pos_eo] = pos_global_odd;
  }
  
  //printf("This was init_idxgauge_mpi().\n");
  
}
*/


/*
void init_idxgauge_mpi() {		// works!

  int pos_eo, pos_global;
  
  for (pos_global = 0; pos_global < (VOLUME+RAND); pos_global++) {
  
    pos_eo = g_lexic2eosub[pos_global];
    
    if (iseven[pos_global] == 1) {
    //if (pos_global%2 == 0) {
      eoidx_even[pos_eo] = pos_global;
    }
    else {
      eoidx_odd[pos_eo]  = pos_global;
    }
      
  }
  
  //printf("This was init_idxgauge_mpi().\n");
  
}
*/


/*
void init_idxgauge_mpi() {		// works!

  int x, y, z, t;
  int ind;
  int evenpos = 0;
  int oddpos = 0;
  
  for (t = 0; t < T; t++) {
    for (x = 0; x < LX; x++) {
      for (y = 0; y < LY; y++) {
        for (z = 0; z < LZ; z++) {
          ind = g_ipt[t][x][y][z];
          if ((t+x+y+z) % 2 == 0) {
            eoidx_even[evenpos] = ind;
            evenpos++;
          }
          else {
            eoidx_odd[oddpos] = ind;
            oddpos++;
          }
  }}}} // INTERN
  
  
  		t = T;
  		  for (x = 0; x < LX; x++) {
  		    for (y = 0; y < LY; y++) {
  		      for (z = 0; z < LZ; z++) {
  		        ind = VOLUME + z + LZ*y + LZ*LY*x;
  		        //if (iseven[ind] == 1) {
  		        if ((t+x+y+z) % 2 == 0) {
  		          eoidx_even[evenpos] = ind;
  		          evenpos++;
  		        }
  		        else {
  		          eoidx_odd[oddpos] = ind;
  		          oddpos++;
  		        }
  		}}} // EXTERN
  
  
  				t = -1;
  				  for (x = 0; x < LX; x++) {
  				    for (y = 0; y < LY; y++) {
  				      for (z = 0; z < LZ; z++) {
  				        ind = VOLUME + LX*LY*LZ + z + LZ*y + LZ*LY*x;
  				        //if (iseven[ind] == 1) {
  				        if ((t+x+y+z) % 2 == 0) {
  				          eoidx_even[evenpos] = ind;
  				          evenpos++;
  				        }
  				        else {
  				          eoidx_odd[oddpos] = ind;
  				          oddpos++;
  				        }
  				}}} // EXTERN
  				
  //printf("This was init_idxgauge_mpi().\n");
  
}
*/


#endif	// MPI






/////////////////////
// initializations //
/////////////////////


#ifdef MPI

// allocates memory for the fields for the alternative way of addressing positions in dev_Hopping_Matrix_alternate()

void init_gpu_indexfields() {
  
  size_t size;
  
  		// debug
  		/*
  		printf("Test: %p ?= %p\n", *g_iup, g_iup[0]);
  		printf("Test: %p ?= %p\n", *g_idn, g_idn[0]);
  		printf("Test: %p ?= %p\n", g_lexic2eo, &g_lexic2eo[0]);
  		printf("Test: %p ?= %p\n", g_lexic2eosub, &g_lexic2eosub[0]);
  		printf("Test: %p ?= %p\n", g_eo2lexic, &g_eo2lexic[0]);
  		printf("Test: %p ?= %p ?= %p\n", ***g_ipt, &g_ipt[0][0][0][0], g_ipt[0][0][0]);
  		*/
  
  size = 4*(VOLUME+RAND)*sizeof(int);
  cudaMalloc((void **) &dev_g_iup, size);
  cudaMalloc((void **) &dev_g_idn, size);
  cudaMemcpy(dev_g_iup, g_iup[0], size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_g_idn, g_idn[0], size, cudaMemcpyHostToDevice);
  
  size = (VOLUME+RAND)*sizeof(int);
  cudaMalloc((void **) &dev_g_lexic2eo, size);
  cudaMalloc((void **) &dev_g_lexic2eosub, size);
  cudaMemcpy(dev_g_lexic2eo, g_lexic2eo, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_g_lexic2eosub, g_lexic2eosub, size, cudaMemcpyHostToDevice);
  
  size = (VOLUME+RAND)*sizeof(int);
  cudaMalloc((void **) &dev_g_eo2lexic, size);
  cudaMemcpy(dev_g_eo2lexic, g_eo2lexic, size, cudaMemcpyHostToDevice);
  
  size = VOLUME*sizeof(int);
  cudaMalloc((void **) &dev_g_ipt, size);
  cudaMemcpy(dev_g_ipt, g_ipt[0][0][0], size, cudaMemcpyHostToDevice);
  
}




// frees the memory

void free_gpu_indexfields() {

  cudaFree(dev_g_iup);
  cudaFree(dev_g_idn);
  
  cudaFree(dev_g_lexic2eo);
  cudaFree(dev_g_lexic2eosub);
  
  cudaFree(dev_g_eo2lexic);
  
  cudaFree(dev_g_ipt);

}

#endif	// MPI








