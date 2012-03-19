#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "global.h"
#include "GPU/cudadefs.h"
#include "su3.h"
#include "geometry_eo.h"
#include "start.h"
#include "temporalgauge.h"
#include "stdio.h"
#include "stdlib.h"
#ifdef MPI
  #include<mpi.h>
  #include "mpi_init.h"
#endif



su3 * g_trafo;
su3 * tempgauge_field = NULL;

su3 * left;
su3 * right;



static su3 unit_su3 (void)
{
   su3 u;
   _su3_one(u);
   return(u);
}



/*copy a complete gauge field*/
/* THINK OF PARALLELIZATION (RAND!!!)*/
void copy_gauge_field (su3 ** to, su3 ** from)
{
  for (int ix = 0; ix < VOLUME; ix++)
  {
    _su3_assign(to[ix][0], from[ix][0]);
    _su3_assign(to[ix][1], from[ix][1]);
    _su3_assign(to[ix][2], from[ix][2]);
    _su3_assign(to[ix][3], from[ix][3]);
  }
}

/*
  Set the trafo field for a temporal gauge 
  g(t=0) == ID 
  other g's are determined recursively from U (gfield) requiering that U^{'}_0 != ID
  => only the U(t=T-1) are not ID!!
*/
int init_temporalgauge_trafo (const int V, su3** gfield) {

#ifndef MPI

   int it, iz, iy, ix;
   
   int pos;
   
   if ((void *)(g_trafo = (su3 *) calloc(V, sizeof(su3))) == NULL ) {
    printf("malloc error in 'init_temporalgauge_trafo'\n"); 
    return(2);
  }
  
  /* initialize first timeslice (t=0) with unit matrices*/
  for (ix = 0; ix < LX; ix++) {
    for (iy = 0; iy < LY; iy++) {
      for (iz = 0; iz < LZ; iz++) {
        g_trafo[g_ipt[0][ix][iy][iz]] = unit_su3();
      }
    }
  }
  
  /* U^{'}_0(x)  g(x) U_0(x) g^{+}(x+0) != ID   =>  g_(x+0) = g(x) U_0(x)  */
  for (it = 1; it < T; it++) {
    for (ix = 0; ix < LX; ix++) {
      for (iy = 0; iy < LY; iy++) {
        for (iz = 0; iz < LZ; iz++) {
          pos = g_ipt[it][ix][iy][iz];
          _su3_times_su3( g_trafo[ g_ipt[it  ][ix][iy][iz] ] ,
                          g_trafo[ g_ipt[it-1][ix][iy][iz] ] ,
                          //gfield [ g_ipt[it-1][ix][iy][iz] ] [0]  );
                          gfield [ g_idn[pos][0]           ] [0]  );
        }
      }
    } 
  }

#else // MPI

  int it, iz, iy, ix;
  
  int pos;
  
  MPI_Status status;
  
  
  
  if ((void *)(left = (su3 *) calloc(LX*LY*LZ, sizeof(su3))) == NULL ) {		// allocates memory for a time-slice of su3-matrices
    printf("malloc error in 'init_temporalgauge_trafo_mpi'\n"); 
    return(-1);
  }
  
  if ((void *)(right = (su3 *) calloc(LX*LY*LZ, sizeof(su3))) == NULL ) {		// allocates memory for a time-slice of su3-matrices
    printf("malloc error in 'init_temporalgauge_trafo_mpi'\n"); 
    return(-1);
  }
  
  
  
  
  if ((void *)(g_trafo = (su3 *) calloc(V, sizeof(su3))) == NULL ) {			// allocates memory for V su3-matrices
    printf("malloc error in 'init_temporalgauge_trafo'\n"); 
    return(2);
  } 
  
  
  
  
  //////////////////////////////////////////////
  // initializing the transformation matrices //
  //////////////////////////////////////////////
  
  
  // first process in t-direction
  
  if (g_cart_id == 0) {
  	
  	/* initialize first timeslice (t=0) with unit matrices*/
  	for (ix = 0; ix < LX; ix++) {
  	  for (iy = 0; iy < LY; iy++) {
  	    for (iz = 0; iz < LZ; iz++) {
  	      g_trafo[g_ipt[0][ix][iy][iz]] = unit_su3();					// g_trafo[0-th time slice]  =  ID
  	    }
  	  }
  	}
  	
  	/* U^{'}_0(x)  =  g(x) U_0(x) g^{+}(x+0)  !=  ID   =>   g_(x+0)  =  g(x) U_0(x)  */
  	for (it = 1; it < T; it++) {
  	  for (ix = 0; ix < LX; ix++) {
  	    for (iy = 0; iy < LY; iy++) {
  	      for (iz = 0; iz < LZ; iz++) {
  	        _su3_times_su3( g_trafo[ g_ipt[it  ][ix][iy][iz] ] , 				// g_trafo[next t-slice]  =  g_trafo[old t-slice]  *  gfield[old t-slice][t-dir.]
  	                        g_trafo[ g_ipt[it-1][ix][iy][iz] ] ,
  	                        gfield [ g_ipt[it-1][ix][iy][iz] ] [0] );
  	        
  	      }
  	    }
  	  } 
  	}
  	
  	
  	// sending
  	MPI_Send((void *)(g_trafo+(T-1)*LX*LY*LZ), LX*LY*LZ, mpi_su3, g_nb_t_up, 0, g_cart_grid);
  	//MPI_Send((void *)(g_trafo+(T-1)*LX*LY*LZ), LX*LY*LZ, mpi_su3, g_cart_id+1, 0, g_cart_grid);
  	
  	printf("g_cart_id = %i has send a message to %i\n", g_cart_id, g_nb_t_up);
  	
  	
  } // first process
  
  
  
  
  // following processes
  
  else {
  	
  	// receiving
  	MPI_Recv((void *)left, LX*LY*LZ, mpi_su3, g_nb_t_dn, 0, g_cart_grid, &status);
  	//MPI_Recv((void *)left, LX*LY*LZ, mpi_su3, g_cart_id-1, 0, g_cart_grid, &status);
  	
  	
  	printf("g_cart_id = %i has received a message from %i\n", g_cart_id, g_nb_t_dn);
  	
  	it = 0;
  	for (ix = 0; ix < LX; ix++) {
  	  for (iy = 0; iy < LY; iy++) {
  	    for (iz = 0; iz < LZ; iz++) {
  	      pos = g_ipt[it][ix][iy][iz];
  	      _su3_times_su3( g_trafo[ g_ipt[it  ][ix][iy][iz] ] ,				// g_trafo[0-th time slice]  =  left[xchanged t-slice]  * gfield[
  	                      left   [ g_ipt[it  ][ix][iy][iz] ] ,
  	                      gfield [ g_idn[pos ][0]          ] [0] );				// notice: have to access the RAND region of the gauge field
  	    }
  	  }
  	}
  	
  	
  	for (it = 1; it < T; it++) {
  	  for (ix = 0; ix < LX; ix++) {
  	    for (iy = 0; iy < LY; iy++) {
  	      for (iz = 0; iz < LZ; iz++) {
  	        _su3_times_su3( g_trafo[ g_ipt[it  ][ix][iy][iz] ] ,
  	                        g_trafo[ g_ipt[it-1][ix][iy][iz] ] ,
  	                        gfield [ g_ipt[it-1][ix][iy][iz] ] [0] );
  	        
  	      }
  	    }
  	  } 
  	}
  	
  	
  	// sending
  	if (g_cart_id != g_nproc-1) {
  	  MPI_Send((void *)(g_trafo+(T-1)*LX*LY*LZ), LX*LY*LZ, mpi_su3, g_nb_t_up, 0, g_cart_grid);
  	  //MPI_Send((void *)(g_trafo+(T-1)*LX*LY*LZ), LX*LY*LZ, mpi_su3, g_cart_id+1, 0, g_cart_grid);
  	  
  	  printf("g_cart_id = %i has send a message to %i\n", g_cart_id, g_nb_t_up);
  	  
  	}
  	
  
  } // following processes
  
  
  
  
  ////////////////////////////////////////////
  // exchanging the transformation matrices //
  ////////////////////////////////////////////
  
  
  MPI_Sendrecv((void *)(g_trafo), LX*LY*LZ, mpi_su3, g_nb_t_dn, 1,
               (void *)(right  ), LX*LY*LZ, mpi_su3, g_nb_t_up, 1,
               g_cart_grid, &status);
  
  printf("g_cart_id = %i has send to %i and received from %i\n", g_cart_id, g_nb_t_dn, g_nb_t_up);


#endif // MPI


  /* 
    allocate and initialize g_tempgauge_field which holds a copy of the 
    global gauge field g_gauge_field which is copied back after the inversion
    when the temporal gauge is undone again
  */
  
  int i = 0;
  
  if ((void *)(g_tempgauge_field = (su3 **) calloc(V, sizeof(su3*))) == NULL ) {
    printf ("malloc error in 'init_temporalgauge_trafo'\n"); 
    return(1);
  }
  if ((void *)(tempgauge_field = (su3 *) calloc(4*V+1, sizeof(su3))) == NULL ) {
    printf ("malloc error in 'init_temporalgauge_trafo'\n"); 
    return(2);
  }
  
  #if (defined SSE || defined SSE2 || defined SSE3)
    g_tempgauge_field[0] = (su3*)(((unsigned long int)(tempgauge_field)+ALIGN_BASE)&~ALIGN_BASE);
  #else
    g_tempgauge_field[0] = tempgauge_field;
  #endif
  
  for(i = 1; i < V; i++){
    g_tempgauge_field[i] = g_tempgauge_field[i-1]+4;
  }

  /* copy the original field */
  copy_gauge_field(g_tempgauge_field, g_gauge_field);
  
  return(0);
  
}




/*

// MPI implementation									// was merged into init_temporalgauge_without_mpi()

#ifdef MPI

int init_temporalgauge_trafo_mpi (const int V, su3 ** gfield) {				// will initialize  g_trafo[]  as the transformation matrices
											//	and  g_tempgauge_field  as a copy of  g_gauge_field
  int it, iz, iy, ix;
  
  int pos;
  
  MPI_Status status;
  
  
  
  if ((void *)(left = (su3 *) calloc(LX*LY*LZ, sizeof(su3))) == NULL ) {		// allocates memory for a time-slice of su3-matrices
    printf("malloc error in 'init_temporalgauge_trafo_mpi'\n"); 
    return(-1);
  }
  
  if ((void *)(right = (su3 *) calloc(LX*LY*LZ, sizeof(su3))) == NULL ) {		// allocates memory for a time-slice of su3-matrices
    printf("malloc error in 'init_temporalgauge_trafo_mpi'\n"); 
    return(-1);
  }
  
  
  
  
  if ((void *)(g_trafo = (su3 *) calloc(V, sizeof(su3))) == NULL ) {			// allocates memory for V su3-matrices
    printf("malloc error in 'init_temporalgauge_trafo'\n"); 
    return(2);
  } 
  
  
  
  
  //////////////////////////////////////////////
  // initializing the transformation matrices //
  //////////////////////////////////////////////
  
  
  // first process in t-direction
  
  if (g_cart_id == 0) {
  	
  	// initialize first timeslice (t=0) with unit matrices
  	for (ix = 0; ix < LX; ix++) {
  	  for (iy = 0; iy < LY; iy++) {
  	    for (iz = 0; iz < LZ; iz++) {
  	      g_trafo[g_ipt[0][ix][iy][iz]] = unit_su3();					// g_trafo[0-th time slice]  =  ID
  	    }
  	  }
  	}
  	
  	// U^{'}_0(x)  =  g(x) U_0(x) g^{+}(x+0)  !=  ID   =>   g_(x+0)  =  g(x) U_0(x)
  	for (it = 1; it < T; it++) {
  	  for (ix = 0; ix < LX; ix++) {
  	    for (iy = 0; iy < LY; iy++) {
  	      for (iz = 0; iz < LZ; iz++) {
  	        _su3_times_su3( g_trafo[ g_ipt[it  ][ix][iy][iz] ] , 				// g_trafo[next t-slice]  =  g_trafo[old t-slice]  *  gfield[old t-slice][t-dir.]
  	                        g_trafo[ g_ipt[it-1][ix][iy][iz] ] ,
  	                        gfield [ g_ipt[it-1][ix][iy][iz] ] [0] );
  	        
  	      }
  	    }
  	  } 
  	}
  	
  	
  	// sending
  	MPI_Send((void *)(g_trafo+(T-1)*LX*LY*LZ), LX*LY*LZ, mpi_su3, g_nb_t_up, 0, g_cart_grid);
  	//MPI_Send((void *)(g_trafo+(T-1)*LX*LY*LZ), LX*LY*LZ, mpi_su3, g_cart_id+1, 0, g_cart_grid);
  	
  	printf("g_cart_id = %i has send a message to %i\n", g_cart_id, g_nb_t_up);
  	
  	
  } // first process
  
  
  
  
  // following processes
  
  else {
  	
  	// receiving
  	MPI_Recv((void *)left, LX*LY*LZ, mpi_su3, g_nb_t_dn, 0, g_cart_grid, &status);
  	//MPI_Recv((void *)left, LX*LY*LZ, mpi_su3, g_cart_id-1, 0, g_cart_grid, &status);
  	
  	
  	printf("g_cart_id = %i has received a message from %i\n", g_cart_id, g_nb_t_dn);
  	
  	it = 0;
  	for (ix = 0; ix < LX; ix++) {
  	  for (iy = 0; iy < LY; iy++) {
  	    for (iz = 0; iz < LZ; iz++) {
  	      pos = g_ipt[it][ix][iy][iz];
  	      _su3_times_su3( g_trafo[ g_ipt[it  ][ix][iy][iz] ] ,				// g_trafo[0-th time slice]  =  left[xchanged t-slice]  * gfield[
  	                      left   [ g_ipt[it  ][ix][iy][iz] ] ,
  	                      gfield [ g_idn[pos ][0]          ] [0] );				// notice: have to access the RAND region of the gauge field
  	    }
  	  }
  	}
  	
  	
  	for (it = 1; it < T; it++) {
  	  for (ix = 0; ix < LX; ix++) {
  	    for (iy = 0; iy < LY; iy++) {
  	      for (iz = 0; iz < LZ; iz++) {
  	        _su3_times_su3( g_trafo[ g_ipt[it  ][ix][iy][iz] ] ,
  	                        g_trafo[ g_ipt[it-1][ix][iy][iz] ] ,
  	                        gfield [ g_ipt[it-1][ix][iy][iz] ] [0] );
  	        
  	      }
  	    }
  	  } 
  	}
  	
  	
  	// sending
  	if (g_cart_id != g_nproc-1) {
  	  MPI_Send((void *)(g_trafo+(T-1)*LX*LY*LZ), LX*LY*LZ, mpi_su3, g_nb_t_up, 0, g_cart_grid);
  	  //MPI_Send((void *)(g_trafo+(T-1)*LX*LY*LZ), LX*LY*LZ, mpi_su3, g_cart_id+1, 0, g_cart_grid);
  	  
  	  printf("g_cart_id = %i has send a message to %i\n", g_cart_id, g_nb_t_up);
  	  
  	}
  	
  
  } // following processes
  
  
  
  
  ////////////////////////////////////////////
  // exchanging the transformation matrices //
  ////////////////////////////////////////////
  
  
  MPI_Sendrecv((void *)(g_trafo), LX*LY*LZ, mpi_su3, g_nb_t_dn, 1,
               (void *)(right  ), LX*LY*LZ, mpi_su3, g_nb_t_up, 1,
               g_cart_grid, &status);
  
  printf("g_cart_id = %i has send to %i and received from %i\n", g_cart_id, g_nb_t_dn, g_nb_t_up);
  
  
  
  
  
  // all processes
  
  // copying the gaugefield (for later undoing the transformation)
  
  if ((void *)(g_tempgauge_field = (su3 **) calloc(V, sizeof(su3*))) == NULL ) {	// allocates  V  su3 *
    printf ("malloc error in 'init_temporalgauge_trafo'\n"); 
  return(1);
  }
  if ((void *)(tempgauge_field = (su3 *) calloc(4*V+1, sizeof(su3))) == NULL ) {	// allocates  4*V+1  su3-matrices
    printf ("malloc error in 'init_temporalgauge_trafo'\n"); 
    return(2);
  }
  
  #if (defined SSE || defined SSE2 || defined SSE3)
    g_tempgauge_field[0] = (su3*)(((unsigned long int)(tempgauge_field)+ALIGN_BASE)&~ALIGN_BASE);
  #else
    g_tempgauge_field[0] = tempgauge_field;
  #endif
  
  int i = 0;
  
  for (i = 1; i < V; i++) {
    g_tempgauge_field[i] = g_tempgauge_field[i-1]+4;
  }
  
  // copy the original field
  copy_gauge_field(g_tempgauge_field, g_gauge_field);
  
  
  return(0);
  
  
}//init_temporalgauge_trafo_mpi()

#endif //MPI

*/





void finalize_temporalgauge() {

  free(g_trafo);
  free(tempgauge_field);
  free(g_tempgauge_field);
  
  #ifdef MPI
    free(left);
    free(right);
  #endif

}





/*

//  apply gauge transform to gfield with the trafo stored in trafofield

void apply_gtrafo2 (su3 ** gfield, su3 * trafofield) {

  int it, iz, iy, ix, xpos, mu;
  
  su3 temp1;
  
  if (g_proc_id == 0) {
    printf("Applying gauge transformation...");
  }
  
  for (it = 0; it < T; it++) {
    for (ix = 0; ix < LX; ix++) {
      for (iy = 0; iy < LY; iy++) {
        for (iz = 0; iz < LZ; iz++) {
        
          xpos = g_ipt[it][ix][iy][iz];
          
          for (mu = 0; mu < 4; mu++) {
            // help = g(x) U_mu(x)
            _su3_times_su3( temp1, trafofield[xpos],  gfield[xpos][mu]  );
            // U_mu(x) <- U_mu^{'}(x) = help g^{+}(x+mu)
            _su3_times_su3d( gfield[xpos][mu],temp1, trafofield[ g_iup[xpos][mu]  ]);
          }
          
        }
      }
    } 
  }
  
  if (g_proc_id == 0) {
    printf("done\n");
  }
  
  // update gauge copy fields in the next call to HoppingMatrix
  g_update_gauge_copy = 1;
}

*/





//  apply gauge transform to gfield with the trafo stored in trafofield

void apply_gtrafo (su3 ** gfield, su3 * trafofield) {

  int it, iz, iy, ix;
  int pos;
  int mu;
  
  su3 temp1;
  
  if (g_proc_id == 0) {
    printf("Applying gauge transformation...");
  }
  
  for (it = 0; it < T; it++) {
    for (ix = 0; ix < LX; ix++) {
      for (iy = 0; iy < LY; iy++) {
        for (iz = 0; iz < LZ; iz++) {
        
          #ifdef MPI				// this is the MPI implementation of the GLOBAL TEMPORALGAUGE
          
            pos = g_ipt[it][ix][iy][iz];
            
            for (mu = 0; mu < 4; mu++) {
              if ((it != T-1) || (mu != 0)) {
                /* help = g(x) U_mu(x) */
                _su3_times_su3( temp1, trafofield[pos],  gfield[pos][mu]  );			// temp1  =  trafofield[pos]  *  gfield[pos][mu]
                /* U_mu(x) <- U_mu^{'}(x) = help g^{+}(x+mu)*/
                _su3_times_su3d( gfield[pos][mu],temp1, trafofield[ g_iup[pos][mu]  ]);		// gfield[pos][mu]  =  temp1  *  trafofield[ g_iup[pos][mu] ] _ {dagger}
              }											//                  =  trafofield[pos] * gfield[pos][mu] * trafofield[ g_iup[pos][mu] ]_{dagger}
              else {
                _su3_times_su3( temp1, trafofield[pos],  gfield[pos][mu]  );
                _su3_times_su3d( gfield[pos][mu],temp1, right[ g_ipt[0][ix][iy][iz]  ]);	// "rightest" transf. matrices are stored in right[]
              }
            }
            
          #else					// in case of using this version with MPI this is
          					// a LOCAL version of TEMPORALGAUGE
            pos = g_ipt[it][ix][iy][iz];
            
            for (mu = 0; mu < 4; mu++) {
              if ((it != T-1) || (mu != 0)) {
                /* help = g(x) U_mu(x) */
                _su3_times_su3( temp1, trafofield[pos],  gfield[pos][mu]  );
                /* U_mu(x) <- U_mu^{'}(x) = help g^{+}(x+mu)*/
                _su3_times_su3d( gfield[pos][mu],temp1, trafofield[ g_iup[pos][mu]  ]);
              }
              else { // (it = T-1) && (mu = 0)
                _su3_times_su3( temp1, trafofield[pos],  gfield[pos][mu]  );
                _su3_times_su3d( gfield[pos][mu],temp1, trafofield[ g_ipt[0][ix][iy][iz]  ]);	// "rightest" transf. matrices are the first (periodic) and are initialized to ID
              }
            }
            
          #endif
          
        }
      }
    } 
  }
  
  if (g_proc_id == 0) {
    printf("done\n");
  }
  
  /* update gauge copy fields in the next call to HoppingMatrix */
  g_update_gauge_copy = 1;
 
}//apply_gtrafo()




/*
  apply the inverse gauge transform to gfield with the trafo stored in trafofield
*/

// this is not really needed, instead we are copying the original gauge field

void apply_inv_gtrafo (su3 ** gfield, su3 * trafofield) {

 int it, iz, iy, ix;
 int xpos;
 int mu;
  
 su3 temp1, temp2;
 
 if(g_proc_id == 0) {
   printf("Applying INVERSE gauge transformation...");
 }
 
 for (it = 0; it < T; it++) {
   for (ix = 0; ix < LX; ix++) {
     for (iy = 0; iy < LY; iy++) {
       for (iz = 0; iz < LZ; iz++) {
       
         xpos = g_ipt[it][ix][iy][iz];
         
         for (mu = 0; mu < 4; mu++) {
            /*
            _su3d_times_su3( temp1, trafofield[xpos],  gfield[xpos][mu]  );

            _su3_times_su3( gfield[xpos][mu],temp1, trafofield[ g_iup[xpos][mu]  ]);            
            */
           
           /* help = U^{'}_mu(x) g(x+mu)*/
            _su3_times_su3( temp1,  gfield[xpos][mu], trafofield[ g_iup[xpos][mu]]  );	// temp1  =  gfield[xpos][mu]  *  trafofield[ g_iup[xpos][mu] ]

            /* U_mu(x) <- g^{+}(x) help */
            _su3_dagger(temp2, trafofield[xpos]  )					// temp2  =  trafofield[xpos]_{dagger}
            _su3_times_su3( gfield[xpos][mu], temp2, temp1);				// gfield[xpos][mu]  =  temp2  * temp1
											//                   =  trafofield[xpos]_{dagger}  *  gfield[xpos][mu]  *  trafofield[ g_iup[xpos][mu] ]
          }  
  }}}}
  
  if(g_proc_id == 0) {
    printf("done\n");
  }
  
  /* update gauge copy fields in the next call to HoppingMatrix */
  g_update_gauge_copy = 1;
  
}



/* 
  apply inverse gauge transform to spinor 
  U_0(x) = g^{+}(x) U^{'}_0(x) g(x+0)
  => psi(x) = g^{+}(x) psi^{'}(x)
  (the primed (^{'}) quantities are the gauge transformed fields)
*/

void apply_inv_gtrafo_spinor (spinor * spin, su3 * trafofield) {

 int it, iz, iy, ix;
 int pos;
  
 spinor temp;
 
 if(g_proc_id == 0) {
   printf("Applying INVERSE gauge transformation to spinor...");
 }
 
 for (it = 0; it < T; it++) {
  for (ix = 0; ix < LX; ix++) {
    for (iy = 0; iy < LY; iy++) {
      for (iz = 0; iz < LZ; iz++) {
      
        pos = g_ipt[it][ix][iy][iz];
        
        _su3_inverse_multiply(temp.s0, trafofield[pos], spin[pos].s0);
        _su3_inverse_multiply(temp.s1, trafofield[pos], spin[pos].s1);
        _su3_inverse_multiply(temp.s2, trafofield[pos], spin[pos].s2);
        _su3_inverse_multiply(temp.s3, trafofield[pos], spin[pos].s3);
        
        _vector_assign(spin[pos].s0,temp.s0);
        _vector_assign(spin[pos].s1,temp.s1);
        _vector_assign(spin[pos].s2,temp.s2);
        _vector_assign(spin[pos].s3,temp.s3);
      
        }
      }
    } 
  }
  
  if (g_proc_id == 0) {
    printf("done\n");
  }
  
}




/* 
  apply gauge transform to spinor 
  U^{'}_0(x) = g(x) U_0(x) g^{+}(x+0)
  => psi^{'}(x) = g(x) psi(x)
  (the primed (^{'}) quantities are the gauge transformed fields)
*/

void apply_gtrafo_spinor (spinor * spin, su3 * trafofield) {

 int it, iz, iy, ix;
 int pos; 
 spinor temp;
 
 if(g_proc_id == 0) {
   printf("Applying gauge transformation to spinor...");
 }
 
  for (it = 0; it < T; it++) {
    for (ix = 0; ix < LX; ix++) {
      for (iy = 0; iy < LY; iy++) {
        for (iz = 0; iz < LZ; iz++) {
        
          pos = g_ipt[it][ix][iy][iz];
          
          _su3_multiply(temp.s0, trafofield[pos], spin[pos].s0);
          _su3_multiply(temp.s1, trafofield[pos], spin[pos].s1);
          _su3_multiply(temp.s2, trafofield[pos], spin[pos].s2);
          _su3_multiply(temp.s3, trafofield[pos], spin[pos].s3);
          
          _vector_assign(spin[pos].s0,temp.s0);
          _vector_assign(spin[pos].s1,temp.s1);
          _vector_assign(spin[pos].s2,temp.s2);
          _vector_assign(spin[pos].s3,temp.s3);
        }
      }
    } 
  }
  
  if(g_proc_id == 0) {
    printf("done\n");
  }
  
}




/* 
  apply gauge transform to ODD spinor 
  U^{'}_0(x) = g(x) U_0(x) g^{+}(x+0)
  => psi^{'}(x) = g(x) psi(x)
  (the primed (^{'}) quantities are the gauge transformed fields)
*/

void apply_gtrafo_spinor_odd (spinor * spin, su3 * trafofield) {

  int it, iz, iy, ix;
  int pos;
  int oddpos; 
  spinor temp;

  if (g_proc_id == 0) {
    printf("Applying  gauge transformation to odd spinor...");
  }
  
  for (it = 0; it < T; it++) {
    for (ix = 0; ix < LX; ix++) { 
      for (iy = 0; iy < LY; iy++) {
        for (iz = 0; iz < LZ; iz++) {
        
          if ((it + ix + iy + iz) % 2 != 0) {
            /* odd positions */
            pos = g_ipt[it][ix][iy][iz];
            oddpos = g_lexic2eosub[ pos  ];

            _su3_multiply(temp.s0, trafofield[pos], spin[oddpos].s0);
            _su3_multiply(temp.s1, trafofield[pos], spin[oddpos].s1);
            _su3_multiply(temp.s2, trafofield[pos], spin[oddpos].s2);
            _su3_multiply(temp.s3, trafofield[pos], spin[oddpos].s3);
            
            _vector_assign(spin[oddpos].s0, temp.s0);
            _vector_assign(spin[oddpos].s1, temp.s1);
            _vector_assign(spin[oddpos].s2, temp.s2);
            _vector_assign(spin[oddpos].s3, temp.s3);
            
          }
        }
      }
    } 
  }
  
  if (g_proc_id == 0) {
    printf("done\n");
  }
  
}




/* 
  apply inverse gauge transform to ODD spinor 
  U_0(x) = g^{+}(x) U^{'}_0(x) g(x+0)
  => psi(x) = g^{+}(x) psi^{'}(x)
  (the primed (^{'}) quantities are the gauge ttemp.s0ransformed fields)
*/

void apply_inv_gtrafo_spinor_odd (spinor * spin, su3 * trafofield) {

  int it, iz, iy, ix;
  int pos;
  int oddpos;
  
  spinor temp;
  
  if (g_proc_id == 0) {
    printf("Applying INVERSE gauge transformation to odd spinor...");
  }
  for (it = 0; it < T; it++) {
    for (ix = 0; ix < LX; ix++) {
      for (iy = 0; iy < LY; iy++) {
        for (iz = 0; iz < LZ; iz++) {
        
          if ((it + ix + iy + iz) % 2 != 0) {
          
            /* odd positions */
            pos = g_ipt[it][ix][iy][iz];
            oddpos = g_lexic2eosub[ pos ];
            
            _su3_inverse_multiply(temp.s0, trafofield[pos], spin[oddpos].s0);
            _su3_inverse_multiply(temp.s1, trafofield[pos], spin[oddpos].s1);
            _su3_inverse_multiply(temp.s2, trafofield[pos], spin[oddpos].s2);
            _su3_inverse_multiply(temp.s3, trafofield[pos], spin[oddpos].s3);
            
            _vector_assign(spin[oddpos].s0, temp.s0);
            _vector_assign(spin[oddpos].s1, temp.s1);
            _vector_assign(spin[oddpos].s2, temp.s2);
            _vector_assign(spin[oddpos].s3, temp.s3);
            
          }
        
        }
      }
    } 
  }
  
  if (g_proc_id == 0) {
    printf("done\n");
  }
  
}




/* 
  apply gauge transform to EVENspinor 
  U^{'}_0(x) = g(x) U_0(x) g^{+}(x+0)
  => psi^{'}(x) = g(x) psi(x)
  (the primed (^{'}) quantities are the gauge transformed fields)
*/

void apply_gtrafo_spinor_even (spinor * spin, su3 * trafofield) {

  int it, iz, iy, ix;
  int pos;
  int evenpos;
  
  spinor temp;
  
  if (g_proc_id == 0) {
    printf("Applying  gauge transformation to even spinor...");
  }
  
  for (it = 0; it < T; it++) {
    for (ix = 0; ix < LX; ix++) { 
      for (iy = 0; iy < LY; iy++) {
        for (iz = 0; iz < LZ; iz++) {
         
          if ((it + ix + iy + iz) % 2 == 0) {
          
            /* even positions */
            pos = g_ipt[it][ix][iy][iz];
            evenpos = g_lexic2eosub[ pos  ];

            _su3_multiply(temp.s0, trafofield[pos], spin[evenpos].s0);
            _su3_multiply(temp.s1, trafofield[pos], spin[evenpos].s1);
            _su3_multiply(temp.s2, trafofield[pos], spin[evenpos].s2);
            _su3_multiply(temp.s3, trafofield[pos], spin[evenpos].s3);
            
            _vector_assign(spin[evenpos].s0, temp.s0);
            _vector_assign(spin[evenpos].s1, temp.s1);
            _vector_assign(spin[evenpos].s2, temp.s2);
            _vector_assign(spin[evenpos].s3, temp.s3);
            
          }
        }
      }
    } 
  }
  
  if (g_proc_id == 0) {
    printf("done\n");
  }
  
}



/* 
  apply inverse gauge transform to EVEN spinor 
  U_0(x) = g^{+}(x) U^{'}_0(x) g(x+0)
  => psi(x) = g^{+}(x) psi^{'}(x)
  (the primed (^{'}) quantities are the gauge transformed fields)
*/
void apply_inv_gtrafo_spinor_even (spinor * spin, su3 * trafofield) {

  int it, iz, iy, ix;
  int xpos;
  int evenpos;
  
  spinor temp;
  
  if (g_proc_id == 0) {
    printf("Applying INVERSE gauge transformation to even spinor...");
  }
  for (it = 0; it < T; it++) {
    for (ix = 0; ix < LX; ix++) {
      for (iy = 0; iy < LY; iy++) {
        for (iz = 0; iz < LZ; iz++) {
        
          if ((it+ix+iy+iz)%2 == 0) {
            /* even positions */
            xpos = g_ipt[it][ix][iy][iz];
            evenpos = g_lexic2eosub[ xpos ];
            
            _su3_inverse_multiply(temp.s0, trafofield[xpos], spin[evenpos].s0);
            _su3_inverse_multiply(temp.s1, trafofield[xpos], spin[evenpos].s1);
            _su3_inverse_multiply(temp.s2, trafofield[xpos], spin[evenpos].s2);
            _su3_inverse_multiply(temp.s3, trafofield[xpos], spin[evenpos].s3);
            
            _vector_assign(spin[evenpos].s0,temp.s0);
            _vector_assign(spin[evenpos].s1,temp.s1);
            _vector_assign(spin[evenpos].s2,temp.s2);
            _vector_assign(spin[evenpos].s3,temp.s3);
            
          }
        
        }
      }
    } 
  }
  
  if (g_proc_id == 0) {
    printf("done\n");
  }
  
}








