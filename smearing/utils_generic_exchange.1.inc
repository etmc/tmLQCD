#    if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  
  MPI_Isend(buffer[gI_0_0_0_0], 1, slice_T_cont_type, g_nb_t_dn, 83,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_L_0_0_0], 1, slice_T_cont_type, g_nb_t_up, 83, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Isend(buffer[gI_Lm1_0_0_0], 1, slice_T_cont_type, g_nb_t_up, 84,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_m1_0_0_0], 1, slice_T_cont_type, g_nb_t_dn, 84, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    /* t2-Rand */
    MPI_Isend(buffer[gI_p1_0_0_0], 1, slice_T_cont_type, g_nb_t_dn, 85, 
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_Lp1_0_0_0], 1, slice_T_cont_type, g_nb_t_up, 85, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    /* t2-Rand */
    MPI_Isend(buffer[gI_Lm2_0_0_0], 1, slice_T_cont_type, g_nb_t_up, 86, 
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_m2_0_0_0], 1, slice_T_cont_type, g_nb_t_dn, 86, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
#    endif
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Isend(buffer[gI_0_0_0_0], 1, slice_X_gath_type, g_nb_x_dn, 87,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_L_0_0], 1, slice_X_cont_type, g_nb_x_up, 87, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* x-Rand */
  MPI_Isend(buffer[gI_0_Lm1_0_0], 1, slice_X_gath_type, g_nb_x_up, 88,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_m1_0_0], 1, slice_X_cont_type, g_nb_x_dn, 88,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2-Rand */
    MPI_Isend(buffer[gI_0_p1_0_0], 1, slice_X_gath_type, g_nb_x_dn, 89,
	     g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_Lp1_0_0], 1, slice_X_cont_type, g_nb_x_up, 89, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2-Rand */
    MPI_Isend(buffer[gI_0_Lm2_0_0], 1, slice_X_gath_type, g_nb_x_up, 90,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_m2_0_0], 1, slice_X_cont_type, g_nb_x_dn, 90,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
#    endif

  MPI_Waitall(cntr, request, status);
  cntr=0;

  /* Communications of the xt (x2t and t2x) edges are done by using the previously 
     communicated x-borders whose t-borders are now exchanged in t directions [ORD!] */
  /* In this case the code cannot be completely independent of the definition in Index, 
     since edge_XT_gath_type are defined by joining together the x=L and the x=-1 parts.
     For this reason we need to know that x=L comes before x=-1 in the definition of 
     Index() and hence we need to refer to the starting point gI_0_L_0_0 . [DEP!] */

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  /* is on the x-Rand: xt-edge */
  MPI_Isend(buffer[gI_0_L_0_0], 1, edge_XT_gath_type, g_nb_t_dn, 100,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_L_L_0_0], 1, edge_XT_cont_type, g_nb_t_up, 100, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  /* xt-edge */
  MPI_Isend(buffer[gI_Lm1_L_0_0], 1, edge_XT_gath_type, g_nb_t_up, 101,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_m1_L_0_0], 1, edge_XT_cont_type, g_nb_t_dn, 101,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* t2x-edge */
    MPI_Isend(buffer[gI_p1_L_0_0], 1, edge_XT_gath_type, g_nb_t_dn, 102,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_Lp1_L_0_0], 1, edge_XT_cont_type, g_nb_t_up, 102, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* t2x-edge */
    MPI_Isend(buffer[gI_Lm2_L_0_0], 1, edge_XT_gath_type, g_nb_t_up, 103,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_m2_L_0_0], 1, edge_XT_cont_type, g_nb_t_dn, 103,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* x2t-edge */   /* x=L+1 comes before x=-2. see [DEP!]  */
    MPI_Isend(buffer[gI_0_Lp1_0_0], 1, edge_XT_gath_type, g_nb_t_dn, 104,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_L_Lp1_0_0], 1, edge_XT_cont_type, g_nb_t_up, 104, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* x2t-edge */
    MPI_Isend(buffer[gI_Lm1_Lp1_0_0], 1, edge_XT_gath_type, g_nb_t_up, 105,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_m1_Lp1_0_0], 1, edge_XT_cont_type, g_nb_t_dn, 105,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
#    endif

#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend(buffer[gI_0_0_0_0], 1, slice_Y_gath_type, g_nb_y_dn, 106,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_0_L_0], 1, slice_Y_cont_type, g_nb_y_up, 106, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  MPI_Isend(buffer[gI_0_0_Lm1_0], 1, slice_Y_gath_type, g_nb_y_up, 107,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_0_m1_0], 1, slice_Y_cont_type, g_nb_y_dn, 107,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* y2-Rand */
    MPI_Isend(buffer[gI_0_0_p1_0], 1, slice_Y_gath_type, g_nb_y_dn, 108,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_0_Lp1_0], 1, slice_Y_cont_type, g_nb_y_up, 108, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* y2-Rand */
    MPI_Isend(buffer[gI_0_0_Lm2_0], 1, slice_Y_gath_type, g_nb_y_up, 109,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_0_m2_0], 1, slice_Y_cont_type, g_nb_y_dn, 109,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
#    endif

  MPI_Waitall(cntr, request, status);
  cntr=0;

  /* see [ORD!] above, where now x plays the role of t and y the role of x */
  /* see [DEP!] above, where now y=L comes before y=-1 */ 

#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the y-Rand -> yx-edge*/
  MPI_Isend(buffer[gI_0_0_L_0], 1, edge_XY_gath_type, g_nb_x_dn, 110,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_L_L_0], 1, edge_XY_cont_type, g_nb_x_up, 110, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* yx-edge */
  MPI_Isend(buffer[gI_0_Lm1_L_0], 1, edge_XY_gath_type, g_nb_x_up, 111,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_m1_L_0],  1, edge_XY_cont_type, g_nb_x_dn, 111,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

#    endif

  /* see [ORD!] above, where now y plays the role of t and t the role of x */
  /* see [DEP!] above, where now t=L comes before t=-1 */ 

#    if (defined PARALLELXYT || defined PARALLELXYZT )

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the t-Rand -> ty-edge*/
  MPI_Isend(buffer[gI_L_0_0_0], 1, edge_YT_gath_type, g_nb_y_dn, 112,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_L_0_L_0], 1, edge_YT_cont_type, g_nb_y_up, 112, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* ty-edge */
  MPI_Isend(buffer[gI_L_0_Lm1_0], 1, edge_YT_gath_type, g_nb_y_up, 113,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_L_0_m1_0], 1, edge_YT_cont_type, g_nb_y_dn, 113,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

#    endif

  if(g_dbw2rand > 0) {

#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* x2y edge */  /* y=L comes before y=-1 */
    MPI_Isend(buffer[gI_0_p1_L_0], 1, edge_XY_gath_type, g_nb_x_dn, 114,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_Lp1_L_0], 1, edge_XY_cont_type, g_nb_x_up, 114, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* x2y-edge */
    MPI_Isend(buffer[gI_0_Lm2_L_0], 1, edge_XY_gath_type, g_nb_x_up, 115,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_m2_L_0], 1, edge_XY_cont_type, g_nb_x_dn, 115,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;


    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2x -edge */
    MPI_Isend(buffer[gI_0_0_Lp1_0], 1, edge_XY_gath_type, g_nb_x_dn, 116,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_L_Lp1_0], 1, edge_XY_cont_type, g_nb_x_up, 116, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2x edge */
    MPI_Isend(buffer[gI_0_Lm1_Lp1_0], 1, edge_XY_gath_type, g_nb_x_up, 117,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_m1_Lp1_0], 1, edge_XY_cont_type, g_nb_x_dn, 117,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

#    endif
#    if (defined PARALLELXYT || defined PARALLELXYZT )

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* t2y-edge */
    MPI_Isend(buffer[gI_Lp1_0_0_0], 1, edge_YT_gath_type, g_nb_y_dn, 118,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_Lp1_0_L_0], 1, edge_YT_cont_type, g_nb_y_up, 118, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* t2y edge */
    MPI_Isend(buffer[gI_Lp1_0_Lm1_0], 1, edge_YT_gath_type, g_nb_y_up, 119,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_Lp1_0_m1_0], 1, edge_YT_cont_type, g_nb_y_dn, 119,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2t edge */
    MPI_Isend(buffer[gI_L_0_p1_0], 1, edge_YT_gath_type, g_nb_y_dn, 120,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_L_0_Lp1_0], 1, edge_YT_cont_type, g_nb_y_up, 120, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2t-edge */
    MPI_Isend(buffer[gI_L_0_Lm2_0], 1, edge_YT_gath_type, g_nb_y_up, 121,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_L_0_m2_0], 1, edge_YT_cont_type, g_nb_y_dn, 121,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
#    endif
  }
#    if (defined PARALLELXYZT || defined PARALLELXYZ )
  /* z-Rand */
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend(buffer[gI_0_0_0_0], 1, slice_Z_gath_type, g_nb_z_dn, 122,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_0_0_L], 1, slice_Z_cont_type, g_nb_z_up, 122, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Isend(buffer[gI_0_0_0_Lm1], 1, slice_Z_gath_type, g_nb_z_up, 123,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_0_0_m1], 1, slice_Z_cont_type, g_nb_z_dn, 123,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2-Rand */
    MPI_Isend(buffer[gI_0_0_0_p1], 1, slice_Z_gath_type, g_nb_z_dn, 124,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_0_0_Lp1], 1, slice_Z_cont_type, g_nb_z_up, 124, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2-Rand */
    MPI_Isend(buffer[gI_0_0_0_Lm2], 1, slice_Z_gath_type, g_nb_z_up, 125,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_0_0_m2], 1, slice_Z_cont_type, g_nb_z_dn, 125,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
#    endif
  MPI_Waitall(cntr, request, status);
  cntr=0;

  /* see [ORD!] above, where now x plays the role of t and z the role of x */
  /* see [DEP!] above, where now z=L comes before z=-1 */ 

#    if (defined PARALLELXYZT || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the z-Rand -> zx-edge*/
  MPI_Isend(buffer[gI_0_0_0_L], 1, edge_XZ_gath_type, g_nb_x_dn, 126,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_L_0_L], 1, edge_XZ_cont_type, g_nb_x_up, 126, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* zx-edge */
  MPI_Isend(buffer[gI_0_Lm1_0_L], 1, edge_XZ_gath_type, g_nb_x_up, 127,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_m1_0_L],
	    1, edge_XZ_cont_type, g_nb_x_dn, 127,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

#    endif

  /* see [ORD!] above, where now z plays the role of t and t the role of x */
  /* see [DEP!] above, where now t=L comes before t=-1 */ 

#    if (defined PARALLELXYZT)
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  /* is on the t-Rand -> tz-edge*/
  MPI_Isend(buffer[gI_L_0_0_0], 1, edge_ZT_gath_type, g_nb_z_dn, 128,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_L_0_0_L], 1, edge_ZT_cont_type, g_nb_z_up, 128, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;
  
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  /* tz-edge */
  MPI_Isend(buffer[gI_L_0_0_Lm1], 1, edge_ZT_gath_type, g_nb_z_up, 129,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_L_0_0_m1], 1, edge_ZT_cont_type, g_nb_z_dn, 129,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

#    endif

  /* see [ORD!] above, where now y plays the role of t and z the role of x */
  /* see [DEP!] above, where now z=L comes before z=-1 */ 

#    if (defined PARALLELXYZT || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the z-Rand -> zy-edge*/
  MPI_Isend(buffer[gI_0_0_0_L], 1, edge_YZ_gath_type, g_nb_y_dn, 130,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_0_L_L], 1, edge_YZ_cont_type, g_nb_y_up, 130, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* zy-edge */
  MPI_Isend(buffer[gI_0_0_Lm1_L], 1, edge_YZ_gath_type, g_nb_y_up, 131,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[gI_0_0_m1_L], 1, edge_YZ_cont_type, g_nb_y_dn, 131,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

#    endif

  if(g_dbw2rand > 0) {

#    if (defined PARALLELXYZT)
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* t2z edge */ /* t=L+1 comes before t=-2*/
    MPI_Isend(buffer[gI_Lp1_0_0_0], 1, edge_ZT_gath_type, g_nb_z_dn, 132,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_Lp1_0_0_L], 1, edge_ZT_cont_type, g_nb_z_up, 132, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* t2z-edge */
    MPI_Isend(buffer[gI_Lp1_0_0_Lm1], 1, edge_ZT_gath_type, g_nb_z_up, 133,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_Lp1_0_0_m1], 1, edge_ZT_cont_type, g_nb_z_dn, 133,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;


    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2t -edge */
    MPI_Isend(buffer[gI_L_0_0_p1], 1, edge_ZT_gath_type, g_nb_z_dn, 134,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_L_0_0_Lp1], 1, edge_ZT_cont_type, g_nb_z_up, 134, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2t edge */
    MPI_Isend(buffer[gI_L_0_0_Lm2], 1, edge_ZT_gath_type, g_nb_z_up, 135,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_L_0_0_m2], 1, edge_ZT_cont_type, g_nb_z_dn, 135,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

#    endif
#    if (defined PARALLELXYZT || defined PARALLELXYZ )
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* z2x-edge */
    MPI_Isend(buffer[gI_0_0_0_Lp1], 1, edge_XZ_gath_type, g_nb_x_dn, 136,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_L_0_Lp1], 1, edge_XZ_cont_type, g_nb_x_up, 136, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* z2x edge */
    MPI_Isend(buffer[gI_0_Lm1_0_Lp1], 1, edge_XZ_gath_type, g_nb_x_up, 137,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_m1_0_Lp1], 1, edge_XZ_cont_type, g_nb_x_dn, 137,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2z edge */
    MPI_Isend(buffer[gI_0_p1_0_L], 1, edge_XZ_gath_type, g_nb_x_dn, 138,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_Lp1_0_L], 1, edge_XZ_cont_type, g_nb_x_up, 138, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2z-edge */
    MPI_Isend(buffer[gI_0_Lm2_0_L], 1, edge_XZ_gath_type, g_nb_x_up, 139,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_m2_0_L], 1, edge_XZ_cont_type, g_nb_x_dn, 139,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

#    endif
#    if (defined PARALLELXYZT || defined PARALLELXYZ )
    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* z2y-edge */ /* z=L+1 comes before z=-2 */
    MPI_Isend(buffer[gI_0_0_0_Lp1], 1, edge_YZ_gath_type, g_nb_y_dn, 140,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_0_L_Lp1], 1, edge_YZ_cont_type, g_nb_y_up, 140, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* z2y edge */
    MPI_Isend(buffer[gI_0_0_Lm1_Lp1], 1, edge_YZ_gath_type, g_nb_y_up, 141,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_0_m1_Lp1], 1, edge_YZ_cont_type, g_nb_y_dn, 141,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2z edge */ /* z=L comes before z=-1 */
    MPI_Isend(buffer[gI_0_0_p1_L], 1, edge_YZ_gath_type, g_nb_y_dn, 142,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_0_Lp1_L], 1, edge_YZ_cont_type, g_nb_y_up, 142, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2z-edge */
    MPI_Isend(buffer[gI_0_0_Lm2_L], 1, edge_YZ_gath_type, g_nb_y_up, 143,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[gI_0_0_m2_L], 1, edge_YZ_cont_type, g_nb_y_dn, 143,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
#    endif
  }

#    if (defined PARALLELXYZT || defined PARALLELXYZ )
  MPI_Waitall(cntr, request, status);
#    endif
