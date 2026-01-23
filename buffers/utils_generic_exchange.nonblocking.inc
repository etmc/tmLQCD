  MPI_Isend(buffer[0],          1, slice_T_cont_type, g_nb_t_dn, 83,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME], 1, slice_T_cont_type, g_nb_t_up, 83, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Isend(buffer[(T-1)*LX*LY*LZ], 1, slice_T_cont_type, g_nb_t_up, 84,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[(T+1)*LX*LY*LZ], 1, slice_T_cont_type, g_nb_t_dn, 84, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    /* t2-Rand */
    MPI_Isend(buffer[1*LX*LY*LZ],     1, slice_T_cont_type, g_nb_t_dn, 85, 
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND], 1, slice_T_cont_type, g_nb_t_up, 85, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    /* t2-Rand */
    MPI_Isend(buffer[(T-2)*LX*LY*LZ],          1, slice_T_cont_type, g_nb_t_up, 86, 
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND+LX*LY*LZ], 1, slice_T_cont_type, g_nb_t_dn, 86, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
  
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Isend(buffer[0],              1, slice_X_gath_type, g_nb_x_dn, 87,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[(T+2)*LX*LY*LZ], 1, slice_X_cont_type, g_nb_x_up, 87, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* x-Rand */
  MPI_Isend(buffer[(LX-1)*LY*LZ],             1, slice_X_gath_type, g_nb_x_up, 88,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[(T+2)*LX*LY*LZ + T*LY*LZ], 1, slice_X_cont_type, g_nb_x_dn, 88,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2-Rand */
    MPI_Isend(buffer[LY*LZ],                     1, slice_X_gath_type, g_nb_x_dn, 89,
	     g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND+2*LX*LY*LZ], 1, slice_X_cont_type, g_nb_x_up, 89, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2-Rand */
    MPI_Isend(buffer[(LX-2)*LY*LZ],                        1, slice_X_gath_type, g_nb_x_up, 90,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND+2*LX*LY*LZ + T*LY*LZ], 1, slice_X_cont_type, g_nb_x_dn, 90,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
#    endif
  MPI_Waitall(cntr, request, status);
  cntr=0;
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* The edges */

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  /* is on the x-Rand: xt-edge */
  MPI_Isend(buffer[(T+2)*LX*LY*LZ], 1, edge_XT_gath_type, g_nb_t_dn, 100,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND],  1, edge_XT_cont_type, g_nb_t_up, 100, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  /* xt-edge */
  MPI_Isend(buffer[(T+2)*LX*LY*LZ + (T-1)*LY*LZ], 1, edge_XT_gath_type, g_nb_t_up, 101,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 2*LY*LZ],      1, edge_XT_cont_type, g_nb_t_dn, 101,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* t2x-edge */
    MPI_Isend(buffer[(T+2)*LX*LY*LZ + LY*LZ],
	      1, edge_XT_gath_type, g_nb_t_dn, 102,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND],
	      1, edge_XT_cont_type, g_nb_t_up, 102, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* t2x-edge */
    MPI_Isend(buffer[(T+2)*LX*LY*LZ + (T-2)*LY*LZ],
	      1, edge_XT_gath_type, g_nb_t_up, 103,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 2*LY*LZ],
	      1, edge_XT_cont_type, g_nb_t_dn, 103,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* x2t-edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ],
	      1, edge_XT_gath_type, g_nb_t_dn, 104,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 4*LY*LZ],
	      1, edge_XT_cont_type, g_nb_t_up, 104, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* x2t-edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + (T-1)*LY*LZ],
	      1, edge_XT_gath_type, g_nb_t_up, 105,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 6*LY*LZ],
	      1, edge_XT_cont_type, g_nb_t_dn, 105,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
  /* end of if defined PARALLELXT || PARALLELXYT || PARALLELXYZT*/
#    endif

#    if (defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend(buffer[0],                            1, slice_Y_gath_type, g_nb_y_dn, 106,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + 2*LZ*(LX*LY + T*LY)], 1, slice_Y_cont_type, g_nb_y_up, 106, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  MPI_Isend(buffer[(LY-1)*LZ],                              1, slice_Y_gath_type, g_nb_y_up, 107,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + T*LX*LZ], 1, slice_Y_cont_type, g_nb_y_dn, 107,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* y2-Rand */
    MPI_Isend(buffer[LZ],                              1, slice_Y_gath_type, g_nb_y_dn, 108,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ], 1, slice_Y_cont_type, g_nb_y_up, 108, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* y2-Rand */
    MPI_Isend(buffer[(LY-2)*LZ],                                 1, slice_Y_gath_type, g_nb_y_up, 109,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + T*LX*LZ], 1, slice_Y_cont_type, g_nb_y_dn, 109,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
#    endif
  MPI_Waitall(cntr, request, status);
  cntr=0;
#    if (defined PARALLELXYT || defined PARALLELXYZT)

  /* jetzt wirds richtig eklig ... */

  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the y-Rand -> yx-edge*/
  MPI_Isend(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ], 1, edge_XY_gath_type, g_nb_x_dn, 110,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 4*LY*LZ],      1, edge_XY_cont_type, g_nb_x_up, 110, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* yx-edge */
  MPI_Isend(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + (LX-1)*LZ], 1, edge_XY_gath_type, g_nb_x_up, 111,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 4*LY*LZ + 2*T*LZ],         1, edge_XY_cont_type, g_nb_x_dn, 111,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the t-Rand -> ty-edge*/
  MPI_Isend(buffer[VOLUME],                           1, edge_YT_gath_type, g_nb_y_dn, 112,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ], 1, edge_YT_cont_type, g_nb_y_up, 112, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* ty-edge */
  MPI_Isend(buffer[VOLUME + (LY-1)*LZ],                         1, edge_YT_gath_type, g_nb_y_up, 113,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 2*LX*LZ], 1, edge_YT_cont_type, g_nb_y_dn, 113,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;



  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* x2y edge */
    MPI_Isend(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + LZ],
	      1, edge_XY_gath_type, g_nb_x_dn, 114,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ],
	      1, edge_XY_cont_type, g_nb_x_up, 114, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* x2y-edge */
    MPI_Isend(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-2)*LZ],
	      1, edge_XY_gath_type, g_nb_x_up, 115,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 2*T*LZ],
	      1, edge_XY_cont_type, g_nb_x_dn, 115,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;


    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2x -edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ],
	      1, edge_XY_gath_type, g_nb_x_dn, 116,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 4*T*LZ],
	      1, edge_XY_cont_type, g_nb_x_up, 116, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2x edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-1)*LZ],
	      1, edge_XY_gath_type, g_nb_x_up, 117,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 6*T*LZ],
	      1, edge_XY_cont_type, g_nb_x_dn, 117,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* t2y-edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND],
	      1, edge_YT_gath_type, g_nb_y_dn, 118,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ],
	      1, edge_YT_cont_type, g_nb_y_up, 118, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* t2y edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND + (LY-1)*LZ], 
	      1, edge_YT_gath_type, g_nb_y_up, 119,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 2*LX*LZ],
	      1, edge_YT_cont_type, g_nb_y_dn, 119,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2t edge */
    MPI_Isend(buffer[VOLUME + LZ],
	      1, edge_YT_gath_type, g_nb_y_dn, 120,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 4*LX*LZ],
	      1, edge_YT_cont_type, g_nb_y_up, 120, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2t-edge */
    MPI_Isend(buffer[VOLUME + (LY-2)*LZ],
	      1, edge_YT_gath_type, g_nb_y_up, 121,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 6*LX*LZ],
	      1, edge_YT_cont_type, g_nb_y_dn, 121,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }

  /* end of if defined PARALLELXYT || PARALLELXYZT */
#    endif
#    if defined PARALLELXYZT
  /* z-Rand */
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend(buffer[0],
	    1, slice_Z_gath_type, g_nb_z_dn, 122,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*LZ*T*LX], 
	    1, slice_Z_cont_type, g_nb_z_up, 122, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Isend(buffer[LZ-1],
	    1, slice_Z_gath_type, g_nb_z_up, 123,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + T*LX*LY],
	    1, slice_Z_cont_type, g_nb_z_dn, 123,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2-Rand */
    MPI_Isend(buffer[1],
	      1, slice_Z_gath_type, g_nb_z_dn, 124,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ],
	      1, slice_Z_cont_type, g_nb_z_up, 124, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2-Rand */
    MPI_Isend(buffer[LZ-2],
	      1, slice_Z_gath_type, g_nb_z_up, 125,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ + T*LX*LY],
	      1, slice_Z_cont_type, g_nb_z_dn, 125,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
#    endif
  MPI_Waitall(cntr, request, status);
#    if defined PARALLELXYZT
  cntr=0;
  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the z-Rand -> zx-edge*/
  MPI_Isend(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ], 
	    1, edge_XZ_gath_type, g_nb_x_dn, 126,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ],
	    1, edge_XZ_cont_type, g_nb_x_up, 126, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* zx-edge */
  MPI_Isend(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LX-1)*LY], 
	    1, edge_XZ_gath_type, g_nb_x_up, 127,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 2*T*LY],
	    1, edge_XZ_cont_type, g_nb_x_dn, 127,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  /* is on the t-Rand -> tz-edge*/
  MPI_Isend(buffer[VOLUME],
	    1, edge_ZT_gath_type, g_nb_z_dn, 128,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY], 
	    1, edge_ZT_cont_type, g_nb_z_up, 128, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;
  
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  /* tz-edge */
  MPI_Isend(buffer[VOLUME + (LZ-1)],
	    1, edge_ZT_gath_type, g_nb_z_up, 129,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY], 
	    1, edge_ZT_cont_type, g_nb_z_dn, 129,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the z-Rand -> zy-edge*/
  MPI_Isend(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ],
	    1, edge_YZ_gath_type, g_nb_y_dn, 130,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY], 
	    1, edge_YZ_cont_type, g_nb_y_up, 130, 
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* zy-edge */
  MPI_Isend(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LY-1)],
	    1, edge_YZ_gath_type, g_nb_y_up, 131,
	    g_cart_grid, &request[cntr]);
  MPI_Irecv(buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 2*T*LX], 
	    1, edge_YZ_cont_type, g_nb_y_dn, 131,
	    g_cart_grid, &request[cntr+1]);
  cntr=cntr+2;

  /* rectangular gauge action Stuff! */
  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* t2z edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND],
	      1, edge_ZT_gath_type, g_nb_z_dn, 132,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ],
	      1, edge_ZT_cont_type, g_nb_z_up, 132, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* t2z-edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND + (LZ-1)],
	      1, edge_ZT_gath_type, g_nb_z_up, 133,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 2*LX*LY],
	      1, edge_ZT_cont_type, g_nb_z_dn, 133,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;


    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2t -edge */
    MPI_Isend(buffer[VOLUME + 1],
	      1, edge_ZT_gath_type, g_nb_z_dn, 134,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 4*LX*LY],
	      1, edge_ZT_cont_type, g_nb_z_up, 134, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2t edge */
    MPI_Isend(buffer[VOLUME + (LZ-2)],
	      1, edge_ZT_gath_type, g_nb_z_up, 135,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 6*LX*LY],
	      1, edge_ZT_cont_type, g_nb_z_dn, 135,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* z2x-edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ],
	      1, edge_XZ_gath_type, g_nb_x_dn, 136,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY],
	      1, edge_XZ_cont_type, g_nb_x_up, 136, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* z2x edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LX-1)*LY], 
	      1, edge_XZ_gath_type, g_nb_x_up, 137,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 2*T*LY],
	      1, edge_XZ_cont_type, g_nb_x_dn, 137,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2z edge */
    MPI_Isend(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + LY],
	      1, edge_XZ_gath_type, g_nb_x_dn, 138,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 4*T*LY],
	      1, edge_XZ_cont_type, g_nb_x_up, 138, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2z-edge */
    MPI_Isend(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LX-2)*LY],
	      1, edge_XZ_gath_type, g_nb_x_up, 139,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 6*T*LY],
	      1, edge_XZ_cont_type, g_nb_x_dn, 139,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* z2y-edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ],
	      1, edge_YZ_gath_type, g_nb_y_dn, 140,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY],
	      1, edge_YZ_cont_type, g_nb_y_up, 140, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* z2y edge */
    MPI_Isend(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LY-1)], 
	      1, edge_YZ_gath_type, g_nb_y_up, 141,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 2*T*LX],
	      1, edge_YZ_cont_type, g_nb_y_dn, 141,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2z edge */
    MPI_Isend(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + 1],
	      1, edge_YZ_gath_type, g_nb_y_dn, 142,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 4*T*LX],
	      1, edge_YZ_cont_type, g_nb_y_up, 142, 
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2z-edge */
    MPI_Isend(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LY-2)],
	      1, edge_YZ_gath_type, g_nb_y_up, 143,
	      g_cart_grid, &request[cntr]);
    MPI_Irecv(buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 6*T*LX],
	      1, edge_YZ_cont_type, g_nb_y_dn, 143,
	      g_cart_grid, &request[cntr+1]);
    cntr=cntr+2;
  }
  MPI_Waitall(cntr, request, status);
#endif