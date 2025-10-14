  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv(buffer[0],          1, slice_T_cont_type, g_nb_t_dn, 83, 
	       buffer[VOLUME], 1, slice_T_cont_type, g_nb_t_up, 83, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv(buffer[(T-1)*LX*LY*LZ], 1, slice_T_cont_type, g_nb_t_up, 84, 
	       buffer[(T+1)*LX*LY*LZ], 1, slice_T_cont_type, g_nb_t_dn, 84, 
	       g_cart_grid, &status);

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    /* t2-Rand */
    MPI_Sendrecv(buffer[1*LX*LY*LZ],     1, slice_T_cont_type, g_nb_t_dn, 85, 
		 buffer[VOLUMEPLUSRAND], 1, slice_T_cont_type, g_nb_t_up, 85, 
		 g_cart_grid, &status);

    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    /* t2-Rand */
    MPI_Sendrecv(buffer[(T-2)*LX*LY*LZ],          1, slice_T_cont_type, g_nb_t_up, 86, 
		 buffer[VOLUMEPLUSRAND+LX*LY*LZ], 1, slice_T_cont_type, g_nb_t_dn, 86, 
		 g_cart_grid, &status);
  }
  
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv(buffer[0],              1, slice_X_gath_type, g_nb_x_dn, 93,
	       buffer[(T+2)*LX*LY*LZ], 1, slice_X_cont_type, g_nb_x_up, 93, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* x2-Rand */
  MPI_Sendrecv(buffer[(LX-1)*LY*LZ],             1, slice_X_gath_type, g_nb_x_up, 94,
	       buffer[(T+2)*LX*LY*LZ + T*LY*LZ], 1, slice_X_cont_type, g_nb_x_dn, 94,
	       g_cart_grid, &status);

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2-Rand */
    MPI_Sendrecv(buffer[LY*LZ],                     1, slice_X_gath_type, g_nb_x_dn, 95,
		 buffer[VOLUMEPLUSRAND+2*LX*LY*LZ], 1, slice_X_cont_type, g_nb_x_up, 95, 
		 g_cart_grid, &status);
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2-Rand */
    MPI_Sendrecv(buffer[(LX-2)*LY*LZ],                        1, slice_X_gath_type, g_nb_x_up, 96,
		 buffer[VOLUMEPLUSRAND+2*LX*LY*LZ + T*LY*LZ], 1, slice_X_cont_type, g_nb_x_dn, 96,
		 g_cart_grid, &status);
  }

  /* The edges */

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  /* is on the x-Rand: xt-edge */
  MPI_Sendrecv(buffer[(T+2)*LX*LY*LZ], 1, edge_XT_gath_type, g_nb_t_dn, 95,
	       buffer[VOLUME + RAND],  1, edge_XT_cont_type, g_nb_t_up, 95, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  /* xt-edge */
  MPI_Sendrecv(buffer[(T+2)*LX*LY*LZ + (T-1)*LY*LZ], 1, edge_XT_gath_type, g_nb_t_up, 96,
	       buffer[VOLUME + RAND + 2*LY*LZ],      1, edge_XT_cont_type, g_nb_t_dn, 96,
	       g_cart_grid, &status);

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* t2x-edge */
    MPI_Sendrecv(buffer[(T+2)*LX*LY*LZ + LY*LZ],
		 1, edge_XT_gath_type, g_nb_t_dn, 97,
		 buffer[VOLUMEPLUSRAND + RAND],
		 1, edge_XT_cont_type, g_nb_t_up, 97, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* t2x-edge */
    MPI_Sendrecv(buffer[(T+2)*LX*LY*LZ + (T-2)*LY*LZ],
		 1, edge_XT_gath_type, g_nb_t_up, 98,
		 buffer[VOLUMEPLUSRAND + RAND + 2*LY*LZ],
		 1, edge_XT_cont_type, g_nb_t_dn, 98,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* x2t-edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ],
		 1, edge_XT_gath_type, g_nb_t_dn, 97,
		 buffer[VOLUMEPLUSRAND + RAND + 4*LY*LZ],
		 1, edge_XT_cont_type, g_nb_t_up, 97, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* x2t-edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + (T-1)*LY*LZ],
		 1, edge_XT_gath_type, g_nb_t_up, 98,
		 buffer[VOLUMEPLUSRAND + RAND + 6*LY*LZ],
		 1, edge_XT_cont_type, g_nb_t_dn, 98,
		 g_cart_grid, &status);
  }
  /* end of if defined PARALLELXT || PARALLELXYT || PARALLELXYZT*/
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv(buffer[0],                            1, slice_Y_gath_type, g_nb_y_dn, 103,
	       buffer[VOLUME + 2*LZ*(LX*LY + T*LY)], 1, slice_Y_cont_type, g_nb_y_up, 103, 
	       g_cart_grid, &status);
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  MPI_Sendrecv(buffer[(LY-1)*LZ],                              1, slice_Y_gath_type, g_nb_y_up, 104,
	       buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + T*LX*LZ], 1, slice_Y_cont_type, g_nb_y_dn, 104,
	       g_cart_grid, &status);

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* y2-Rand */
    MPI_Sendrecv(buffer[LZ],                              1, slice_Y_gath_type, g_nb_y_dn, 105,
		 buffer[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ], 1, slice_Y_cont_type, g_nb_y_up, 105, 
		 g_cart_grid, &status);
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* y2-Rand */
    MPI_Sendrecv(buffer[(LY-2)*LZ],                                 1, slice_Y_gath_type, g_nb_y_up, 106,
		 buffer[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + T*LX*LZ], 1, slice_Y_cont_type, g_nb_y_dn, 106,
		 g_cart_grid, &status);
  }

  /* jetzt wirds richtig eklig ... */

  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the y-Rand -> yx-edge*/
  MPI_Sendrecv(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ], 1, edge_XY_gath_type, g_nb_x_dn, 107,
	       buffer[VOLUME + RAND + 4*LY*LZ],      1, edge_XY_cont_type, g_nb_x_up, 107, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* yx-edge */
  MPI_Sendrecv(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + (LX-1)*LZ], 1, edge_XY_gath_type, g_nb_x_up, 108,
	       buffer[VOLUME + RAND + 4*LY*LZ + 2*T*LZ],         1, edge_XY_cont_type, g_nb_x_dn, 108,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the t-Rand -> ty-edge*/
  MPI_Sendrecv(buffer[VOLUME],                           1, edge_YT_gath_type, g_nb_y_dn, 109,
	       buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ], 1, edge_YT_cont_type, g_nb_y_up, 109, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* ty-edge */
  MPI_Sendrecv(buffer[VOLUME + (LY-1)*LZ],                         1, edge_YT_gath_type, g_nb_y_up, 110,
	       buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 2*LX*LZ], 1, edge_YT_cont_type, g_nb_y_dn, 110,
	       g_cart_grid, &status);



  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* x2y edge */
    MPI_Sendrecv(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + LZ],
		 1, edge_XY_gath_type, g_nb_x_dn, 97,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ],
		 1, edge_XY_cont_type, g_nb_x_up, 97, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* x2y-edge */
    MPI_Sendrecv(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-2)*LZ],
		 1, edge_XY_gath_type, g_nb_x_up, 98,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 2*T*LZ],
		 1, edge_XY_cont_type, g_nb_x_dn, 98,
		 g_cart_grid, &status);


    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2x -edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ],
		 1, edge_XY_gath_type, g_nb_x_dn, 97,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 4*T*LZ],
		 1, edge_XY_cont_type, g_nb_x_up, 97, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2x edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-1)*LZ],
		 1, edge_XY_gath_type, g_nb_x_up, 98,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 6*T*LZ],
		 1, edge_XY_cont_type, g_nb_x_dn, 98,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* t2y-edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND],
		 1, edge_YT_gath_type, g_nb_y_dn, 197,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ],
		 1, edge_YT_cont_type, g_nb_y_up, 197, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* t2y edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND + (LY-1)*LZ], 
		 1, edge_YT_gath_type, g_nb_y_up, 198,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 2*LX*LZ],
		 1, edge_YT_cont_type, g_nb_y_dn, 198,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2t edge */
    MPI_Sendrecv(buffer[VOLUME + LZ],
		 1, edge_YT_gath_type, g_nb_y_dn, 297,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 4*LX*LZ],
		 1, edge_YT_cont_type, g_nb_y_up, 297, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2t-edge */
    MPI_Sendrecv(buffer[VOLUME + (LY-2)*LZ],
		 1, edge_YT_gath_type, g_nb_y_up, 298,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 6*LX*LZ],
		 1, edge_YT_cont_type, g_nb_y_dn, 298,
		 g_cart_grid, &status);
  }

  /* end of if defined PARALLELXYT || PARALLELXYZT */
#  endif
#  if defined PARALLELXYZT
  /* z-Rand */
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Sendrecv(buffer[0],
	       1, slice_Z_gath_type, g_nb_z_dn, 303,
	       buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*LZ*T*LX], 
	       1, slice_Z_cont_type, g_nb_z_up, 303, 
	       g_cart_grid, &status);
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Sendrecv(buffer[LZ-1],
	       1, slice_Z_gath_type, g_nb_z_up, 304,
	       buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + T*LX*LY],
	       1, slice_Z_cont_type, g_nb_z_dn, 304,
	       g_cart_grid, &status);

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2-Rand */
    MPI_Sendrecv(buffer[1],
		 1, slice_Z_gath_type, g_nb_z_dn, 305,
		 buffer[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ],
		 1, slice_Z_cont_type, g_nb_z_up, 305, 
		 g_cart_grid, &status);
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2-Rand */
    MPI_Sendrecv(buffer[LZ-2],
		 1, slice_Z_gath_type, g_nb_z_up, 306,
		 buffer[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ + T*LX*LY],
		 1, slice_Z_cont_type, g_nb_z_dn, 306,
		 g_cart_grid, &status);
  }

  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the z-Rand -> zx-edge*/
  MPI_Sendrecv(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ], 
	       1, edge_XZ_gath_type, g_nb_x_dn, 307,
	       buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ],
	       1, edge_XZ_cont_type, g_nb_x_up, 307, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* zx-edge */
  MPI_Sendrecv(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LX-1)*LY], 
	       1, edge_XZ_gath_type, g_nb_x_up, 308,
	       buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 2*T*LY],
	       1, edge_XZ_cont_type, g_nb_x_dn, 308,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  /* is on the t-Rand -> tz-edge*/
  MPI_Sendrecv(buffer[VOLUME],
	       1, edge_ZT_gath_type, g_nb_z_dn, 309,
	       buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY], 
	       1, edge_ZT_cont_type, g_nb_z_up, 309, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  /* tz-edge */
  MPI_Sendrecv(buffer[VOLUME + (LZ-1)],
	       1, edge_ZT_gath_type, g_nb_z_up, 310,
	       buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY], 
	       1, edge_ZT_cont_type, g_nb_z_dn, 310,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the z-Rand -> zy-edge*/
  MPI_Sendrecv(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ],
	       1, edge_YZ_gath_type, g_nb_y_dn, 310,
	       buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY], 
	       1, edge_YZ_cont_type, g_nb_y_up, 310, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* zy-edge */
  MPI_Sendrecv(buffer[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LY-1)],
	       1, edge_YZ_gath_type, g_nb_y_up, 310,
	       buffer[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 2*T*LX], 
	       1, edge_YZ_cont_type, g_nb_y_dn, 310,
	       g_cart_grid, &status);

  /* rectangular gauge action Stuff! */
  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* t2z edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND],
		 1, edge_ZT_gath_type, g_nb_z_dn, 500,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ],
		 1, edge_ZT_cont_type, g_nb_z_up, 500, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* t2z-edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND + (LZ-1)],
		 1, edge_ZT_gath_type, g_nb_z_up, 501,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 2*LX*LY],
		 1, edge_ZT_cont_type, g_nb_z_dn, 501,
		 g_cart_grid, &status);


    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2t -edge */
    MPI_Sendrecv(buffer[VOLUME + 1],
		 1, edge_ZT_gath_type, g_nb_z_dn, 502,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 4*LX*LY],
		 1, edge_ZT_cont_type, g_nb_z_up, 502, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2t edge */
    MPI_Sendrecv(buffer[VOLUME + (LZ-2)],
		 1, edge_ZT_gath_type, g_nb_z_up, 503,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 6*LX*LY],
		 1, edge_ZT_cont_type, g_nb_z_dn, 503,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* z2x-edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ],
		 1, edge_XZ_gath_type, g_nb_x_dn, 504,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY],
		 1, edge_XZ_cont_type, g_nb_x_up, 504, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* z2x edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LX-1)*LY], 
		 1, edge_XZ_gath_type, g_nb_x_up, 504,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 2*T*LY],
		 1, edge_XZ_cont_type, g_nb_x_dn, 504,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2z edge */
    MPI_Sendrecv(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + LY],
		 1, edge_XZ_gath_type, g_nb_x_dn, 505,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 4*T*LY],
		 1, edge_XZ_cont_type, g_nb_x_up, 505, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2z-edge */
    MPI_Sendrecv(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LX-2)*LY],
		 1, edge_XZ_gath_type, g_nb_x_up, 506,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 6*T*LY],
		 1, edge_XZ_cont_type, g_nb_x_dn, 506,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* z2y-edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ],
		 1, edge_YZ_gath_type, g_nb_y_dn, 507,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY],
		 1, edge_YZ_cont_type, g_nb_y_up, 507, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* z2y edge */
    MPI_Sendrecv(buffer[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LY-1)], 
		 1, edge_YZ_gath_type, g_nb_y_up, 508,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 2*T*LX],
		 1, edge_YZ_cont_type, g_nb_y_dn, 508,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2z edge */
    MPI_Sendrecv(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + 1],
		 1, edge_YZ_gath_type, g_nb_y_dn, 509,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 4*T*LX],
		 1, edge_YZ_cont_type, g_nb_y_up, 509, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2z-edge */
    MPI_Sendrecv(buffer[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LY-2)],
		 1, edge_YZ_gath_type, g_nb_y_up, 510,
		 buffer[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 6*T*LX],
		 1, edge_YZ_cont_type, g_nb_y_dn, 510,
		 g_cart_grid, &status);

  }

#endif /* PARALLELXYZT */
