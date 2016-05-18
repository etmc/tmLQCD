/***********************************************************************
 *
 * Copyright (C) 2008 Albert Deuzeman, Siebren Reker, Carsten Urbach
 *               2010 Claude Tadonki, Carsten Urbach
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

void _PSWITCH(little_field_gather)(_C_TYPE * w) {
  int ib;
  _C_TYPE * wt_buf=NULL, * wx_buf=NULL, * wy_buf=NULL, * wz_buf=NULL, 
    * w_buf=NULL, * w_source=NULL, * w_dest=NULL;
  _C_TYPE * wt=NULL, * wx=NULL, * wy=NULL, * wz=NULL;

  /************************************************************************/
  /* This routine has been extended for multi_dimensional blocking        */
  /* by Claude Tadonki (claude.tadonki@u-psud.fr) from PetaQCD project    */
  /* June 2010                                                            */
  /************************************************************************/

  wt = w + ( 0*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction t starts
  wx = w + ( 1*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction x starts
  wy = w + ( 2*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction y starts
  wz = w + ( 3*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction z starts

#ifdef MPI
  int request = 0;
  int err;
  w_buf = calloc(8 * nb_blocks * g_N_s, sizeof(_C_TYPE)); // +-t +-x +-y +-z

  wt_buf = w_buf + ( 0*(2*nb_blocks)) * g_N_s; // Were data in the direction t starts
  wx_buf = w_buf + ( 1*(2*nb_blocks)) * g_N_s; // Were data in the direction x starts
  wy_buf = w_buf + ( 2*(2*nb_blocks)) * g_N_s; // Were data in the direction y starts
  wz_buf = w_buf + ( 3*(2*nb_blocks)) * g_N_s; // Were data in the direction z starts

  /* We first exchange the fields regardless of block considerations                   */
  /* The data need to be received in an intermediate buffer because of later shuffling */

  if(g_nproc_t > 1) {
    /* Send t up */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_t_up, T_UP, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wt_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_t_dn, T_UP, g_cart_grid, &lrequests[request]);
    request++;
    /* Send t down */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_t_dn, T_DN, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wt_buf, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_t_up, T_DN, g_cart_grid, &lrequests[request]);
    request++;
  }
  if(g_nproc_x > 1) {
    /* Send x up */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_x_up, X_UP, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wx_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_x_dn, X_UP, g_cart_grid, &lrequests[request]);
    request++;
    /* Send x down */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_x_dn, X_DN, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wx_buf, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_x_up, X_DN, g_cart_grid, &lrequests[request]);
    request++;
  }
  if(g_nproc_y > 1) {
    /* Send y up */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_y_up, Y_UP, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wy_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_y_dn, Y_UP, g_cart_grid, &lrequests[request]);
    request++;
    /* Send y down */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_y_dn, Y_DN, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wy_buf, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_y_up, Y_DN, g_cart_grid, &lrequests[request]);
    request++;
  }
  if(g_nproc_z > 1) {
    /* Send z up */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_z_up, Z_UP, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wz_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_z_dn, Z_UP, g_cart_grid, &lrequests[request]);
    request++;
    /* Send z down */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_z_dn, Z_DN, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wz_buf, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_z_up, Z_DN, g_cart_grid, &lrequests[request]);
    request++;
  }
  err = MPI_Waitall(request, lrequests, lstatus);

#endif
  
  /* We now correct the field according to block partitionning               */
  /* We could have avoid the previous corresponding MPI communication        */
  /* We proceed like this for code simplicity, maybe will be optimized later */
  
  for(int pm = 0; pm < 8; pm++) {
    for(int bt = 0; bt < nblks_t; bt++) {
      for(int bx = 0; bx < nblks_x; bx++) {
        for(int by = 0; by < nblks_y; by++) {
          for(int bz = 0; bz < nblks_z; bz++) {
            ib = block_index(bt, bx, by, bz) * g_N_s;
            switch(pm){ 
            case T_UP: /* Direction +t */
              w_dest = wt + ib;
              if( bt == nblks_t - 1 ) {
                ib = block_index(0, bx, by, bz) * g_N_s; 
                if(g_nproc_t > 1) w_source = wt_buf + ib;
                else w_source = w + ib;
              }
              // got it from the MPI exchange
              else  {
                ib = block_index(bt + 1, bx, by, bz) * g_N_s; 
                w_source = w + ib;
              }
              // got it from the diagonal block
              break; 
            case T_DN: /* Direction -t */
              w_dest = wt + ib + nb_blocks * g_N_s;
              if( bt == 0 ) {
                ib = block_index(nblks_t - 1, bx, by, bz) * g_N_s; 
                if(g_nproc_t > 1) w_source = wt_buf + ib + nb_blocks * g_N_s;
                else w_source = w + ib;
              }
              // got it from the MPI exchange
              else  {
                ib = block_index(bt - 1, bx, by, bz) * g_N_s;
                w_source = w + ib;
              }
              // got it from the diagonal block
              break; 
            case X_UP: /* Direction +x */
              w_dest = wx + ib;
              if( bx == nblks_x - 1 ) {
                ib = block_index(bt, 0, by, bz) * g_N_s; 
                if(g_nproc_x > 1) w_source = wx_buf + ib;
                else w_source = w + ib;
              }
              // got it from the MPI exchange
              else  {
                ib = block_index(bt, bx + 1, by, bz) * g_N_s; 
                w_source = w + ib;
              }
              // got it from the diagonal block
              break; 
            case X_DN: /* Direction -x */
              w_dest = wx + ib + nb_blocks * g_N_s;
              if( bx == 0 ) {
                ib = block_index(bt, nblks_x - 1, by, bz) * g_N_s; 
                if(g_nproc_x > 1) w_source = wx_buf + ib + nb_blocks * g_N_s;
                else w_source = w + ib;
              }
              // got it from the MPI exchange
              else  {
                ib = block_index(bt, bx - 1, by, bz) * g_N_s;
                w_source = w + ib;
              }
              // got it from the diagonal block
              break; 
            case Y_UP: /* Direction +y */
              w_dest = wy + ib;
              if( by == nblks_y - 1 ) {
                ib = block_index(bt, bx, 0, bz) * g_N_s; 
                if(g_nproc_y > 1) w_source = wy_buf + ib;
                else w_source = w + ib;
              }
              // got it from the MPI exchange
              else  {
                ib = block_index(bt, bx, by + 1, bz) * g_N_s; 
                w_source = w + ib;
              }
              // got it from the diagonal block
              break; 
            case Y_DN: /* Direction -y */
              w_dest = wy + ib + nb_blocks * g_N_s;
              if( by == 0 ) {
                ib = block_index(bt, bx, nblks_y - 1, bz) * g_N_s; 
                if(g_nproc_y > 1) w_source = wy_buf + ib + nb_blocks * g_N_s;
                else w_source = w + ib;
              }
              // got it from the MPI exchange
              else  {
                ib = block_index(bt, bx, by - 1, bz) * g_N_s;
                w_source = w + ib;
              }
              // got it from the diagonal block
              break; 
            case Z_UP: /* Direction +z */
              w_dest = wz + ib;
              if( bz == nblks_z - 1 ) {
                ib = block_index(bt, bx, by, 0) * g_N_s; 
                if(g_nproc_z > 1) w_source = wz_buf + ib;
                else w_source = w + ib;
              }
              // got it from the MPI exchange
              else  {
                ib = block_index(bt, bx, by, bz + 1) * g_N_s; 
                w_source = w + ib;
              }
              // got it from the diagonal block
              break; 
            case Z_DN: /* Direction -z */
              w_dest = wz + ib + nb_blocks * g_N_s;
              if( bz == 0 ) {
                ib = block_index(bt, bx, by, nblks_z - 1) * g_N_s; 
                if(g_nproc_z > 1) w_source = wz_buf + ib + nb_blocks * g_N_s;
                else w_source = w + ib;
              }
              // got it from the MPI exchange
              else  {
                ib = block_index(bt, bx, by, bz - 1) * g_N_s; 
                w_source = w + ib; 
              }
              // got it from the diagonal block
              break; 
              
            default: 
              w_dest = NULL;
              w_source = NULL;
            }
            memcpy(w_dest, w_source, g_N_s * sizeof(_C_TYPE));
          }
        }
      }
    }
  }
  free(w_buf);
  
  return;
}

void _PSWITCH(little_field_gather_eo)(const int eo, _C_TYPE * w) {

  int ib, ib2;
  _C_TYPE *wt = NULL, *wx = NULL, *wy = NULL, *wz = NULL;
  _C_TYPE *wt_buf = NULL, *wx_buf = NULL, *wy_buf = NULL, *wz_buf = NULL, *w_buf = NULL, *w_source = NULL, *w_dest = NULL;

  
  wt = w + ( 0*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction t starts
  wx = w + ( 1*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction x starts
  wy = w + ( 2*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction y starts
  wz = w + ( 3*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction z starts

#ifdef MPI
  int request = 0;
  int err;
  w_buf = calloc(8 * nb_blocks * g_N_s, sizeof(_C_TYPE)); // +-t +-x +-y +-z

  wt_buf = w_buf + ( 0*(2*nb_blocks)) * g_N_s; // Were data in the direction t starts
  wx_buf = w_buf + ( 1*(2*nb_blocks)) * g_N_s; // Were data in the direction x starts
  wy_buf = w_buf + ( 2*(2*nb_blocks)) * g_N_s; // Were data in the direction y starts
  wz_buf = w_buf + ( 3*(2*nb_blocks)) * g_N_s; // Were data in the direction z starts

  /* We first exchange the fields regardless of block considerations                   */
  /* The data need to be received in an intermediate buffer because of later shuffling */

  if(g_nproc_t > 1) {
    /* Send t up */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_t_up, T_UP, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wt_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_t_dn, T_UP, g_cart_grid, &lrequests[request]);
    request++;
    /* Send t down */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_t_dn, T_DN, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wt_buf, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_t_up, T_DN, g_cart_grid, &lrequests[request]);
    request++;
  }
  if(g_nproc_x > 1) {
    /* Send x up */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_x_up, X_UP, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wx_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_x_dn, X_UP, g_cart_grid, &lrequests[request]);
    request++;
    /* Send x down */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_x_dn, X_DN, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wx_buf, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_x_up, X_DN, g_cart_grid, &lrequests[request]);
    request++;
  }
  if(g_nproc_y > 1) {
    /* Send y up */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_y_up, Y_UP, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wy_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_y_dn, Y_UP, g_cart_grid, &lrequests[request]);
    request++;
    /* Send y down */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_y_dn, Y_DN, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wy_buf, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_y_up, Y_DN, g_cart_grid, &lrequests[request]);
    request++;
  }
  if(g_nproc_z > 1) {
    /* Send z up */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_z_up, Z_UP, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wz_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_z_dn, Z_UP, g_cart_grid, &lrequests[request]);
    request++;
    /* Send z down */
    MPI_Isend(w, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_z_dn, Z_DN, g_cart_grid, &lrequests[request]);
    request++;
    MPI_Irecv(wz_buf, nb_blocks * g_N_s, _MPI_C_TYPE, g_nb_z_up, Z_DN, g_cart_grid, &lrequests[request]);
    request++;
  }
  err = MPI_Waitall(request, lrequests, lstatus);

#endif
  
  /* We now correct the field according to block partitionning               */

  for(int pm = 0; pm < 8; pm++) {
    ib2=0;
    for(int bt = 0; bt < nblks_t; bt++) {
      for(int bx = 0; bx < nblks_x; bx++) {
        for(int by = 0; by < nblks_y; by++) {
          for(int bz = 0; bz < nblks_z; bz++) {
            if ((bt+bx+by+bz)%2==eo) {
              ib2 = index_block_eo[block_index(bt, bx, by, bz)] * g_N_s;
              
              switch(pm){ 
              case T_UP: /* Direction +t */
                w_dest = wt + ib2;
                if( bt == nblks_t - 1 ) {
                  ib = index_block_eo[block_index(0,bx, by,bz)] * g_N_s; 
                  if(g_nproc_t > 1) w_source = wt_buf + ib;     
                  else w_source = w + ib;
                }
                // got it from the MPI exchange
                else  {
                  ib = index_block_eo[block_index(bt+1, bx, by, bz)] * g_N_s; 
                  w_source = w + ib;    
                }
                // got it from the diagonal block
                break; 
              case T_DN: /* Direction -t */
                w_dest = wt + ib2 + nb_blocks * g_N_s;
                if( bt == 0) {
                  ib = index_block_eo[block_index(nblks_t-1, bx,by,bz)] * g_N_s; 
                  if(g_nproc_t > 1) w_source = wt_buf + ib + nb_blocks * g_N_s;
                  else w_source = w + ib;
                } // got it from the MPI exchange
                else  {
                  ib = index_block_eo[block_index(bt-1,bx, by, bz)] * g_N_s; 
                  w_source = w + ib; 
                }
                // got it from the diagonal block
                break; 
              case X_UP: /* Direction +x */
                w_dest = wx + ib2;
                if( bx == nblks_x - 1 ) {
                  ib = index_block_eo[block_index(bt, 0, by,bz)] * g_N_s; 
                  if(g_nproc_x > 1) w_source = wx_buf + ib;     
                  else w_source = w + ib;
                }
                // got it from the MPI exchange
                else  {
                  ib = index_block_eo[block_index(bt, bx+1, by, bz)] * g_N_s; 
                  w_source = w + ib;
                }
                // got it from the diagonal block
                break; 
              case X_DN: /* Direction -x */
                w_dest = wx + ib2 + nb_blocks * g_N_s;
                if( bx == 0) {ib = index_block_eo[block_index(bt, nblks_x-1, by,bz)] * g_N_s;
                  if(g_nproc_x > 1) w_source = wx_buf + ib + nb_blocks * g_N_s;
                  else w_source = w + ib;
                }
                // got it from the MPI exchange
                else  {
                  ib = index_block_eo[block_index(bt, bx-1, by, bz)] * g_N_s; 
                  w_source = w + ib;
                }
                // got it from the diagonal block
                break; 
              case Y_UP: /* Direction +y */
                w_dest = wy + ib2;
                if( by == nblks_y - 1 ) {
                  ib = index_block_eo[block_index(bt, bx, 0,bz)] * g_N_s; 
                  if(g_nproc_y > 1) w_source = wy_buf + ib;
                  else w_source = w + ib;
                }
                // got it from the MPI exchange
                else  {
                  ib = index_block_eo[block_index(bt, bx, by+1, bz)] * g_N_s;
                  w_source = w + ib;
                }
                // got it from the diagonal block
                break; 
              case Y_DN: /* Direction -y */
                w_dest = wy + ib2 + nb_blocks * g_N_s;
                if( by == 0) {
                  ib = index_block_eo[block_index(bt, bx, nblks_y-1, bz)] * g_N_s;
                  if(g_nproc_y > 1) w_source = wy_buf + ib + nb_blocks * g_N_s;
                  else w_source = w + ib;
                }
                // got it from the MPI exchange
                else  {
                  ib = index_block_eo[block_index(bt, bx, by-1, bz)] * g_N_s;
                  w_source = w + ib; 
                }
                // got it from the diagonal block
                break; 
              case Z_UP: /* Direction +z */
                w_dest = wz + ib2;
                if( bz == nblks_z - 1 ) {
                  ib = index_block_eo[block_index(bt, bx, by, 0)] * g_N_s;
                  if(g_nproc_z > 1) w_source = wz_buf + ib;
                  else w_source = w + ib;
                }
                // got it from the MPI exchange
                else  {
                  ib = index_block_eo[block_index(bt, bx, by, bz + 1)] * g_N_s;
                  w_source = w + ib;
                }
                // got it from the diagonal block
                break; 
              case Z_DN: /* Direction -z */
                w_dest = wz + ib2 + nb_blocks * g_N_s;
                if( bz == 0) {
                  ib = index_block_eo[block_index(bt, bx, by, nblks_z - 1)] * g_N_s;
                  if(g_nproc_z > 1) w_source = wz_buf + ib + nb_blocks * g_N_s;
                  else w_source = w + ib;
                }
                // got it from the MPI exchange
                else  {
                  ib = index_block_eo[block_index(bt, bx, by, bz - 1)] * g_N_s;
                  w_source = w + ib;
                }
                // got it from the diagonal block
                break; 
              default:
                w_dest = NULL;
                w_source = NULL;
              }
              memcpy(w_dest, w_source, g_N_s * sizeof(_C_TYPE));
            }
          }
        }
      }
    }
  }
  free(w_buf);
  return;
}
