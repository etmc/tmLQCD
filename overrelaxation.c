/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasenbusch
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
/*---------------------------------------------------------------------- 
! 
! Author
! M. Hasenbusch 2001
! Martin.Hasenbusch@desy.de
! This file provides the fuctions heatbath_sweep and overrel_sweep
! to update the gauge field.
! In the main program one merely says
!
!    .
!    .
!    heatbath_sweep();
!    overrel_sweep();
!     
! without any arguments. The random number generator
! is assumed to be initialized, e.g. by a call to rcinit from
! the main program. Similarly, the geometry has to be defined
! by ``call geometry". A sweep through the lattice proceeds
! in sequential order in a given time slice. 
! The updating procedure uses three Cabibbo-Marinari 
! subgroups for both, the over relaxation and the heatbath.
! For the latter we employ the procedure by Fabricius and  Haan.
! Details and references can be found in the notes by Peter Weisz.
!
! The code is based on the F code provide by 
! Stefan Sint   15/8/95 and Stefano Capitani   -   Jan/Feb 1997
! New:
! fuctions heatbath_sweep_adj and overrel_sweep_adj
! to update the gauge field with a mixed fundamental and adjoint
! action.
!
!   S_G  =  -\beta_f  \sum_P 1/N   Re Tr_f U_P
!           -\beta_a  \sum_P 1/N^2 ( Tr_f U_P^* )( Tr_f U_P )
!
!   for simplicity we have fixed betap = 6.0 here.
!
!  new Wed Oct  1 11:01:48 MEST 2003:    1-dim parallelisation for the 
!  gauge-update (pure Wilson only) by M.Hasenbusch
!  in the present version, it is not expected that the auxiliary fields 
!  for the boundaries are set consistently before the update is called.
!  However, after the update, the auxiliary fields are not at their proper
!  values, and xchange_gaugefield(); has to be called before e.g. D_psi();
!  can be used.
*
*
*  Checking that 1-dim. parallelisation (x-direction) works, 
*   (Check and correction in 
*     bin/pure_gauge.c, bin/local_update.c, bin/geometry.c, 
*     message-passing/xchange_gaugefield.c, observable/plaquette.c)
*    done by Kei-ichi Nagai  
*
*/
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "monomial/moment_energy.h"
#include "ranlxd.h"
#include "sse.h"
#include "get_staples.h"
#include "overrelaxation.h"

/****************************************************************/
/*
  flip_subgroup
  input: int ix, int mu, su3 vv, int i
*/
/****************************************************************/

void flip_subgroup(int ix, int mu, su3 vv, int i){
  static double vv0,vv1,vv2,vv3,aa0,aa1,aa2,aa3;
  static double aux,norm_vv_sq; 
  
  static su3 a,w,v;
  su3 *z;
  _su3_assign(v,vv);
  _su3_one(a);
  z=&g_gauge_field[ix][mu];     
  _su3_times_su3d(w,*z,v);

  /*
    According to Peter's notes ``A Cabibbo-Marinari SU(3)....", eqs. (A.14-A.17)
    we have */
  if(i==1)
  {
    vv0 =  creal(w.c00) + creal(w.c11);   
    vv3 = -cimag(w.c00) + cimag(w.c11);
    vv1 = -cimag(w.c01) - cimag(w.c10);
    vv2 = -creal(w.c01) + creal(w.c10);
  }
  else if(i==2)
  {
    vv0 =  creal(w.c00) + creal(w.c22);   
    vv3 = -cimag(w.c00) + cimag(w.c22);
    vv1 = -cimag(w.c02) - cimag(w.c20);
    vv2 = -creal(w.c02) + creal(w.c20);
  }
  else
  {
    vv0 =  creal(w.c11) + creal(w.c22);   
    vv3 = -cimag(w.c11) + cimag(w.c22);
    vv1 = -cimag(w.c12) - cimag(w.c21);
    vv2 = -creal(w.c12) + creal(w.c21);
  }

  norm_vv_sq= vv0 * vv0 + vv1 * vv1 + vv2 * vv2 + vv3 * vv3;

  aux= 2.0 * vv0 / norm_vv_sq;
  aa0 = aux * vv0-1.0;
  aa1 = aux * vv1;
  aa2 = aux * vv2;
  aa3 = aux * vv3;

  /*  aa is embedded in the SU(3) matrix (a) which can be multiplied on
      the link variable using the su3_type operator * . */
      
  if(i==1)
  {
    a.c00 = aa0 + aa3 * I;
    a.c11 = conj(a.c00);
    a.c01 = aa2 + aa1 * I;
    a.c10 = -conj(a.c01);
  }
  else if(i==2)
  {
    a.c00 = aa0 + aa3 * I;
    a.c22 = conj(a.c00);
    a.c02 = aa2 + aa1 * I;
    a.c20 = -conj(a.c02);
  }
  else
  {
    a.c11 = aa0 + aa3 * I;
    a.c22 = conj(a.c11);
    a.c12 = aa2 + aa1 * I;
    a.c21 = -conj(a.c12);
  }

  _su3_times_su3(w,a,*z);
  *z=w;
}

#if defined PARALLEL1
void overrel_sweep(){
  int x0,x1,x2,x3;
  int mu,ix;
  static su3 v;
  if(LX<2) {printf("LX is smaller than 2 \n"); exit(0);}
/* xchange the gauge-field */
  xchange_gaugefield(g_gauge_field);
/* update the left half of the sublattice */
  for(x1=0;x1<LX/2;x1++)
    {
    for(x0=0;x0<T;x0++)
      {
      for(x2=0;x2<LY;x2++)
        {
        for(x3=0;x3<LZ;x3++)
          {
          ix=g_ipt[x0][x1][x2][x3];
          for(mu=0;mu<4;mu++){
            get_staples(&v,ix,mu,g_gauge_field);
            flip_subgroup(ix,mu,v,1);
            flip_subgroup(ix,mu,v,2);
            flip_subgroup(ix,mu,v,3);
            }
          }
        }
      }
    }
/* xchange the gauge-field */
  xchange_gaugefield(g_gauge_field);
/* update the right half of the sub lattice */
  for(x1=LX/2;x1<LX;x1++)
    {
    for(x0=0;x0<T;x0++)
      {
      for(x2=0;x2<LY;x2++)
        {
        for(x3=0;x3<LZ;x3++)
          {
          ix=g_ipt[x0][x1][x2][x3];
          for(mu=0;mu<4;mu++){
            get_staples(&v,ix,mu,g_gauge_field);
            flip_subgroup(ix,mu,v,1);
            flip_subgroup(ix,mu,v,2);
            flip_subgroup(ix,mu,v,3);
            }
          }
        }
      }
    }
}
#else
void overrel_sweep(){
  int mu,ix;
  static su3 v;
  for(mu=0;mu<4;mu++){
    for(ix=0;ix<VOLUME;ix++){
      get_staples(&v,ix,mu,g_gauge_field);
      flip_subgroup(ix,mu,v,1);
      flip_subgroup(ix,mu,v,2);
      flip_subgroup(ix,mu,v,3);
    }
  }
}
#endif

