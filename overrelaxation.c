/*---------------------------------------------------------------------- 
! $Id$
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
!           -\beta_a  \sum_P 1/N^2 ( Tr_f U_P^* ) ( Tr_f U_P )
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
#include "hybrid_update.h"
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
  if(i==1){
    vv0 =  w.c00.re+w.c11.re;   
    vv3 = -w.c00.im+w.c11.im;
    vv1 = -w.c01.im-w.c10.im;
    vv2 = -w.c01.re+w.c10.re;
  }
  else if(i==2){
    vv0 =  w.c00.re+w.c22.re;   
    vv3 = -w.c00.im+w.c22.im;
    vv1 = -w.c02.im-w.c20.im;
    vv2 = -w.c02.re+w.c20.re;
  }
  else{
    vv0 =  w.c11.re+w.c22.re;   
    vv3 = -w.c11.im+w.c22.im;
    vv1 = -w.c12.im-w.c21.im;
    vv2 = -w.c12.re+w.c21.re;
  }

  norm_vv_sq=vv0*vv0+vv1*vv1+vv2*vv2+vv3*vv3;

  aux= 2.0*vv0/norm_vv_sq;
  aa0 = aux*vv0-1.0;
  aa1 = aux*vv1;
  aa2 = aux*vv2;
  aa3 = aux*vv3;

  /*  aa is embedded in the SU(3) matrix (a) which can be multiplied on
      the link variable using the su3_type operator * . */
      
  if(i==1){
    a.c00.re =  aa0;
    a.c00.im =  aa3;
    a.c11.re =  aa0;
    a.c11.im = -aa3;
    a.c01.re =  aa2;
    a.c01.im =  aa1;
    a.c10.re = -aa2;
    a.c10.im =  aa1;
  }
  else if(i==2){
    a.c00.re =  aa0;
    a.c00.im =  aa3;
    a.c22.re =  aa0;
    a.c22.im = -aa3;
    a.c02.re =  aa2;
    a.c02.im =  aa1;
    a.c20.re = -aa2;
    a.c20.im =  aa1;
  }
  else{
    a.c11.re =  aa0;
    a.c11.im =  aa3;
    a.c22.re =  aa0;
    a.c22.im = -aa3;
    a.c12.re =  aa2;
    a.c12.im =  aa1;
    a.c21.re = -aa2;
    a.c21.im =  aa1;
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
  xchange_gaugefield();
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
            v=get_staples(ix,mu);
            flip_subgroup(ix,mu,v,1);
            flip_subgroup(ix,mu,v,2);
            flip_subgroup(ix,mu,v,3);
            }
          }
        }
      }
    }
/* xchange the gauge-field */
  xchange_gaugefield();
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
            v=get_staples(ix,mu);
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
      v=get_staples(ix,mu);
      flip_subgroup(ix,mu,v,1);
      flip_subgroup(ix,mu,v,2);
      flip_subgroup(ix,mu,v,3);
    }
  }
}
#endif

