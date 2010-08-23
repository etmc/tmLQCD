#include "global.h"
#include "GPU/cudadefs.h"
#include "su3.h"
#include "geometry_eo.h"
#include "start.h"
#include "temporalgauge.h"
#include "stdio.h"
#include "stdlib.h"

su3 * g_trafo;
su3 * tempgauge_field = NULL;



static su3 unit_su3(void)
{
   su3 u;

   u.c00.re=1.0;
   u.c00.im=0.0;
   u.c01.re=0.0;
   u.c01.im=0.0;
   u.c02.re=0.0;
   u.c02.im=0.0;

   u.c10.re=0.0;
   u.c10.im=0.0;
   u.c11.re=1.0;
   u.c11.im=0.0;
   u.c12.re=0.0;
   u.c12.im=0.0;

   u.c20.re=0.0;
   u.c20.im=0.0;
   u.c21.re=0.0;
   u.c21.im=0.0;
   u.c22.re=1.0;
   u.c22.im=0.0;

   return(u);
}



/*copy a complete gauge field*/
/* THINK OF PARALLELIZATION (RAND!!!)*/
void copy_gauge_field(su3** to, su3** from){
  int ix;
  for(ix=0; ix<VOLUME; ix++){
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
int init_temporalgauge_trafo(const int V, su3** gfield){
   int it, iz, iy, ix;

   
   if((void*)(g_trafo = (su3*)calloc(V, sizeof(su3))) == NULL) {
    printf("malloc error in 'init_temporalgauge_trafo'\n"); 
    return(2);
  }
  
  
  /* initialize first timeslice (t=0) with unit matrices*/
  for(ix=0; ix<LX; ix++){
    for(iy=0; iy<LY; iy++){
      for(iz=0; iz<LZ; iz++){
        g_trafo[g_ipt[0][ix][iy][iz]] = unit_su3();
      }
    }
  }
  
  
  /* U^{'}_0(x)  g(x) U_0(x) g^{+}(x+0) != ID   =>  g_(x+0) = g(x) U_0(x)  */
  for(it=1; it<T; it++){
    for(ix=0; ix<LX; ix++){
      for(iy=0; iy<LY; iy++){
        for(iz=0; iz<LZ; iz++){
          _su3_times_su3( g_trafo[g_ipt[it][ix][iy][iz] ] , 
                          g_trafo[ g_ipt[it-1][ix][iy][iz] ] ,
                          gfield[g_ipt[it-1][ix][iy][iz]][0]  
                        );
          
        }
      }
    } 
  }


  /* 
    allocate and initialize g_tempgauge_field which holds a copy of the 
    global gauge field g_gauge_field which is copied back after the inversion
    when the temporal gauge is undone again
  */
  int i=0;
  if((void*)(g_tempgauge_field = (su3**)calloc(V, sizeof(su3*))) == NULL) {
    printf ("malloc error in 'init_temporalgauge_trafo'\n"); 
    return(1);
  }
  if((void*)(tempgauge_field = (su3*)calloc(4*V+1, sizeof(su3))) == NULL) {
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



void finalize_temporalgauge(){
  free(g_trafo);
  free(tempgauge_field);
  free(g_tempgauge_field);
}



/*
  apply gauge transform to gfield with the trafo stored in trafofield
*/
void apply_gtrafo(su3 ** gfield, su3 * trafofield){
 int it, iz, iy, ix, xpos, mu; 
 su3 temp1;
 if(g_proc_id == 0) {
   printf("Applying gauge transformation...");
 }
  for(it=0; it<T; it++){
    for(ix=0; ix<LX; ix++){
      for(iy=0; iy<LY; iy++){
        for(iz=0; iz<LZ; iz++){
          xpos = g_ipt[it][ix][iy][iz];
          for(mu=0;mu<4;mu++){
            /* help = g(x) U_mu(x) */
            _su3_times_su3( temp1, trafofield[xpos],  gfield[xpos][mu]  );
            /* U_mu(x) <- U_mu^{'}(x) = help g^{+}(x+mu)*/
            _su3_times_su3d( gfield[xpos][mu],temp1, trafofield[ g_iup[xpos][mu]  ]);
          }
        }
      }
    } 
  }
   if(g_proc_id == 0) {
   printf("done\n");
 }
 /* update gauge copy fields in the next call to HoppingMatrix */
 g_update_gauge_copy = 1;
}




/*
  apply the inverse gauge transform to gfield with the trafo stored in trafofield
*/
void apply_inv_gtrafo(su3 ** gfield, su3 * trafofield){
 int it, iz, iy, ix, xpos, mu; 
 su3 temp1, temp2;
 if(g_proc_id == 0) {
   printf("Applying INVERSE gauge transformation...");
 }
  for(it=0; it<T; it++){
    for(ix=0; ix<LX; ix++){
      for(iy=0; iy<LY; iy++){
        for(iz=0; iz<LZ; iz++){
          xpos = g_ipt[it][ix][iy][iz];
          for(mu=0;mu<4;mu++){
            /*
            _su3d_times_su3( temp1, trafofield[xpos],  gfield[xpos][mu]  );

            _su3_times_su3( gfield[xpos][mu],temp1, trafofield[ g_iup[xpos][mu]  ]);            
            */
           
           /* help = U^{'}_mu(x) g(x+mu)*/
            _su3_times_su3( temp1,  gfield[xpos][mu], trafofield[ g_iup[xpos][mu]]  );

            /* U_mu(x) <- g^{+}(x) help */
            _su3_dagger(temp2, trafofield[xpos]  )
            _su3_times_su3( gfield[xpos][mu], temp2, temp1);


          }
        }
      }
    } 
  }
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
void apply_inv_gtrafo_spinor(spinor * spin, su3 * trafofield){
 int it, iz, iy, ix, xpos; 
 spinor temp;
 if(g_proc_id == 0) {
   printf("Applying INVERSE gauge transformation to spinor...");
 }
  for(it=0; it<T; it++){
    for(ix=0; ix<LX; ix++){
      for(iy=0; iy<LY; iy++){
        for(iz=0; iz<LZ; iz++){
          xpos = g_ipt[it][ix][iy][iz];
          _su3_inverse_multiply(temp.s0, trafofield[xpos], spin[xpos].s0);
          _su3_inverse_multiply(temp.s1, trafofield[xpos], spin[xpos].s1);
          _su3_inverse_multiply(temp.s2, trafofield[xpos], spin[xpos].s2);
          _su3_inverse_multiply(temp.s3, trafofield[xpos], spin[xpos].s3);
          
          _vector_assign(spin[xpos].s0,temp.s0);
          _vector_assign(spin[xpos].s1,temp.s1);
          _vector_assign(spin[xpos].s2,temp.s2);
          _vector_assign(spin[xpos].s3,temp.s3);
        }
      }
    } 
  }
   if(g_proc_id == 0) {
   printf("done\n");
 }
}






/* 
  apply gauge transform to spinor 
  U^{'}_0(x) = g(x) U_0(x) g^{+}(x+0)
  => psi^{'}(x) = g(x) psi(x)
  (the primed (^{'}) quantities are the gauge transformed fields)
*/
void apply_gtrafo_spinor(spinor * spin, su3 * trafofield){
 int it, iz, iy, ix, xpos; 
 spinor temp;
 
 if(g_proc_id == 0) {
   printf("Applying gauge transformation to spinor...");
 }
  for(it=0; it<T; it++){
    for(ix=0; ix<LX; ix++){
      for(iy=0; iy<LY; iy++){
        for(iz=0; iz<LZ; iz++){
          xpos = g_ipt[it][ix][iy][iz];
          _su3_multiply(temp.s0, trafofield[xpos], spin[xpos].s0);
          _su3_multiply(temp.s1, trafofield[xpos], spin[xpos].s1);
          _su3_multiply(temp.s2, trafofield[xpos], spin[xpos].s2);
          _su3_multiply(temp.s3, trafofield[xpos], spin[xpos].s3);
          
            _vector_assign(spin[xpos].s0,temp.s0);
            _vector_assign(spin[xpos].s1,temp.s1);
            _vector_assign(spin[xpos].s2,temp.s2);
            _vector_assign(spin[xpos].s3,temp.s3);
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
void apply_gtrafo_spinor_odd(spinor * spin, su3 * trafofield){
 int it, iz, iy, ix, xpos, oddpos; 
 spinor temp;

 if(g_proc_id == 0) {
   printf("Applying  gauge transformation to odd spinor...");
 }
  for(it=0; it<T; it++){
    for(ix=0; ix<LX; ix++){ 
      for(iy=0; iy<LY; iy++){
        for(iz=0; iz<LZ; iz++){
          if(((it+ix+iy+iz)%2 != 0)){
            /* odd positions */
            xpos = g_ipt[it][ix][iy][iz];
            oddpos = g_lexic2eosub[ xpos  ];

            _su3_multiply(temp.s0, trafofield[xpos], spin[oddpos].s0);
            _su3_multiply(temp.s1, trafofield[xpos], spin[oddpos].s1);
            _su3_multiply(temp.s2, trafofield[xpos], spin[oddpos].s2);
            _su3_multiply(temp.s3, trafofield[xpos], spin[oddpos].s3);
            
            _vector_assign(spin[oddpos].s0,temp.s0);
            _vector_assign(spin[oddpos].s1,temp.s1);
            _vector_assign(spin[oddpos].s2,temp.s2);
            _vector_assign(spin[oddpos].s3,temp.s3);
            
          }
        }
      }
    } 
  }
   if(g_proc_id == 0) {
   printf("done\n");
 }
}



/* 
  apply inverse gauge transform to ODD spinor 
  U_0(x) = g^{+}(x) U^{'}_0(x) g(x+0)
  => psi(x) = g^{+}(x) psi^{'}(x)
  (the primed (^{'}) quantities are the gauge ttemp.s0ransformed fields)
*/
void apply_inv_gtrafo_spinor_odd(spinor * spin, su3 * trafofield){
 int it, iz, iy, ix, xpos, oddpos; 
 spinor temp;
 if(g_proc_id == 0) {
   printf("Applying INVERSE gauge transformation to odd spinor...");
 }
  for(it=0; it<T; it++){
    for(ix=0; ix<LX; ix++){
      for(iy=0; iy<LY; iy++){
        for(iz=0; iz<LZ; iz++){
          if(((it+ix+iy+iz)%2 != 0)){
            /* odd positions */
            xpos = g_ipt[it][ix][iy][iz];
            oddpos = g_lexic2eosub[ xpos ];
            
            _su3_inverse_multiply(temp.s0, trafofield[xpos], spin[oddpos].s0);
            _su3_inverse_multiply(temp.s1, trafofield[xpos], spin[oddpos].s1);
            _su3_inverse_multiply(temp.s2, trafofield[xpos], spin[oddpos].s2);
            _su3_inverse_multiply(temp.s3, trafofield[xpos], spin[oddpos].s3);
            
            _vector_assign(spin[oddpos].s0,temp.s0);
            _vector_assign(spin[oddpos].s1,temp.s1);
            _vector_assign(spin[oddpos].s2,temp.s2);
            _vector_assign(spin[oddpos].s3,temp.s3);
            
            
          }
        }
      }
    } 
  }
   if(g_proc_id == 0) {
   printf("done\n");
 }
}














/* 
  apply gauge transform to EVENspinor 
  U^{'}_0(x) = g(x) U_0(x) g^{+}(x+0)
  => psi^{'}(x) = g(x) psi(x)
  (the primed (^{'}) quantities are the gauge transformed fields)
*/
void apply_gtrafo_spinor_even(spinor * spin, su3 * trafofield){
 int it, iz, iy, ix, xpos, evenpos; 
 spinor temp;

 if(g_proc_id == 0) {
   printf("Applying  gauge transformation to even spinor...");
 }
  for(it=0; it<T; it++){
    for(ix=0; ix<LX; ix++){ 
      for(iy=0; iy<LY; iy++){
        for(iz=0; iz<LZ; iz++){ 
          if(((it+ix+iy+iz)%2 == 0)){
            /* even positions */
            xpos = g_ipt[it][ix][iy][iz];
            evenpos = g_lexic2eosub[ xpos  ];

            _su3_multiply(temp.s0, trafofield[xpos], spin[evenpos].s0);
            _su3_multiply(temp.s1, trafofield[xpos], spin[evenpos].s1);
            _su3_multiply(temp.s2, trafofield[xpos], spin[evenpos].s2);
            _su3_multiply(temp.s3, trafofield[xpos], spin[evenpos].s3);
            
            _vector_assign(spin[evenpos].s0,temp.s0);
            _vector_assign(spin[evenpos].s1,temp.s1);
            _vector_assign(spin[evenpos].s2,temp.s2);
            _vector_assign(spin[evenpos].s3,temp.s3);
            
          }
        }
      }
    } 
  }
   if(g_proc_id == 0) {
   printf("done\n");
 }
}



/* 
  apply inverse gauge transform to EVEN spinor 
  U_0(x) = g^{+}(x) U^{'}_0(x) g(x+0)
  => psi(x) = g^{+}(x) psi^{'}(x)
  (the primed (^{'}) quantities are the gauge transformed fields)
*/
void apply_inv_gtrafo_spinor_even(spinor * spin, su3 * trafofield){
 int it, iz, iy, ix, xpos, evenpos; 
 spinor temp;
 if(g_proc_id == 0) {
   printf("Applying INVERSE gauge transformation to even spinor...");
 }
  for(it=0; it<T; it++){
    for(ix=0; ix<LX; ix++){
      for(iy=0; iy<LY; iy++){
        for(iz=0; iz<LZ; iz++){
          if(((it+ix+iy+iz)%2 == 0)){
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
   if(g_proc_id == 0) {
   printf("done\n");
 }
}








