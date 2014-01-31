



#define _vector_tensor_vector_add_dev(t, u, v, w, z) \
   (t)[0][0].re=(u).c0.re*(v).c0.re+(u).c0.im*(v).c0.im + (w).c0.re*(z).c0.re+(w).c0.im*(z).c0.im; \
   (t)[0][0].im=(u).c0.im*(v).c0.re-(u).c0.re*(v).c0.im + (w).c0.im*(z).c0.re-(w).c0.re*(z).c0.im; \
   (t)[0][1].re=(u).c0.re*(v).c1.re+(u).c0.im*(v).c1.im + (w).c0.re*(z).c1.re+(w).c0.im*(z).c1.im; \
   (t)[0][1].im=(u).c0.im*(v).c1.re-(u).c0.re*(v).c1.im + (w).c0.im*(z).c1.re-(w).c0.re*(z).c1.im; \
   (t)[0][2].re=(u).c0.re*(v).c2.re+(u).c0.im*(v).c2.im + (w).c0.re*(z).c2.re+(w).c0.im*(z).c2.im; \
   (t)[0][2].im=(u).c0.im*(v).c2.re-(u).c0.re*(v).c2.im + (w).c0.im*(z).c2.re-(w).c0.re*(z).c2.im; \
   (t)[1][0].re=(u).c1.re*(v).c0.re+(u).c1.im*(v).c0.im + (w).c1.re*(z).c0.re+(w).c1.im*(z).c0.im; \
   (t)[1][0].im=(u).c1.im*(v).c0.re-(u).c1.re*(v).c0.im + (w).c1.im*(z).c0.re-(w).c1.re*(z).c0.im; \
   (t)[1][1].re=(u).c1.re*(v).c1.re+(u).c1.im*(v).c1.im + (w).c1.re*(z).c1.re+(w).c1.im*(z).c1.im; \
   (t)[1][1].im=(u).c1.im*(v).c1.re-(u).c1.re*(v).c1.im + (w).c1.im*(z).c1.re-(w).c1.re*(z).c1.im; \
   (t)[1][2].re=(u).c1.re*(v).c2.re+(u).c1.im*(v).c2.im + (w).c1.re*(z).c2.re+(w).c1.im*(z).c2.im; \
   (t)[1][2].im=(u).c1.im*(v).c2.re-(u).c1.re*(v).c2.im + (w).c1.im*(z).c2.re-(w).c1.re*(z).c2.im; \
   (t)[2][0].re=(u).c2.re*(v).c0.re+(u).c2.im*(v).c0.im + (w).c2.re*(z).c0.re+(w).c2.im*(z).c0.im; \
   (t)[2][0].im=(u).c2.im*(v).c0.re-(u).c2.re*(v).c0.im + (w).c2.im*(z).c0.re-(w).c2.re*(z).c0.im; \
   (t)[2][1].re=(u).c2.re*(v).c1.re+(u).c2.im*(v).c1.im + (w).c2.re*(z).c1.re+(w).c2.im*(z).c1.im; \
   (t)[2][1].im=(u).c2.im*(v).c1.re-(u).c2.re*(v).c1.im + (w).c2.im*(z).c1.re-(w).c2.re*(z).c1.im; \
   (t)[2][2].re=(u).c2.re*(v).c2.re+(u).c2.im*(v).c2.im + (w).c2.re*(z).c2.re+(w).c2.im*(z).c2.im; \
   (t)[2][2].im=(u).c2.im*(v).c2.re-(u).c2.re*(v).c2.im + (w).c2.im*(z).c2.re-(w).c2.re*(z).c2.im; 


#define _vector_add_dev(r,s1,s2) \
   (*r).c0.re=(*s1).c0.re+(*s2).c0.re; \
   (*r).c0.im=(*s1).c0.im+(*s2).c0.im; \
   (*r).c1.re=(*s1).c1.re+(*s2).c1.re; \
   (*r).c1.im=(*s1).c1.im+(*s2).c1.im; \
   (*r).c2.re=(*s1).c2.re+(*s2).c2.re; \
   (*r).c2.im=(*s1).c2.im+(*s2).c2.im;

#define _vector_i_add_dev(r,s1,s2) \
   (*r).c0.re=(*s1).c0.re-(*s2).c0.im; \
   (*r).c0.im=(*s1).c0.im+(*s2).c0.re; \
   (*r).c1.re=(*s1).c1.re-(*s2).c1.im; \
   (*r).c1.im=(*s1).c1.im+(*s2).c1.re; \
   (*r).c2.re=(*s1).c2.re-(*s2).c2.im; \
   (*r).c2.im=(*s1).c2.im+(*s2).c2.re;

#define _vector_sub_dev(r,s1,s2) \
   (*r).c0.re=(*s1).c0.re-(*s2).c0.re; \
   (*r).c0.im=(*s1).c0.im-(*s2).c0.im; \
   (*r).c1.re=(*s1).c1.re-(*s2).c1.re; \
   (*r).c1.im=(*s1).c1.im-(*s2).c1.im; \
   (*r).c2.re=(*s1).c2.re-(*s2).c2.re; \
   (*r).c2.im=(*s1).c2.im-(*s2).c2.im;

#define _vector_i_sub_dev(r,s1,s2)	    \
   (*r).c0.re=(*s1).c0.re+(*s2).c0.im; \
   (*r).c0.im=(*s1).c0.im-(*s2).c0.re; \
   (*r).c1.re=(*s1).c1.re+(*s2).c1.im; \
   (*r).c1.im=(*s1).c1.im-(*s2).c1.re; \
   (*r).c2.re=(*s1).c2.re+(*s2).c2.im; \
   (*r).c2.im=(*s1).c2.im-(*s2).c2.re;

#define _complex_times_su3_dev(r,c,s) \
   (r)[0][0].re=(c).re*(s)[0][0].re-(c).im*(s)[0][0].im; \
   (r)[0][0].im=(c).re*(s)[0][0].im+(c).im*(s)[0][0].re; \
   (r)[0][1].re=(c).re*(s)[0][1].re-(c).im*(s)[0][1].im; \
   (r)[0][1].im=(c).re*(s)[0][1].im+(c).im*(s)[0][1].re; \
   (r)[0][2].re=(c).re*(s)[0][2].re-(c).im*(s)[0][2].im; \
   (r)[0][2].im=(c).re*(s)[0][2].im+(c).im*(s)[0][2].re; \
   (r)[1][0].re=(c).re*(s)[1][0].re-(c).im*(s)[1][0].im; \
   (r)[1][0].im=(c).re*(s)[1][0].im+(c).im*(s)[1][0].re; \
   (r)[1][1].re=(c).re*(s)[1][1].re-(c).im*(s)[1][1].im; \
   (r)[1][1].im=(c).re*(s)[1][1].im+(c).im*(s)[1][1].re; \
   (r)[1][2].re=(c).re*(s)[1][2].re-(c).im*(s)[1][2].im; \
   (r)[1][2].im=(c).re*(s)[1][2].im+(c).im*(s)[1][2].re; \
   (r)[2][0].re=(c).re*(s)[2][0].re-(c).im*(s)[2][0].im; \
   (r)[2][0].im=(c).re*(s)[2][0].im+(c).im*(s)[2][0].re; \
   (r)[2][1].re=(c).re*(s)[2][1].re-(c).im*(s)[2][1].im; \
   (r)[2][1].im=(c).re*(s)[2][1].im+(c).im*(s)[2][1].re; \
   (r)[2][2].re=(c).re*(s)[2][2].re-(c).im*(s)[2][2].im; \
   (r)[2][2].im=(c).re*(s)[2][2].im+(c).im*(s)[2][2].re; 


#define _trace_lambda_add_assign_dev(r,a) \
(r).d1+= (-(a)[1][0].im-(a)[0][1].im); \
(r).d2+= (+(a)[1][0].re-(a)[0][1].re); \
(r).d3+= (-(a)[0][0].im+(a)[1][1].im); \
(r).d4+= (-(a)[2][0].im-(a)[0][2].im); \
(r).d5+= (+(a)[2][0].re-(a)[0][2].re); \
(r).d6+= (-(a)[2][1].im-(a)[1][2].im); \
(r).d7+= (+(a)[2][1].re-(a)[1][2].re); \
(r).d8+= ((-(a)[0][0].im-(a)[1][1].im + 2.0*a[2][2].im)*0.577350269189625);


#define _trace_lambda_mul_add_assign_dev(r,c,a) \
(r).d1+= c*(-(a)[1][0].im-(a)[0][1].im); \
(r).d2+= c*(+(a)[1][0].re-(a)[0][1].re); \
(r).d3+= c*(-(a)[0][0].im+(a)[1][1].im); \
(r).d4+= c*(-(a)[2][0].im-(a)[0][2].im); \
(r).d5+= c*(+(a)[2][0].re-(a)[0][2].re); \
(r).d6+= c*(-(a)[2][1].im-(a)[1][2].im); \
(r).d7+= c*(+(a)[2][1].re-(a)[1][2].re); \
(r).d8+= c*((-(a)[0][0].im-(a)[1][1].im + 2.0*a[2][2].im)*0.577350269189625);


//this adds two su3 vectors where a and b are given as spinors 
//and their associated color vector is given as the component s_a and s_b
__device__ void dev_vector_add_spinorcomponent(dev_su3_vec_d * out, dev_spinor_d * a, int s_a, dev_spinor_d * b, int s_b){

 dev_su3_vec_d *col_a, * col_b;
 //convert into double pointer and add 6(no of doubles in dev_su3_vec_d) * (spinor number)
 col_a = (dev_su3_vec_d*) ( ( (double*) a ) + 6*s_a);
 col_b = (dev_su3_vec_d*) ( ( (double*) b ) + 6*s_b);
 
 _vector_add_dev(out, col_a, col_b);
}

//this adds two su3 vectors where a and b are given as spinors 
//and their associated color vector is given as the component s_a and s_b
//imaginary unit included
__device__ void dev_vector_i_add_spinorcomponent(dev_su3_vec_d * out, dev_spinor_d * a, int s_a, dev_spinor_d * b, int s_b){

 dev_su3_vec_d *col_a, * col_b;
 //convert into double pointer and add 6(no of doubles in dev_su3_vec_d) * (spinor number)
 col_a = (dev_su3_vec_d*) ( ( (double*) a ) + 6*s_a);
 col_b = (dev_su3_vec_d*) ( ( (double*) b ) + 6*s_b);
 
 _vector_i_add_dev(out, col_a, col_b);
}

//this subtracts two su3 vectors where a and b are given as spinors 
//and their associated color vector is given as the component s_a and s_b
__device__ void dev_vector_sub_spinorcomponent(dev_su3_vec_d * out, dev_spinor_d * a, int s_a, dev_spinor_d * b, int s_b){

 dev_su3_vec_d *col_a, * col_b;
 //convert into double pointer and add 6(no of doubles in dev_su3_vec_d) * (spinor number)
 col_a = (dev_su3_vec_d*) ( ( (double*) a ) + 6*s_a);
 col_b = (dev_su3_vec_d*) ( ( (double*) b ) + 6*s_b);
 
 _vector_sub_dev(out, col_a, col_b);
}

//this subtracts two su3 vectors where a and b are given as spinors 
//and their associated color vector is given as the component s_a and s_b
//imaginary unit included
__device__ void dev_vector_i_sub_spinorcomponent(dev_su3_vec_d * out, dev_spinor_d * a, int s_a, dev_spinor_d * b, int s_b){

 dev_su3_vec_d *col_a, * col_b;
 //convert into double pointer and add 6(no of doubles in dev_su3_vec_d) * (spinor number)
 col_a = (dev_su3_vec_d*) ( ( (double*) a ) + 6*s_a);
 col_b = (dev_su3_vec_d*) ( ( (double*) b ) + 6*s_b);
 
 _vector_i_sub_dev(out, col_a, col_b);
}



__global__ void dev_deriv_Sb(dev_su3adj* df0, dev_su3_2v_d * gf,dev_spinor_d * l, dev_spinor_d * k, double factor, 
                              const int * gfindex_site, const int* gfindex_nextsite, 
                              const int * nn_evenodd, const int eo, int start, int size) {



  dev_su3_vec_d psia,psib,phia,phib;
  double4 rr[6];
  double4 sp[6];
  double4 sm[6];

  dev_su3_d v1,v2;
  dev_su3_d up, um;


  int pos,hoppos;
  pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;

  if(pos < (start + size)){

    //ix=g_eo2lexic[icx];
    //rr = (*(l + (icx-ioff)));


    /*multiply the left vector with gamma5*/
    //_vector_minus_assign(rr.s2, rr.s2);
    //_vector_minus_assign(rr.s3, rr.s3);
    dev_read_spinor_d(&(rr[0]), &(l[pos]));
    dev_Gamma5_d(&(rr[0]));


    /*********************** direction +0 ********************/

    //iy=g_iup[ix][0]; icy=g_lexic2eosub[iy];
    //sp = k + icy;
    hoppos = nn_evenodd[8*pos];
    dev_read_spinor_d(&(sp[0]), &(k[hoppos]));

    //up=&g_gauge_field[ix][0];
    dev_reconstructgf_2vtexref_d(gf,4*gfindex_site[pos]+0,&(up));

    
    //_vector_add(psia,(*sp).s0,(*sp).s2);
    //_vector_add(psib,(*sp).s1,(*sp).s3);
    dev_vector_add_spinorcomponent(&(psia), &(sp[0]), 0, &(sp[0]), 2);
    dev_vector_add_spinorcomponent(&(psib), &(sp[0]), 1, &(sp[0]), 3);  

    //_vector_add(phia,rr.s0,rr.s2);
    //_vector_add(phib,rr.s1,rr.s3);
    dev_vector_add_spinorcomponent(&(phia), &(rr[0]), 0, &(rr[0]), 2);
    dev_vector_add_spinorcomponent(&(phib), &(rr[0]), 1, &(rr[0]), 3);


    //_vector_tensor_vector_add_dev(v1, phia, psia, phib, psib);
    //_su3_times_su3d(v2,*up,v1);
    //_complex_times_su3(v1,ka0,v2);
    //_trace_lambda_add_assign(df0[ix][0], v1);
    _vector_tensor_vector_add_dev(v1, phia, psia, phib, psib);
    dev_su3_ti_su3d_d(&(v2), &(up), &(v1));
    _complex_times_su3_dev(v1,dev_k0c_d,v2);
    _trace_lambda_mul_add_assign_dev(df0[4*gfindex_site[pos]+0], 2.0*factor, v1);

    /************** direction -0 ****************************/

    //iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];
    //sm = k + icy;
    hoppos = nn_evenodd[8*pos+4];
    dev_read_spinor_d(&(sm[0]), &(k[hoppos]));

    //um=&g_gauge_field[iy][0];
    dev_reconstructgf_2vtexref_d(gf,4*gfindex_nextsite[hoppos]+0,&(um));
 
    //_vector_sub(psia,(*sm).s0,(*sm).s2);
    //_vector_sub(psib,(*sm).s1,(*sm).s3);
    dev_vector_sub_spinorcomponent(&(psia), &(sm[0]), 0, &(sm[0]), 2);
    dev_vector_sub_spinorcomponent(&(psib), &(sm[0]), 1, &(sm[0]), 3);

    //_vector_sub(phia,rr.s0,rr.s2);
    //_vector_sub(phib,rr.s1,rr.s3);
    dev_vector_sub_spinorcomponent(&(phia), &(rr[0]), 0, &(rr[0]), 2);
    dev_vector_sub_spinorcomponent(&(phib), &(rr[0]), 1, &(rr[0]), 3);

    //_vector_tensor_vector_add(v1, psia, phia, psib, phib);
    //_su3_times_su3d(v2,*um,v1);
    //_complex_times_su3(v1,ka0,v2);
    //_trace_lambda_add_assign(df0[iy][0], v1);
    _vector_tensor_vector_add_dev(v1, psia, phia, psib, phib);
    dev_su3_ti_su3d_d(&(v2), &(um), &(v1));
    _complex_times_su3_dev(v1,dev_k0c_d,v2);
    _trace_lambda_mul_add_assign_dev(df0[4*gfindex_nextsite[hoppos]+0], 2.0*factor, v1);


    /*************** direction +1 **************************/

    //iy=g_iup[ix][1]; icy=g_lexic2eosub[iy];
    //sp = k + icy;
    hoppos = nn_evenodd[8*pos+1];
    dev_read_spinor_d(&(sp[0]), &(k[hoppos]));

    //up=&g_gauge_field[ix][1];      
    dev_reconstructgf_2vtexref_d(gf,4*gfindex_site[pos]+1,&(up));

    //_vector_i_add(psia,(*sp).s0,(*sp).s3);
    //_vector_i_add(psib,(*sp).s1,(*sp).s2);
    dev_vector_i_add_spinorcomponent(&(psia), &(sp[0]), 0, &(sp[0]), 3);
    dev_vector_i_add_spinorcomponent(&(psib), &(sp[0]), 1, &(sp[0]), 2); 

    //_vector_i_add(phia,rr.s0,rr.s3);
    //_vector_i_add(phib,rr.s1,rr.s2);
    dev_vector_i_add_spinorcomponent(&(phia), &(rr[0]), 0, &(rr[0]), 3);
    dev_vector_i_add_spinorcomponent(&(phib), &(rr[0]), 1, &(rr[0]), 2); 

    //_vector_tensor_vector_add(v1, phia, psia, phib, psib);
    //_su3_times_su3d(v2,*up,v1);
    //_complex_times_su3(v1,ka1,v2);
    //_trace_lambda_add_assign(df0[ix][1], v1);
    _vector_tensor_vector_add_dev(v1, phia, psia, phib, psib);
    dev_su3_ti_su3d_d(&(v2), &(up), &(v1));
    _complex_times_su3_dev(v1,dev_k1c_d,v2);
    _trace_lambda_mul_add_assign_dev(df0[4*gfindex_site[pos]+1], 2.0*factor, v1);


    
    /**************** direction -1 *************************/

    //iy=g_idn[ix][1]; icy=g_lexic2eosub[iy];
    //sm = k + icy;
    hoppos = nn_evenodd[8*pos+5];
    dev_read_spinor_d(&(sm[0]), &(k[hoppos])); 
  
    //um=&g_gauge_field[iy][1];
    dev_reconstructgf_2vtexref_d(gf,4*gfindex_nextsite[hoppos]+1,&(um));

    //_vector_i_sub(psia,(*sm).s0,(*sm).s3);
    //_vector_i_sub(psib,(*sm).s1,(*sm).s2);
    dev_vector_i_sub_spinorcomponent(&(psia), &(sm[0]), 0, &(sm[0]), 3);
    dev_vector_i_sub_spinorcomponent(&(psib), &(sm[0]), 1, &(sm[0]), 2); 

    //_vector_i_sub(phia,rr.s0,rr.s3);
    //_vector_i_sub(phib,rr.s1,rr.s2);
    dev_vector_i_sub_spinorcomponent(&(phia), &(rr[0]), 0, &(rr[0]), 3);
    dev_vector_i_sub_spinorcomponent(&(phib), &(rr[0]), 1, &(rr[0]), 2); 

    //_vector_tensor_vector_add(v1, psia, phia, psib, phib);
    //_su3_times_su3d(v2,*um,v1);
    //_complex_times_su3(v1,ka1,v2);
    //_trace_lambda_add_assign(df0[iy][1], v1);
    _vector_tensor_vector_add_dev(v1, psia, phia, psib, phib);
    dev_su3_ti_su3d_d(&(v2), &(um), &(v1));
    _complex_times_su3_dev(v1,dev_k1c_d,v2);
    _trace_lambda_mul_add_assign_dev(df0[4*gfindex_nextsite[hoppos]+1], 2.0*factor, v1);


    /*************** direction +2 **************************/

    //iy=g_iup[ix][2]; icy=g_lexic2eosub[iy];
    //sp = k + icy;
    hoppos = nn_evenodd[8*pos+2];
    dev_read_spinor_d(&(sp[0]), &(k[hoppos]));

    //up=&g_gauge_field[ix][2];
    dev_reconstructgf_2vtexref_d(gf,4*gfindex_site[pos]+2,&(up));

    //_vector_add(psia,(*sp).s0,(*sp).s3);
    //_vector_sub(psib,(*sp).s1,(*sp).s2);
    dev_vector_add_spinorcomponent(&(psia), &(sp[0]), 0, &(sp[0]), 3);
    dev_vector_sub_spinorcomponent(&(psib), &(sp[0]), 1, &(sp[0]), 2);  
     
    //_vector_add(phia,rr.s0,rr.s3);
    //_vector_sub(phib,rr.s1,rr.s2);
    dev_vector_add_spinorcomponent(&(phia), &(rr[0]), 0, &(rr[0]), 3);
    dev_vector_sub_spinorcomponent(&(phib), &(rr[0]), 1, &(rr[0]), 2); 

    //_vector_tensor_vector_add(v1, phia, psia, phib, psib);
    //_su3_times_su3d(v2,*up,v1);
    //_complex_times_su3(v1,ka2,v2);
    //_trace_lambda_add_assign(df0[ix][2], v1);
    _vector_tensor_vector_add_dev(v1, phia, psia, phib, psib);
    dev_su3_ti_su3d_d(&(v2), &(up), &(v1));
    _complex_times_su3_dev(v1,dev_k2c_d,v2);
    _trace_lambda_mul_add_assign_dev(df0[4*gfindex_site[pos]+2], 2.0*factor, v1);

    /***************** direction -2 ************************/

    //iy=g_idn[ix][2]; icy=g_lexic2eosub[iy];
    //sm = k + icy;
    hoppos = nn_evenodd[8*pos+6];
    dev_read_spinor_d(&(sm[0]), &(k[hoppos])); 

    //um=&g_gauge_field[iy][2];
    dev_reconstructgf_2vtexref_d(gf,4*gfindex_nextsite[hoppos]+2,&(um));

    //_vector_sub(psia,(*sm).s0,(*sm).s3);
    //_vector_add(psib,(*sm).s1,(*sm).s2);
    dev_vector_sub_spinorcomponent(&(psia), &(sm[0]), 0, &(sm[0]), 3);
    dev_vector_add_spinorcomponent(&(psib), &(sm[0]), 1, &(sm[0]), 2); 

    //_vector_sub(phia,rr.s0,rr.s3);
    //_vector_add(phib,rr.s1,rr.s2);
    dev_vector_sub_spinorcomponent(&(phia), &(rr[0]), 0, &(rr[0]), 3);
    dev_vector_add_spinorcomponent(&(phib), &(rr[0]), 1, &(rr[0]), 2); 

    //_vector_tensor_vector_add(v1, psia, phia, psib, phib);
    //_su3_times_su3d(v2,*um,v1);
    //_complex_times_su3(v1,ka2,v2);
    //_trace_lambda_add_assign(df0[iy][2], v1);
    _vector_tensor_vector_add_dev(v1, psia, phia, psib, phib);
    dev_su3_ti_su3d_d(&(v2), &(um), &(v1));
    _complex_times_su3_dev(v1,dev_k2c_d,v2);
    _trace_lambda_mul_add_assign_dev(df0[4*gfindex_nextsite[hoppos]+2], 2.0*factor, v1);


    /****************** direction +3 ***********************/

    //iy=g_iup[ix][3]; icy=g_lexic2eosub[iy];
    //sp = k + icy;
    hoppos = nn_evenodd[8*pos+3];
    dev_read_spinor_d(&(sp[0]), &(k[hoppos]));

    //up=&g_gauge_field[ix][3];
    dev_reconstructgf_2vtexref_d(gf,4*gfindex_site[pos]+3,&(up)); 

    //_vector_i_add(psia,(*sp).s0,(*sp).s2);
    //_vector_i_sub(psib,(*sp).s1,(*sp).s3);
    dev_vector_i_add_spinorcomponent(&(psia), &(sp[0]), 0, &(sp[0]), 2);
    dev_vector_i_sub_spinorcomponent(&(psib), &(sp[0]), 1, &(sp[0]), 3);

    //_vector_i_add(phia,rr.s0,rr.s2);
    //_vector_i_sub(phib,rr.s1,rr.s3);
    dev_vector_i_add_spinorcomponent(&(phia), &(rr[0]), 0, &(rr[0]), 2);
    dev_vector_i_sub_spinorcomponent(&(phib), &(rr[0]), 1, &(rr[0]), 3);

    //_vector_tensor_vector_add(v1, phia, psia, phib, psib);
    //_su3_times_su3d(v2,*up,v1);
    //_complex_times_su3(v1,ka3,v2);
    //_trace_lambda_add_assign(df0[ix][3], v1);
    _vector_tensor_vector_add_dev(v1, phia, psia, phib, psib);
    dev_su3_ti_su3d_d(&(v2), &(up), &(v1));
    _complex_times_su3_dev(v1,dev_k3c_d,v2);
    _trace_lambda_mul_add_assign_dev(df0[4*gfindex_site[pos]+3], 2.0*factor, v1);


    /***************** direction -3 ************************/

    //iy=g_idn[ix][3]; icy=g_lexic2eosub[iy];
    //sm = k + icy;
    hoppos = nn_evenodd[8*pos+7];
    dev_read_spinor_d(&(sm[0]), &(k[hoppos])); 

    //um=&g_gauge_field[iy][3];
    dev_reconstructgf_2vtexref_d(gf,4*gfindex_nextsite[hoppos]+3,&(um));

    //_vector_i_sub(psia,(*sm).s0,(*sm).s2);
    //_vector_i_add(psib,(*sm).s1,(*sm).s3);
    dev_vector_i_sub_spinorcomponent(&(psia), &(sm[0]), 0, &(sm[0]), 2);
    dev_vector_i_add_spinorcomponent(&(psib), &(sm[0]), 1, &(sm[0]), 3); 

    //_vector_i_sub(phia,rr.s0,rr.s2);
    //_vector_i_add(phib,rr.s1,rr.s3);
    dev_vector_i_sub_spinorcomponent(&(phia), &(rr[0]), 0, &(rr[0]), 2);
    dev_vector_i_add_spinorcomponent(&(phib), &(rr[0]), 1, &(rr[0]), 3); 

    //_vector_tensor_vector_add(v1, psia, phia, psib, phib);
    //_su3_times_su3d(v2,*um,v1);
    //_complex_times_su3(v1,ka3,v2);
    //_trace_lambda_add_assign(df0[iy][3], v1);
    _vector_tensor_vector_add_dev(v1, psia, phia, psib, phib);
    dev_su3_ti_su3d_d(&(v2), &(um), &(v1));
    _complex_times_su3_dev(v1,dev_k3c_d,v2);
    _trace_lambda_mul_add_assign_dev(df0[4*gfindex_nextsite[hoppos]+3], 2.0*factor, v1); 
  
  }
}











extern "C" void gpu_deriv_Sb(const int ieo, spinor * const l, spinor * const k, 
                             hamiltonian_field_t * const hf, const double factor){
cudaError_t cudaerr;
printf("GPU deriv_Sb...\n");
int host_check_VOL2, Vol;


  #ifdef MPI
   Vol = (VOLUME + RAND)/2;
  #else
   Vol = VOLUME/2;
  #endif


  if((cudaerr=cudaMemcpyToSymbol(dev_VOL2, &(VOLUME), sizeof(int)))!=cudaSuccess){
    printf("gpu_deriv_Sb(): Could not copy dev_VOL2 to device. Aborting...\n");
    exit(200);
  } 
  cudaMemcpyFromSymbol(&host_check_VOL2, dev_VOL2, sizeof(int)); 
  printf("\tOn device:\n");
  printf("\tdev_VOL2 = %i\n", host_check_VOL2);  

  // FLIP IMAGINARY OF ka0
  // this is needed as update_constants_d gives an 
  // extra sign to the imaginary part of all kappas
  // this is compensated by a sign in GPU Hopping matrices
  ka0 = conj(ka0);

  //update constants, gauge field and momentum field
  update_constants_d(dev_grid);
  update_gpu_fields(hf->gaugefield, hf->derivative,1);


//   dev_complex_d h0;
//   h0.re = ka0.re;    h0.im = ka0.im;
//   // try using constant mem for kappas
//   printf("ka0: %.8f  %.8f\n", ka0.re, ka0.im);
//   printf("ka1: %.8f  %.8f\n", ka1.re, ka1.im);
//   if((cudaerr=cudaMemcpyToSymbol(dev_k0_d, &(h0), sizeof(dev_complex_d)))!=cudaSuccess){
//     printf("gpu_deriv_Sb(): Could not copy dev_k0c_d to device. Aborting...\n");
//     exit(200);
//   } 
// 
//   dev_complex host_check_k0, host_check_k1;
//   cudaMemcpyFromSymbol(&host_check_k0, dev_k0_d, sizeof(dev_complex_d));
//   if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
//         printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
//   }
//   printf ("GPU k0.re = %.8f k0.im = %.8f\n",host_check_k0.re, host_check_k0.im); 
// 
// 
//   dev_complex_d h1;
//   h1.re = (double)ka1.re;    h1.im = -(double)ka1.im;
//   // try using constant mem for kappas
//   cudaMemcpyToSymbol("dev_k1_d", &h1, sizeof(dev_complex_d)) ;
//   cudaMemcpyFromSymbol(&host_check_k1, dev_k1_d, sizeof(dev_complex_d));
//   printf ("GPU k1.re = %.8f k1.im = %.8f\n",host_check_k1.re, host_check_k1.im); 


  if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
        printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
  }
//printf("griddim (gauge_monomial): %d, blockdim: %d\n", griddimgauge, blockdimgauge);

  //bring spinors to device RE-ORDER!
  size_t dev_spinsize_d = 6*Vol * sizeof(dev_spinor_d); // double4 even-odd !
  
  // l == dev_spin_eo1_d
  order_spin_gpu(l, h2d_spin_d);
  cudaMemcpy(dev_spin_eo1_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice); 

  // k == dev_spin_eo2_d
  order_spin_gpu(k, h2d_spin_d);
  cudaMemcpy(dev_spin_eo2_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice); 


if(ieo==0){
//cudaFuncSetCacheConfig(dev_gauge_derivative, cudaFuncCachePreferL1);
dev_deriv_Sb<<<griddimgauge, blockdimgauge >>>(dev_df0_d, dev_gf_d, 
                                               dev_spin_eo1_d, dev_spin_eo2_d, factor,
                                               dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, Vol);   
}
if(ieo==1){
//cudaFuncSetCacheConfig(dev_gauge_derivative, cudaFuncCachePreferL1);
dev_deriv_Sb<<<griddimgauge, blockdimgauge >>>(dev_df0_d, dev_gf_d, 
                                               dev_spin_eo1_d, dev_spin_eo2_d, factor, 
                                               dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, Vol);   
}

if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
}           


  cudaMemcpyFromSymbol(&host_check_VOL2, dev_VOL2, sizeof(int)); 
  printf("\tOn device:\n");
  printf("\tdev_VOL2 = %i\n", host_check_VOL2); 

  //bring momentum field back to host
  to_host_mom(hf);

  if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
        printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
        printf("Error code is: %f\n",cudaerr);
  }


  // FLIP BACK IMAGINARY OF ka0
  ka0= conj(ka0);

printf("finished: GPU deriv_Sb.\n");

}



//
extern "C" void gpu_H_deriv_Sb(const int ieo, spinor * const l, spinor * const k,
                               hamiltonian_field_t * const hf, const double factor){
cudaError_t cudaerr;
printf("GPU H_deriv_Sb...\n");
int host_check_VOL2, Vol;


  #ifdef MPI
   Vol = (VOLUME + RAND)/2;
  #else
   Vol = VOLUME/2;
  #endif


  if((cudaerr=cudaMemcpyToSymbol(dev_VOL2, &(VOLUME), sizeof(int)))!=cudaSuccess){
    printf("gpu_deriv_Sb(): Could not copy dev_VOL2 to device. Aborting...\n");
    exit(200);
  } 
  cudaMemcpyFromSymbol(&host_check_VOL2, dev_VOL2, sizeof(int)); 
  printf("\tOn device:\n");
  printf("\tdev_VOL2 = %i\n", host_check_VOL2);  


   size_t dev_spinsize_d = 6*VOLUME/2 * sizeof(dev_spinor_d); // double4 even-odd !   
   int gridsize;
     //this is the partitioning for the HoppingMatrix kernel
     int blockdim3 = BLOCKD;
     if( VOLUME/2 % blockdim3 == 0){
       gridsize = (int) VOLUME/2/blockdim3;
     }
     else{
       gridsize = (int) VOLUME/2/blockdim3 + 1;
     }
     int griddim3 = gridsize;
   
     //this is the partitioning for dev_mul_one_pm...
     int blockdim4 = BLOCK2D;
     if( VOLUME/2 % blockdim4 == 0){
       gridsize = (int) VOLUME/2/blockdim4;
     }
     else{
       gridsize = (int) VOLUME/2/blockdim4 + 1;
     }
     int griddim4 = gridsize; 



  if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
        printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
  }
//printf("griddim (gauge_monomial): %d, blockdim: %d\n", griddimgauge, blockdimgauge);

  //bring spinors to device RE-ORDER!
  dev_spinsize_d = 6*Vol * sizeof(dev_spinor_d); // double4 even-odd !
  
  order_spin_gpu(l, h2d_spin_d);
  cudaMemcpy(dev_spin0_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice); 


  order_spin_gpu(k, h2d_spin_d);
  cudaMemcpy(dev_spin1_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice); 


if(ieo==0){
//apply H to k -> eo1
dev_H_eo_tm_inv_psi_d(dev_spin1_d, dev_spin_eo1_d, dev_spin_eo2_d, 
		      griddim3, blockdim3, griddim4, blockdim4,
		      dev_eoidx_even, dev_eoidx_odd, 
		      dev_nn_eo, dev_nn_oe, 0, 1.0);
  // FLIP IMAGINARY OF ka0
  // this is needed as update_constants_d gives an 
  // extra sign to the imaginary part of all kappas
  // this is compensated by a sign in GPU Hopping matrices
  ka0 = conj(ka0);

  //update constants, gauge field and momentum field
  update_constants_d(dev_grid);
  update_gpu_fields(hf->gaugefield, hf->derivative,1);
dev_deriv_Sb<<<griddimgauge, blockdimgauge >>>(dev_df0_d, dev_gf_d, 
                                               dev_spin_eo1_d, dev_spin0_d, factor,
                                               dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, Vol);   
}
if(ieo==1){
//apply H to k -> eo1
dev_H_eo_tm_inv_psi_d(dev_spin1_d, dev_spin_eo1_d, dev_spin_eo2_d, 
		      griddim3, blockdim3, griddim4, blockdim4,
		      dev_eoidx_even, dev_eoidx_odd, 
		      dev_nn_eo, dev_nn_oe, 1, -1.0);
  // FLIP IMAGINARY OF ka0
  // this is needed as update_constants_d gives an 
  // extra sign to the imaginary part of all kappas
  // this is compensated by a sign in GPU Hopping matrices
  ka0 = conj(ka0);
  //update constants, gauge field and momentum field
  update_constants_d(dev_grid);
  update_gpu_fields(hf->gaugefield, hf->derivative,1);
  dev_deriv_Sb<<<griddimgauge, blockdimgauge >>>(dev_df0_d, dev_gf_d, 
                                               dev_spin_eo1_d, dev_spin0_d, factor,  
                                               dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, Vol);   
}

if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
}           


  cudaMemcpyFromSymbol(&host_check_VOL2, dev_VOL2, sizeof(int)); 
  printf("\tOn device:\n");
  printf("\tdev_VOL2 = %i\n", host_check_VOL2); 

  //bring momentum field back to host
  to_host_mom(hf);

  if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
        printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
        printf("Error code is: %f\n",cudaerr);
  }


  // FLIP BACK IMAGINARY OF ka0
  ka0 = conj(ka0);

printf("finished: GPU H_deriv_Sb.\n");

}










