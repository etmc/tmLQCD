/***********************************************************************
 *
 * Copyright (C) 2010 Florian Burger
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
 *
 *  
 * File: linalg.cuh
 *
 * CUDA linear algebra functions and implementation of gamma-multiplication
 *
 * 
 *
 **************************************************************************/
 
 
 
template<class RealT> 
__device__ inline dev_complexT<RealT> dev_cconj (dev_complexT<RealT> c){ /*konjugiert komplexe Zahl*/
 dev_complexT<RealT> erg;
 erg.re = c.re;
 erg.im = -1.0*c.im;
return erg;
}

template<class RealTVon,class RealTNach>
__device__ inline void dev_ccopy(dev_complexT<RealTVon>* von, dev_complexT<RealTNach>* nach){/*kopiert complex von nach complex nach*/
  nach->re = RealTNach(von->re);
  nach->im = RealTNach(von->im);
}

template<class RealT>
__device__ inline RealT dev_cabssquare (dev_complexT<RealT> c){ /*gibt abs^2 einer komplexen Zahl zurück*/
 return c.re*c.re + c.im*c.im;
}

template<class RealT>
__device__ inline RealT dev_cabsolute (dev_complexT<RealT> c){/*gibt Betrag einer kompl. zahl zurück*/
 return sqrt(c.re*c.re + c.im*c.im);
}


template<class RealT>
__device__ inline  dev_complexT<RealT> dev_crealmult(dev_complexT<RealT> c1, RealT real){ /*multipliziert c1 mit reeller zahl re*/
  dev_complexT<RealT> erg;
  erg.re = real*c1.re;
  erg.im = real*c1.im;
return erg;
}

template<class RealT>
__device__ inline dev_complexT<RealT> dev_cmult (dev_complexT<RealT> c1, dev_complexT<RealT> c2){ /*multiplizier zwei komplexe Zahlen*/
  dev_complexT<RealT> erg;
  erg.re = c1.re * c2.re - c1.im * c2.im;
  erg.im = c1.re * c2.im + c1.im * c2.re;
return erg;
}

template<class RealT>
__device__ inline dev_complexT<RealT> dev_cadd (dev_complexT<RealT> c1, dev_complexT<RealT> c2){ /*addiert zwei komplexe Zahlen */
  dev_complexT<RealT> erg;
  erg.re = c1.re + c2.re;
  erg.im = c1.im + c2.im;
return erg;
}


template<class RealT>
__device__ inline dev_complexT<RealT> dev_cdiv(dev_complexT<RealT> c1, dev_complexT<RealT> c2) { /* dividiert c1 durch c2 */
  dev_complexT<RealT> erg;
  RealT oneovernenner = 1.0/(c2.re*c2.re + c2.im*c2.im);
  erg.re = oneovernenner*(c1.re*c2.re + c1.im*c2.im);
  erg.im = oneovernenner*(c1.im*c2.re - c1.re*c2.im);
return erg;
}


template<class RealT>
__device__ inline dev_complexT<RealT> dev_csub(dev_complexT<RealT> c1, dev_complexT<RealT> c2){
   dev_complexT<RealT> erg;
   erg.re = c1.re - c2.re;
   erg.im = c1.im - c2.im;
return erg;
}


template<class RealT>
__device__ inline dev_complexT<RealT> dev_initcomplex(RealT re, RealT im){/* gibt komplexe Zahl mit Realt re und Imt im zurück*/
    dev_complexT<RealT> erg;
    erg.re = re;
    erg.im = im;
return (erg);
}





template<class RealT1,class RealT2>
__device__ inline void dev_copy_spinor(typename dev_spinorT<RealT1>::type *i1, typename dev_spinorT<RealT2>::type *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+i)).x = RealT2((*(i1+i)).x);
    (*(i2+i)).y = RealT2((*(i1+i)).y);
    (*(i2+i)).z = RealT2((*(i1+i)).z);
    (*(i2+i)).w = RealT2((*(i1+i)).w);
  }
}

template<class RealT>
__device__ inline void dev_zero_spinor(typename dev_spinorT<RealT>::type *sin){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(sin+i)).x = 0.0;
    (*(sin+i)).y = 0.0;
    (*(sin+i)).z = 0.0;
    (*(sin+i)).w = 0.0;
  }
}






//out = in + lambda in2
template<class RealT>
__device__ inline void dev_skalarmult_add_assign_spinor
(
  typename dev_spinorT<RealT>::type *in, 
  RealT lambda,
  typename dev_spinorT<RealT>::type * in2,
  typename dev_spinorT<RealT>::type * out
){
  int i; 
  #pragma unroll 6
for(i=0;i<6;i++){ //color + spin
    (*(out+i)).x = (*(in+i)).x + lambda* (*(in2+i)).x;
    (*(out+i)).y = (*(in+i)).y + lambda* (*(in2+i)).y;
    (*(out+i)).z = (*(in+i)).z + lambda* (*(in2+i)).z;
    (*(out+i)).w = (*(in+i)).w + lambda* (*(in2+i)).w;
  }
}




//out = in + lambda in2
template<class RealT>
__device__ inline void dev_complexmult_add_assign_spinor
(
  typename dev_spinorT<RealT>::type* in,
  dev_complexT<RealT> lambda,
  typename dev_spinorT<RealT>::type* in2,
  typename dev_spinorT<RealT>::type* out
){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(out+i)).x = (*(in+i)).x + ((*(in2+i)).x*lambda.re - (*(in2+i)).y*lambda.im);
    (*(out+i)).y = (*(in+i)).y + ((*(in2+i)).x*lambda.im + (*(in2+i)).y*lambda.re);
    (*(out+i)).z = (*(in+i)).z + ((*(in2+i)).z*lambda.re - (*(in2+i)).w*lambda.im);
    (*(out+i)).w = (*(in+i)).w + ((*(in2+i)).z*lambda.im + (*(in2+i)).w*lambda.re);
  }
}




//out = in + (lambda)* in2
template<class RealT>
__device__ inline void dev_complexcgmult_add_assign_spinor
(
  typename dev_spinorT<RealT>::type * in,
  dev_complexT<RealT> lambda,
  typename dev_spinorT<RealT>::type* in2,
  typename dev_spinorT<RealT>::type* out
){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(out+i)).x = (*(in+i)).x + ((*(in2+i)).x*lambda.re + (*(in2+i)).y*lambda.im);
    (*(out+i)).y = (*(in+i)).y + (-(*(in2+i)).x*lambda.im + (*(in2+i)).y*lambda.re);
    (*(out+i)).z = (*(in+i)).z + ((*(in2+i)).z*lambda.re + (*(in2+i)).w*lambda.im);
    (*(out+i)).w = (*(in+i)).w + (-(*(in2+i)).z*lambda.im + (*(in2+i)).w*lambda.re);
  }
}



template<class RealT>
__device__ void inline dev_skalarmult_spinor
(
  typename dev_spinorT<RealT>::type* in,
  dev_complexT<RealT> lambda,
  typename dev_spinorT<RealT>::type* out
){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    //out[i] = dev_cmult(in[i],lambda);
    
    (*(out+i)).x = (*(in+i)).x*lambda.re - (*(in+i)).y*lambda.im;
    (*(out+i)).y = (*(in+i)).y*lambda.re + (*(in+i)).x*lambda.im;
    
    (*(out+i)).z = (*(in+i)).z*lambda.re - (*(in+i)).w*lambda.im;
    (*(out+i)).w = (*(in+i)).w*lambda.re + (*(in+i)).z*lambda.im;
  }
}



/*
__device__ void inline dev_skalarmult_gamma5_spinor(dev_spinor * out, const dev_complex lambda, dev_spinor * in){


 (*(out)).x = (*(in)).x*lambda.re;
 (*(out)).x -= (*(in)).y*lambda.im;
 
 (*(out)).y = (*(in)).y*lambda.re;
 (*(out)).y += (*(in)).x*lambda.im;
 
 (*(out)).z = (*(in)).z*lambda.re;
 (*(out)).z -= (*(in)).w*lambda.im;

 (*(out)).w = (*(in)).w*lambda.re;
 (*(out)).w += (*(in)).z*lambda.im;

 
 (*(out+1)).x = (*(in+1)).x*lambda.re;
 (*(out+1)).x -= (*(in+1)).y*lambda.im;
 
 (*(out+1)).y = (*(in+1)).y*lambda.re;
 (*(out+1)).y += (*(in+1)).x*lambda.im;
 
 (*(out+1)).z = (*(in+1)).z*lambda.re;
 (*(out+1)).z -= (*(in+1)).w*lambda.im;

 (*(out+1)).w = (*(in+1)).w*lambda.re;
 (*(out+1)).w += (*(in+1)).z*lambda.im;


 (*(out+2)).x = (*(in+2)).x*lambda.re;
 (*(out+2)).x -= (*(in+2)).y*lambda.im;

 (*(out+2)).y = (*(in+2)).y*lambda.re;
 (*(out+2)).y += (*(in+2)).x*lambda.im;
 
 (*(out+2)).z = (*(in+2)).z*lambda.re;
 (*(out+2)).z -= (*(in+2)).w*lambda.im;

 (*(out+2)).w = (*(in+2)).w*lambda.re;
 (*(out+2)).w += (*(in+2)).z*lambda.im;


 (*(out+3)).x = (*(in+3)).y*lambda.im;
 (*(out+3)).x -= (*(in+3)).x*lambda.re;

 (*(out+3)).y = - (*(in+3)).x*lambda.im;
 (*(out+3)).y -= (*(in+3)).y*lambda.re;
 
 (*(out+3)).z = (*(in+3)).w*lambda.im;
 (*(out+3)).z -= (*(in+3)).z*lambda.re;

 (*(out+3)).w = -(*(in+3)).z*lambda.im;
 (*(out+3)).w -= (*(in+3)).w*lambda.re;


 (*(out+4)).x = (*(in+4)).y*lambda.im;
 (*(out+4)).x -= (*(in+4)).x*lambda.re;

 (*(out+4)).y = - (*(in+4)).x*lambda.im;
 (*(out+4)).y -= (*(in+4)).y*lambda.re;
 
 (*(out+4)).z = (*(in+4)).w*lambda.im;
 (*(out+4)).z -= (*(in+4)).z*lambda.re;

 (*(out+4)).w = -(*(in+4)).z*lambda.im;
 (*(out+4)).w -= (*(in+4)).w*lambda.re;


 (*(out+5)).x = (*(in+5)).y*lambda.im;
 (*(out+5)).x -= (*(in+5)).x*lambda.re;

 (*(out+5)).y = - (*(in+5)).x*lambda.im;
 (*(out+5)).y -= (*(in+5)).y*lambda.re;
 
 (*(out+5)).z = (*(in+5)).w*lambda.im;
 (*(out+5)).z -= (*(in+5)).z*lambda.re;

 (*(out+5)).w = -(*(in+5)).z*lambda.im;
 (*(out+5)).w -= (*(in+5)).w*lambda.re;

}
*/


template<class RealT>
__device__ void inline dev_skalarmult_gamma5_spinor(typename dev_spinorT<RealT>::type* out, dev_complexT<RealT> lambda, typename dev_spinorT<RealT>::type* in){
int i;
 typename dev_spinorT<RealT>::type shelp, tempout;

shelp = *(in);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;
 
 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out)) = tempout;
 
 
shelp = *(in+1);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;
 
 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out+1)) = tempout; 
 
 
shelp = *(in+2);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;

 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out+2)) = tempout; 

 
shelp = *(in+3);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+3)) = tempout;

shelp = *(in+4);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+4)) = tempout;

shelp = *(in+5);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+5)) = tempout;
}



template<class RealT>
__device__ void inline dev_realmult_spinor(typename dev_spinorT<RealT>::type* in, RealT lambda){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    //in[i] = in[i]*lambda;
    (*(in+i)).x = (*(in+i)).x*lambda;
    (*(in+i)).y = (*(in+i)).y*lambda;
    
    (*(in+i)).z = (*(in+i)).z*lambda;
    (*(in+i)).w = (*(in+i)).w*lambda;
  }
}


template<class RealT>
__device__ void inline dev_realmult_spinor_assign(typename dev_spinorT<RealT>::type* out, RealT lambda, typename dev_spinorT<RealT>::type* in){
int i;
#pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
  //out[i] = in[i]*lambda;
      (*(out+i)).x = (*(in+i)).x*lambda;
	  (*(out+i)).y = (*(in+i)).y*lambda;
      
	  (*(out+i)).z = (*(in+i)).z*lambda;
      (*(out+i)).w = (*(in+i)).w*lambda;
  }
}




template<class RealT>
__device__ void dev_assign_realmult_add_spinor
(
  typename dev_spinorT<RealT>::type* out,
  RealT lambda,
  typename dev_spinorT<RealT>::type* in1,
  typename dev_spinorT<RealT>::type* in2
){
int i;
RealT help;
//out = lambda*(in1 + in2)
#pragma unroll 6
  for(i=0;i<6;i++){ //color + spin

      help = (*(in1+i)).x*lambda;
      help += (*(in2+i)).x*lambda;
      (*(out+i)).x = help;
      
      help = (*(in1+i)).y*lambda;
      help += (*(in2+i)).y*lambda;
      (*(out+i)).y = help;      

      help = (*(in1+i)).z*lambda;
      help += (*(in2+i)).z*lambda;
      (*(out+i)).z = help;

      help = (*(in1+i)).w*lambda;
      help += (*(in2+i)).w*lambda;
      (*(out+i)).w = help;
  } 
}


template<class RealT>
__device__ inline void dev_add_spinor_assign(typename dev_spinorT<RealT>::type * i1, typename dev_spinorT<RealT>::type * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x + (*(i2+i)).x;
    (*(i1+i)).y = (*(i1+i)).y + (*(i2+i)).y;
    (*(i1+i)).z = (*(i1+i)).z + (*(i2+i)).z;
    (*(i1+i)).w = (*(i1+i)).w + (*(i2+i)).w;
  }
}



template<class RealT>
__device__ inline void dev_sub_spinor_assign(typename dev_spinorT<RealT>::type * i1, typename dev_spinorT<RealT>::type * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x - (*(i2+i)).x;
    (*(i1+i)).y = (*(i1+i)).y - (*(i2+i)).y;
    (*(i1+i)).z = (*(i1+i)).z - (*(i2+i)).z;
    (*(i1+i)).w = (*(i1+i)).w - (*(i2+i)).w;
  }
}




/*
//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_spintex(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;

s1 = tex1Dfetch(spin_tex,6*pos);
s2 = tex1Dfetch(spin_tex,6*pos+1);


//(*(out+0)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
//(*(out+0)).y = ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );

(*(out+0)).x =  M[0][0].re*s1.x;
  (*(out+0)).y =  M[0][0].re*s1.y;
(*(out+0)).x -= M[0][0].im*s1.y;
  (*(out+0)).y += M[0][0].im*s1.x;
(*(out+0)).x += M[0][1].re*s1.z;
  (*(out+0)).y += M[0][1].re*s1.w;
(*(out+0)).x -= M[0][1].im*s1.w;
  (*(out+0)).y += M[0][1].im*s1.z;
(*(out+0)).x += M[0][2].re*s2.x;
  (*(out+0)).y += M[0][2].re*s2.y;
(*(out+0)).x -= M[0][2].im*s2.y;
  (*(out+0)).y += M[0][2].im*s2.x;



//(*(out+0)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
//(*(out+0)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+0)).z = M[1][0].re*s1.x;
  (*(out+0)).w = M[1][0].re*s1.y;
(*(out+0)).z -= M[1][0].im*s1.y;
  (*(out+0)).w += M[1][0].im*s1.x;
(*(out+0)).z += M[1][1].re*s1.z;
  (*(out+0)).w += M[1][1].re*s1.w;
(*(out+0)).z -= M[1][1].im*s1.w;
  (*(out+0)).w += M[1][1].im*s1.z;
(*(out+0)).z += M[1][2].re*s2.x;
  (*(out+0)).w += M[1][2].re*s2.y;
(*(out+0)).z -= M[1][2].im*s2.y;
  (*(out+0)).w += M[1][2].im*s2.x;



//(*(out+1)).x = ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
//(*(out+1)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );


(*(out+1)).x = M[2][0].re*s1.x;
  (*(out+1)).y = M[2][0].re*s1.y;
(*(out+1)).x -= M[2][0].im*s1.y;
  (*(out+1)).y += M[2][0].im*s1.x;
(*(out+1)).x += M[2][1].re*s1.z;
  (*(out+1)).y += M[2][1].re*s1.w;
(*(out+1)).x -= M[2][1].im*s1.w;
  (*(out+1)).y += M[2][1].im*s1.z;
(*(out+1)).x += M[2][2].re*s2.x;
  (*(out+1)).y += M[2][2].re*s2.y;
(*(out+1)).x -= M[2][2].im*s2.y;
  (*(out+1)).y += M[2][2].im*s2.x;





s1 = tex1Dfetch(spin_tex,6*pos+2);
(*(out+1)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+1)).w =  ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+2)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+2)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+2)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+2)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );



s1 = tex1Dfetch(spin_tex,6*pos+3);
s2 = tex1Dfetch(spin_tex,6*pos+4);
(*(out+3)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+3)).y =   ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );


(*(out+3)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
(*(out+3)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+4)).x =  ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+4)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );



s1 = tex1Dfetch(spin_tex,6*pos+5);
(*(out+4)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+4)).w =   ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+5)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+5)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+5)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+5)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );


}
*/





//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
template<class RealT>
__device__ void dev_su3MtV_spintex(dev_su3M(RealT) M, int pos, dev_spinorM(RealT) * out){

dev_spinorM(RealT) s1, s2;

#ifndef HALF
 s1 = tex1Dfetch(spin_tex,6*pos);
#else
 s1 = tex1Dfetch(spinhalf_tex,6*pos);
 float norm = tex1Dfetch(spinnormhalf_tex,pos);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif

#ifndef HALF
 s2 = tex1Dfetch(spin_tex,6*pos+1);
#else
 s2 = tex1Dfetch(spinhalf_tex,6*pos+1);
 s2.x *= norm; 
 s2.y *= norm; 
 s2.z *= norm; 
 s2.w *= norm;
#endif

(*(out+0)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+0)).y = ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );



(*(out+0)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
(*(out+0)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+1)).x = ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+1)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );


#ifndef HALF
 s1 = tex1Dfetch(spin_tex,6*pos+2);
#else
 s1 = tex1Dfetch(spinhalf_tex,6*pos+2);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif

(*(out+1)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+1)).w =  ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+2)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+2)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+2)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+2)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );


#ifndef HALF
 s1 = tex1Dfetch(spin_tex,6*pos+3);
#else
 s1 = tex1Dfetch(spinhalf_tex,6*pos+3);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif

#ifndef HALF
 s2 = tex1Dfetch(spin_tex,6*pos+4);
#else
 s2 = tex1Dfetch(spinhalf_tex,6*pos+4);
 s2.x *= norm; 
 s2.y *= norm; 
 s2.z *= norm; 
 s2.w *= norm;
#endif
(*(out+3)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+3)).y =   ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );


(*(out+3)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
(*(out+3)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+4)).x =  ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+4)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );


#ifndef HALF
 s1 = tex1Dfetch(spin_tex,6*pos+5);
#else
 s1 = tex1Dfetch(spinhalf_tex,6*pos+5);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif
(*(out+4)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+4)).w =   ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+5)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+5)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+5)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+5)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );


}










//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
template<class RealT>
__device__ void dev_su3MtV(typename dev_su3T<RealT>::type M, const typename dev_spinorT<RealT>::type * s, typename dev_spinorT<RealT>::type * out){

(*(out+0)).x =  ( M[0][0].re*(*(s+0)).x - M[0][0].im*(*(s+0)).y ) + ( M[0][1].re*(*(s+0)).z - M[0][1].im*(*(s+0)).w ) + ( M[0][2].re*(*(s+1)).x - M[0][2].im*(*(s+1)).y );
(*(out+0)).y = ( M[0][0].re*(*(s+0)).y + M[0][0].im*(*(s+0)).x ) + ( M[0][1].re*(*(s+0)).w + M[0][1].im*(*(s+0)).z ) + ( M[0][2].re*(*(s+1)).y + M[0][2].im*(*(s+1)).x );


(*(out+0)).z =  ( M[1][0].re*(*(s+0)).x - M[1][0].im*(*(s+0)).y ) + ( M[1][1].re*(*(s+0)).z - M[1][1].im*(*(s+0)).w ) + ( M[1][2].re*(*(s+1)).x - M[1][2].im*(*(s+1)).y );
(*(out+0)).w =  ( M[1][0].re*(*(s+0)).y + M[1][0].im*(*(s+0)).x ) + ( M[1][1].re*(*(s+0)).w + M[1][1].im*(*(s+0)).z ) + ( M[1][2].re*(*(s+1)).y + M[1][2].im*(*(s+1)).x );


(*(out+1)).x = ( M[2][0].re*(*(s+0)).x - M[2][0].im*(*(s+0)).y ) + ( M[2][1].re*(*(s+0)).z - M[2][1].im*(*(s+0)).w ) + ( M[2][2].re*(*(s+1)).x - M[2][2].im*(*(s+1)).y );
(*(out+1)).y =  ( M[2][0].re*(*(s+0)).y + M[2][0].im*(*(s+0)).x ) + ( M[2][1].re*(*(s+0)).w + M[2][1].im*(*(s+0)).z ) + ( M[2][2].re*(*(s+1)).y + M[2][2].im*(*(s+1)).x );


(*(out+1)).z =  ( M[0][0].re*(*(s+1)).z - M[0][0].im*(*(s+1)).w ) + ( M[0][1].re*(*(s+2)).x - M[0][1].im*(*(s+2)).y ) + ( M[0][2].re*(*(s+2)).z - M[0][2].im*(*(s+2)).w );
(*(out+1)).w =  ( M[0][0].re*(*(s+1)).w + M[0][0].im*(*(s+1)).z ) + ( M[0][1].re*(*(s+2)).y + M[0][1].im*(*(s+2)).x ) + ( M[0][2].re*(*(s+2)).w + M[0][2].im*(*(s+2)).z );


(*(out+2)).x = ( M[1][0].re*(*(s+1)).z - M[1][0].im*(*(s+1)).w ) + ( M[1][1].re*(*(s+2)).x - M[1][1].im*(*(s+2)).y ) + ( M[1][2].re*(*(s+2)).z - M[1][2].im*(*(s+2)).w );
(*(out+2)).y =  ( M[1][0].re*(*(s+1)).w + M[1][0].im*(*(s+1)).z ) + ( M[1][1].re*(*(s+2)).y + M[1][1].im*(*(s+2)).x ) + ( M[1][2].re*(*(s+2)).w + M[1][2].im*(*(s+2)).z );


(*(out+2)).z =  ( M[2][0].re*(*(s+1)).z - M[2][0].im*(*(s+1)).w ) + ( M[2][1].re*(*(s+2)).x - M[2][1].im*(*(s+2)).y ) + ( M[2][2].re*(*(s+2)).z - M[2][2].im*(*(s+2)).w );
(*(out+2)).w =  ( M[2][0].re*(*(s+1)).w + M[2][0].im*(*(s+1)).z ) + ( M[2][1].re*(*(s+2)).y + M[2][1].im*(*(s+2)).x ) + ( M[2][2].re*(*(s+2)).w + M[2][2].im*(*(s+2)).z );


(*(out+3)).x =  ( M[0][0].re*(*(s+3)).x - M[0][0].im*(*(s+3)).y ) + ( M[0][1].re*(*(s+3)).z - M[0][1].im*(*(s+3)).w ) + ( M[0][2].re*(*(s+4)).x - M[0][2].im*(*(s+4)).y );
(*(out+3)).y =   ( M[0][0].re*(*(s+3)).y + M[0][0].im*(*(s+3)).x ) + ( M[0][1].re*(*(s+3)).w + M[0][1].im*(*(s+3)).z ) + ( M[0][2].re*(*(s+4)).y + M[0][2].im*(*(s+4)).x );


(*(out+3)).z =  ( M[1][0].re*(*(s+3)).x - M[1][0].im*(*(s+3)).y ) + ( M[1][1].re*(*(s+3)).z - M[1][1].im*(*(s+3)).w ) + ( M[1][2].re*(*(s+4)).x - M[1][2].im*(*(s+4)).y );
(*(out+3)).w =  ( M[1][0].re*(*(s+3)).y + M[1][0].im*(*(s+3)).x ) + ( M[1][1].re*(*(s+3)).w + M[1][1].im*(*(s+3)).z ) + ( M[1][2].re*(*(s+4)).y + M[1][2].im*(*(s+4)).x );


(*(out+4)).x =  ( M[2][0].re*(*(s+3)).x - M[2][0].im*(*(s+3)).y ) + ( M[2][1].re*(*(s+3)).z - M[2][1].im*(*(s+3)).w ) + ( M[2][2].re*(*(s+4)).x - M[2][2].im*(*(s+4)).y );
(*(out+4)).y =  ( M[2][0].re*(*(s+3)).y + M[2][0].im*(*(s+3)).x ) + ( M[2][1].re*(*(s+3)).w + M[2][1].im*(*(s+3)).z ) + ( M[2][2].re*(*(s+4)).y + M[2][2].im*(*(s+4)).x );


(*(out+4)).z =  ( M[0][0].re*(*(s+4)).z - M[0][0].im*(*(s+4)).w ) + ( M[0][1].re*(*(s+5)).x - M[0][1].im*(*(s+5)).y ) + ( M[0][2].re*(*(s+5)).z - M[0][2].im*(*(s+5)).w );
(*(out+4)).w =   ( M[0][0].re*(*(s+4)).w + M[0][0].im*(*(s+4)).z ) + ( M[0][1].re*(*(s+5)).y + M[0][1].im*(*(s+5)).x ) + ( M[0][2].re*(*(s+5)).w + M[0][2].im*(*(s+5)).z );


(*(out+5)).x = ( M[1][0].re*(*(s+4)).z - M[1][0].im*(*(s+4)).w ) + ( M[1][1].re*(*(s+5)).x - M[1][1].im*(*(s+5)).y ) + ( M[1][2].re*(*(s+5)).z - M[1][2].im*(*(s+5)).w );
(*(out+5)).y =  ( M[1][0].re*(*(s+4)).w + M[1][0].im*(*(s+4)).z ) + ( M[1][1].re*(*(s+5)).y + M[1][1].im*(*(s+5)).x ) + ( M[1][2].re*(*(s+5)).w + M[1][2].im*(*(s+5)).z );


(*(out+5)).z =  ( M[2][0].re*(*(s+4)).z - M[2][0].im*(*(s+4)).w ) + ( M[2][1].re*(*(s+5)).x - M[2][1].im*(*(s+5)).y ) + ( M[2][2].re*(*(s+5)).z - M[2][2].im*(*(s+5)).w );
(*(out+5)).w =  ( M[2][0].re*(*(s+4)).w + M[2][0].im*(*(s+4)).z ) + ( M[2][1].re*(*(s+5)).y + M[2][1].im*(*(s+5)).x ) + ( M[2][2].re*(*(s+5)).w + M[2][2].im*(*(s+5)).z );
}






#ifdef HALF
//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_half(dev_su3 M, const dev_spinor_half * s, const float * s_norm, dev_spinor * out){
float norm = * s_norm;

(*(out+0)).x = ( M[0][0].re*half2fl((*(s+0)).x,norm) -
                 M[0][0].im*half2fl((*(s+0)).y,norm) ) 
             + ( M[0][1].re*half2fl((*(s+0)).z,norm) - 
                 M[0][1].im*half2fl((*(s+0)).w,norm) ) 
             + ( M[0][2].re*half2fl((*(s+1)).x,norm) -
                 M[0][2].im*half2fl((*(s+1)).y,norm) );
(*(out+0)).y = ( M[0][0].re*half2fl((*(s+0)).y,norm) +
                 M[0][0].im*half2fl((*(s+0)).x,norm) ) 
             + ( M[0][1].re*half2fl((*(s+0)).w,norm) +
                 M[0][1].im*half2fl((*(s+0)).z,norm) ) 
             + ( M[0][2].re*half2fl((*(s+1)).y,norm) +
                 M[0][2].im*half2fl((*(s+1)).x,norm) );


(*(out+0)).z =  ( M[1][0].re*half2fl((*(s+0)).x,norm) -
                  M[1][0].im*half2fl((*(s+0)).y,norm) ) 
             +  ( M[1][1].re*half2fl((*(s+0)).z,norm) - 
                  M[1][1].im*half2fl((*(s+0)).w,norm) ) 
             +  ( M[1][2].re*half2fl((*(s+1)).x,norm) -
                  M[1][2].im*half2fl((*(s+1)).y,norm) );
(*(out+0)).w =  ( M[1][0].re*half2fl((*(s+0)).y,norm) +
                  M[1][0].im*half2fl((*(s+0)).x,norm) ) 
             +  ( M[1][1].re*half2fl((*(s+0)).w,norm) +
                  M[1][1].im*half2fl((*(s+0)).z,norm) ) 
             +  ( M[1][2].re*half2fl((*(s+1)).y,norm) + 
                  M[1][2].im*half2fl((*(s+1)).x,norm) );


(*(out+1)).x = ( M[2][0].re*half2fl((*(s+0)).x, norm) -
                 M[2][0].im*half2fl((*(s+0)).y, norm) ) 
             + ( M[2][1].re*half2fl((*(s+0)).z, norm) -
                 M[2][1].im*half2fl((*(s+0)).w, norm) ) 
             + ( M[2][2].re*half2fl((*(s+1)).x, norm) -
                 M[2][2].im*half2fl((*(s+1)).y, norm) );
(*(out+1)).y = ( M[2][0].re*half2fl((*(s+0)).y, norm) +
                 M[2][0].im*half2fl((*(s+0)).x, norm) ) 
             + ( M[2][1].re*half2fl((*(s+0)).w, norm) +
                 M[2][1].im*half2fl((*(s+0)).z, norm) ) 
             + ( M[2][2].re*half2fl((*(s+1)).y, norm) +
                 M[2][2].im*half2fl((*(s+1)).x, norm) );


(*(out+1)).z = ( M[0][0].re*half2fl((*(s+1)).z, norm) -
                 M[0][0].im*half2fl((*(s+1)).w, norm) ) 
             + ( M[0][1].re*half2fl((*(s+2)).x, norm) -
                 M[0][1].im*half2fl((*(s+2)).y, norm) ) 
             + ( M[0][2].re*half2fl((*(s+2)).z, norm) -
                 M[0][2].im*half2fl((*(s+2)).w, norm) );
(*(out+1)).w = ( M[0][0].re*half2fl((*(s+1)).w, norm) +
                 M[0][0].im*half2fl((*(s+1)).z, norm) ) 
             + ( M[0][1].re*half2fl((*(s+2)).y, norm) +
                 M[0][1].im*half2fl((*(s+2)).x, norm) ) 
             + ( M[0][2].re*half2fl((*(s+2)).w, norm) +
                 M[0][2].im*half2fl((*(s+2)).z, norm) );


(*(out+2)).x = ( M[1][0].re*half2fl((*(s+1)).z, norm) -
                 M[1][0].im*half2fl((*(s+1)).w, norm) ) 
             + ( M[1][1].re*half2fl((*(s+2)).x, norm) -
                 M[1][1].im*half2fl((*(s+2)).y, norm) ) 
             + ( M[1][2].re*half2fl((*(s+2)).z, norm) -
                 M[1][2].im*half2fl((*(s+2)).w, norm) );
(*(out+2)).y = ( M[1][0].re*half2fl((*(s+1)).w, norm) +
                 M[1][0].im*half2fl((*(s+1)).z, norm) ) 
             + ( M[1][1].re*half2fl((*(s+2)).y, norm) +
                 M[1][1].im*half2fl((*(s+2)).x, norm) ) 
             + ( M[1][2].re*half2fl((*(s+2)).w, norm) +
                 M[1][2].im*half2fl((*(s+2)).z, norm) );


(*(out+2)).z = ( M[2][0].re*half2fl((*(s+1)).z, norm) -
                 M[2][0].im*half2fl((*(s+1)).w, norm) ) 
             + ( M[2][1].re*half2fl((*(s+2)).x, norm) -
                 M[2][1].im*half2fl((*(s+2)).y, norm) ) 
             + ( M[2][2].re*half2fl((*(s+2)).z, norm) -
                 M[2][2].im*half2fl((*(s+2)).w, norm) );
(*(out+2)).w = ( M[2][0].re*half2fl((*(s+1)).w, norm) +
                 M[2][0].im*half2fl((*(s+1)).z, norm) ) 
             + ( M[2][1].re*half2fl((*(s+2)).y, norm) +
                 M[2][1].im*half2fl((*(s+2)).x, norm) ) 
             + ( M[2][2].re*half2fl((*(s+2)).w, norm) +
                 M[2][2].im*half2fl((*(s+2)).z, norm) );


(*(out+3)).x = ( M[0][0].re*half2fl((*(s+3)).x, norm) -
                 M[0][0].im*half2fl((*(s+3)).y, norm) ) 
             + ( M[0][1].re*half2fl((*(s+3)).z, norm) -
                 M[0][1].im*half2fl((*(s+3)).w, norm) ) 
             + ( M[0][2].re*half2fl((*(s+4)).x, norm) -
                 M[0][2].im*half2fl((*(s+4)).y, norm) );
(*(out+3)).y = ( M[0][0].re*half2fl((*(s+3)).y, norm) +
                 M[0][0].im*half2fl((*(s+3)).x, norm) ) 
             + ( M[0][1].re*half2fl((*(s+3)).w, norm) +
                 M[0][1].im*half2fl((*(s+3)).z, norm) ) 
             + ( M[0][2].re*half2fl((*(s+4)).y, norm) +
                 M[0][2].im*half2fl((*(s+4)).x, norm) );


(*(out+3)).z = ( M[1][0].re*half2fl((*(s+3)).x, norm) -
                 M[1][0].im*half2fl((*(s+3)).y, norm) ) 
             + ( M[1][1].re*half2fl((*(s+3)).z, norm) -
                 M[1][1].im*half2fl((*(s+3)).w, norm) ) 
             + ( M[1][2].re*half2fl((*(s+4)).x, norm) -
                 M[1][2].im*half2fl((*(s+4)).y, norm) );
(*(out+3)).w = ( M[1][0].re*half2fl((*(s+3)).y, norm) +
                 M[1][0].im*half2fl((*(s+3)).x, norm) ) 
             + ( M[1][1].re*half2fl((*(s+3)).w, norm) +
                 M[1][1].im*half2fl((*(s+3)).z, norm) ) 
             + ( M[1][2].re*half2fl((*(s+4)).y, norm) +
                 M[1][2].im*half2fl((*(s+4)).x, norm) );


(*(out+4)).x = ( M[2][0].re*half2fl((*(s+3)).x, norm) -
                 M[2][0].im*half2fl((*(s+3)).y, norm) ) 
             + ( M[2][1].re*half2fl((*(s+3)).z, norm) -
                 M[2][1].im*half2fl((*(s+3)).w, norm) ) 
             + ( M[2][2].re*half2fl((*(s+4)).x, norm) -
                 M[2][2].im*half2fl((*(s+4)).y, norm) );
(*(out+4)).y = ( M[2][0].re*half2fl((*(s+3)).y, norm) +
                 M[2][0].im*half2fl((*(s+3)).x, norm) ) 
             + ( M[2][1].re*half2fl((*(s+3)).w, norm) +
                 M[2][1].im*half2fl((*(s+3)).z, norm) ) 
             + ( M[2][2].re*half2fl((*(s+4)).y, norm) +
                 M[2][2].im*half2fl((*(s+4)).x, norm) );


(*(out+4)).z = ( M[0][0].re*half2fl((*(s+4)).z, norm) -
                 M[0][0].im*half2fl((*(s+4)).w, norm) ) 
             + ( M[0][1].re*half2fl((*(s+5)).x, norm) - 
                 M[0][1].im*half2fl((*(s+5)).y, norm) ) 
             + ( M[0][2].re*half2fl((*(s+5)).z, norm) -
                 M[0][2].im*half2fl((*(s+5)).w, norm) );
(*(out+4)).w = ( M[0][0].re*half2fl((*(s+4)).w, norm) +
                 M[0][0].im*half2fl((*(s+4)).z, norm) ) 
             + ( M[0][1].re*half2fl((*(s+5)).y, norm) +
                 M[0][1].im*half2fl((*(s+5)).x, norm) ) 
             + ( M[0][2].re*half2fl((*(s+5)).w, norm) +
                 M[0][2].im*half2fl((*(s+5)).z, norm) );


(*(out+5)).x = ( M[1][0].re*half2fl((*(s+4)).z, norm) -
                 M[1][0].im*half2fl((*(s+4)).w, norm) ) 
             + ( M[1][1].re*half2fl((*(s+5)).x, norm) -
                 M[1][1].im*half2fl((*(s+5)).y, norm) ) 
             + ( M[1][2].re*half2fl((*(s+5)).z, norm) -
                 M[1][2].im*half2fl((*(s+5)).w, norm) );
(*(out+5)).y = ( M[1][0].re*half2fl((*(s+4)).w, norm) +
                 M[1][0].im*half2fl((*(s+4)).z, norm) ) 
             + ( M[1][1].re*half2fl((*(s+5)).y, norm) +
                 M[1][1].im*half2fl((*(s+5)).x, norm) ) 
             + ( M[1][2].re*half2fl((*(s+5)).w, norm) +
                 M[1][2].im*half2fl((*(s+5)).z, norm) );


(*(out+5)).z = ( M[2][0].re*half2fl((*(s+4)).z, norm) - 
                 M[2][0].im*half2fl((*(s+4)).w, norm) ) 
             + ( M[2][1].re*half2fl((*(s+5)).x, norm) -
                 M[2][1].im*half2fl((*(s+5)).y, norm) ) 
             + ( M[2][2].re*half2fl((*(s+5)).z, norm) -
                 M[2][2].im*half2fl((*(s+5)).w, norm) );
(*(out+5)).w = ( M[2][0].re*half2fl((*(s+4)).w, norm) +
                 M[2][0].im*half2fl((*(s+4)).z, norm) ) 
             + ( M[2][1].re*half2fl((*(s+5)).y, norm) +
                 M[2][1].im*half2fl((*(s+5)).x, norm) ) 
             + ( M[2][2].re*half2fl((*(s+5)).w, norm) +
                 M[2][2].im*half2fl((*(s+5)).z, norm) );
}
#endif







//multipliziert gedaggerte su3-Matrix mal Spinor im Dirac-Raum  -- generated with codegen
template<class RealT>
__device__ void dev_su3MdaggertV(typename dev_su3T<RealT>::type M, typename dev_spinorT<RealT>::type * s, typename dev_spinorT<RealT>::type * out){
  dev_complexT<RealT> help1;
help1.re = M[0][0].re*(*(s+0)).x + M[0][0].im*(*(s+0)).y + M[1][0].re*(*(s+0)).z + M[1][0].im*(*(s+0)).w + M[2][0].re*(*(s+1)).x + M[2][0].im*(*(s+1)).y;
(*(out+0)).x = help1.re;
help1.im = M[0][0].re*(*(s+0)).y - M[0][0].im*(*(s+0)).x + M[1][0].re*(*(s+0)).w - M[1][0].im*(*(s+0)).z + M[2][0].re*(*(s+1)).y - M[2][0].im*(*(s+1)).x;
(*(out+0)).y = help1.im;

help1.re = M[0][1].re*(*(s+0)).x + M[0][1].im*(*(s+0)).y + M[1][1].re*(*(s+0)).z + M[1][1].im*(*(s+0)).w + M[2][1].re*(*(s+1)).x + M[2][1].im*(*(s+1)).y;
(*(out+0)).z = help1.re;
help1.im = M[0][1].re*(*(s+0)).y - M[0][1].im*(*(s+0)).x + M[1][1].re*(*(s+0)).w - M[1][1].im*(*(s+0)).z + M[2][1].re*(*(s+1)).y - M[2][1].im*(*(s+1)).x;
(*(out+0)).w = help1.im;

help1.re = M[0][2].re*(*(s+0)).x + M[0][2].im*(*(s+0)).y + M[1][2].re*(*(s+0)).z + M[1][2].im*(*(s+0)).w + M[2][2].re*(*(s+1)).x + M[2][2].im*(*(s+1)).y;
(*(out+1)).x = help1.re;
help1.im = M[0][2].re*(*(s+0)).y - M[0][2].im*(*(s+0)).x + M[1][2].re*(*(s+0)).w - M[1][2].im*(*(s+0)).z + M[2][2].re*(*(s+1)).y - M[2][2].im*(*(s+1)).x;
(*(out+1)).y = help1.im;

help1.re = M[0][0].re*(*(s+1)).z + M[0][0].im*(*(s+1)).w + M[1][0].re*(*(s+2)).x + M[1][0].im*(*(s+2)).y + M[2][0].re*(*(s+2)).z + M[2][0].im*(*(s+2)).w;
(*(out+1)).z = help1.re;
help1.im = M[0][0].re*(*(s+1)).w - M[0][0].im*(*(s+1)).z + M[1][0].re*(*(s+2)).y - M[1][0].im*(*(s+2)).x + M[2][0].re*(*(s+2)).w - M[2][0].im*(*(s+2)).z;
(*(out+1)).w = help1.im;

help1.re = M[0][1].re*(*(s+1)).z + M[0][1].im*(*(s+1)).w + M[1][1].re*(*(s+2)).x + M[1][1].im*(*(s+2)).y + M[2][1].re*(*(s+2)).z + M[2][1].im*(*(s+2)).w;
(*(out+2)).x = help1.re;
help1.im = M[0][1].re*(*(s+1)).w - M[0][1].im*(*(s+1)).z + M[1][1].re*(*(s+2)).y - M[1][1].im*(*(s+2)).x + M[2][1].re*(*(s+2)).w - M[2][1].im*(*(s+2)).z;
(*(out+2)).y = help1.im;

help1.re = M[0][2].re*(*(s+1)).z + M[0][2].im*(*(s+1)).w + M[1][2].re*(*(s+2)).x + M[1][2].im*(*(s+2)).y + M[2][2].re*(*(s+2)).z + M[2][2].im*(*(s+2)).w;
(*(out+2)).z = help1.re;
help1.im = M[0][2].re*(*(s+1)).w - M[0][2].im*(*(s+1)).z + M[1][2].re*(*(s+2)).y - M[1][2].im*(*(s+2)).x + M[2][2].re*(*(s+2)).w - M[2][2].im*(*(s+2)).z;
(*(out+2)).w = help1.im;

help1.re = M[0][0].re*(*(s+3)).x + M[0][0].im*(*(s+3)).y + M[1][0].re*(*(s+3)).z + M[1][0].im*(*(s+3)).w + M[2][0].re*(*(s+4)).x + M[2][0].im*(*(s+4)).y;
(*(out+3)).x = help1.re;
help1.im = M[0][0].re*(*(s+3)).y - M[0][0].im*(*(s+3)).x + M[1][0].re*(*(s+3)).w - M[1][0].im*(*(s+3)).z + M[2][0].re*(*(s+4)).y - M[2][0].im*(*(s+4)).x;
(*(out+3)).y = help1.im;

help1.re = M[0][1].re*(*(s+3)).x + M[0][1].im*(*(s+3)).y + M[1][1].re*(*(s+3)).z + M[1][1].im*(*(s+3)).w + M[2][1].re*(*(s+4)).x + M[2][1].im*(*(s+4)).y;
(*(out+3)).z = help1.re;
help1.im = M[0][1].re*(*(s+3)).y - M[0][1].im*(*(s+3)).x + M[1][1].re*(*(s+3)).w - M[1][1].im*(*(s+3)).z + M[2][1].re*(*(s+4)).y - M[2][1].im*(*(s+4)).x;
(*(out+3)).w = help1.im;

help1.re = M[0][2].re*(*(s+3)).x + M[0][2].im*(*(s+3)).y + M[1][2].re*(*(s+3)).z + M[1][2].im*(*(s+3)).w + M[2][2].re*(*(s+4)).x + M[2][2].im*(*(s+4)).y;
(*(out+4)).x = help1.re;
help1.im = M[0][2].re*(*(s+3)).y - M[0][2].im*(*(s+3)).x + M[1][2].re*(*(s+3)).w - M[1][2].im*(*(s+3)).z + M[2][2].re*(*(s+4)).y - M[2][2].im*(*(s+4)).x;
(*(out+4)).y = help1.im;

help1.re = M[0][0].re*(*(s+4)).z + M[0][0].im*(*(s+4)).w + M[1][0].re*(*(s+5)).x + M[1][0].im*(*(s+5)).y + M[2][0].re*(*(s+5)).z + M[2][0].im*(*(s+5)).w;
(*(out+4)).z = help1.re;
help1.im = M[0][0].re*(*(s+4)).w - M[0][0].im*(*(s+4)).z + M[1][0].re*(*(s+5)).y - M[1][0].im*(*(s+5)).x + M[2][0].re*(*(s+5)).w - M[2][0].im*(*(s+5)).z;
(*(out+4)).w = help1.im;

help1.re = M[0][1].re*(*(s+4)).z + M[0][1].im*(*(s+4)).w + M[1][1].re*(*(s+5)).x + M[1][1].im*(*(s+5)).y + M[2][1].re*(*(s+5)).z + M[2][1].im*(*(s+5)).w;
(*(out+5)).x = help1.re;
help1.im = M[0][1].re*(*(s+4)).w - M[0][1].im*(*(s+4)).z + M[1][1].re*(*(s+5)).y - M[1][1].im*(*(s+5)).x + M[2][1].re*(*(s+5)).w - M[2][1].im*(*(s+5)).z;
(*(out+5)).y = help1.im;

help1.re = M[0][2].re*(*(s+4)).z + M[0][2].im*(*(s+4)).w + M[1][2].re*(*(s+5)).x + M[1][2].im*(*(s+5)).y + M[2][2].re*(*(s+5)).z + M[2][2].im*(*(s+5)).w;
(*(out+5)).z = help1.re;
help1.im = M[0][2].re*(*(s+4)).w - M[0][2].im*(*(s+4)).z + M[1][2].re*(*(s+5)).y - M[1][2].im*(*(s+5)).x + M[2][2].re*(*(s+5)).w - M[2][2].im*(*(s+5)).z;
(*(out+5)).w = help1.im;
}




// Gamma t
template<class RealT>
__device__ void dev_Gamma0(typename dev_spinorT<RealT>::type * in){
  RealT tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+3)).x;
     (*(in+0)).y = -1.0*(*(in+3)).y;
     (*(in+3)).x = -1.0*tempre;
     (*(in+3)).y = -1.0*tempim;     
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+3)).z;
     (*(in+0)).w = -1.0*(*(in+3)).w;
     (*(in+3)).z = -1.0*tempre;
     (*(in+3)).w = -1.0*tempim; 
 
 
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+4)).x;
     (*(in+1)).y = -1.0*(*(in+4)).y;
     (*(in+4)).x = -1.0*tempre;
     (*(in+4)).y = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).z;
     (*(in+1)).w = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).x;
     (*(in+2)).y = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
   
   
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).z;
     (*(in+2)).w = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;
}



//Gamma z
template<class RealT>
__device__ void dev_Gamma3(typename dev_spinorT<RealT>::type * in){
  RealT tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+3)).y;
     (*(in+0)).y = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+3)).w;
     (*(in+0)).w = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;    
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+4)).y;
     (*(in+1)).y = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;     
     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).w;
     (*(in+1)).w = (*(in+4)).z;
     (*(in+4)).z  = tempim;
     (*(in+4)).w  = -1.0*tempre;     
     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).y;
     (*(in+2)).y = (*(in+5)).x;
     (*(in+5)).x = tempim;
     (*(in+5)).y = -1.0*tempre;    
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).w;
     (*(in+2)).w = (*(in+5)).z;
     (*(in+5)).z = tempim;
     (*(in+5)).w = -1.0*tempre;

}



//Gamma y
template<class RealT>
__device__ void dev_Gamma2(typename dev_spinorT<RealT>::type * in){
  RealT tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+4)).z;
     (*(in+0)).y = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+5)).x;
     (*(in+0)).w = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+5)).z;
     (*(in+1)).y = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).x;
     (*(in+1)).w = (*(in+3)).y;
     (*(in+3)).x = tempre;
     (*(in+3)).y = tempim;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).z;
     (*(in+2)).y = (*(in+3)).w;
     (*(in+3)).z = tempre;
     (*(in+3)).w = tempim;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).x;
     (*(in+2)).w = (*(in+4)).y;
     (*(in+4)).x = tempre;
     (*(in+4)).y = tempim;
}



//Gamma x
template<class RealT>
__device__ void dev_Gamma1(typename dev_spinorT<RealT>::type * in){
  RealT tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+4)).w;
     (*(in+0)).y = -1.0*(*(in+4)).z;
     (*(in+4)).z  = -1.0*tempim;
     (*(in+4)).w  = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+5)).y;
     (*(in+0)).w = -1.0*(*(in+5)).x;
     (*(in+5)).x = -1.0*tempim;
     (*(in+5)).y = tempre;     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+5)).w;
     (*(in+1)).y = -1.0*(*(in+5)).z;
     (*(in+5)).z = -1.0*tempim;
     (*(in+5)).w = tempre;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).y;
     (*(in+1)).w = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).w;
     (*(in+2)).y = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).y;
     (*(in+2)).w = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;
  
}



template<class RealT>
__device__ void dev_Gamma5(typename dev_spinorT<RealT>::type * in){
          (*(in+3)).x = -1.0*(*(in+3)).x;
          (*(in+3)).y = -1.0*(*(in+3)).y;
          (*(in+3)).z = -1.0*(*(in+3)).z;
          (*(in+3)).w = -1.0*(*(in+3)).w;
          (*(in+4)).x = -1.0*(*(in+4)).x;
          (*(in+4)).y = -1.0*(*(in+4)).y; 

          (*(in+4)).z = -1.0*(*(in+4)).z;
          (*(in+4)).w = -1.0*(*(in+4)).w;
          (*(in+5)).x = -1.0*(*(in+5)).x;
          (*(in+5)).y = -1.0*(*(in+5)).y;
          (*(in+5)).z = -1.0*(*(in+5)).z;
          (*(in+5)).w = -1.0*(*(in+5)).w;  
}


template<class RealT>
__device__ void dev_Gamma5_assign(typename dev_spinorT<RealT>::type* out, typename dev_spinorT<RealT>::type* in){
  (*(out)).x = (*(in)).x;
  (*(out)).y = (*(in)).y;
  (*(out)).z = (*(in)).z;
  (*(out)).w = (*(in)).w;
  (*(out+1)).x = (*(in+1)).x;
  (*(out+1)).y = (*(in+1)).y;

  (*(out+1)).z = (*(in+1)).z;
  (*(out+1)).w = (*(in+1)).w;
  (*(out+2)).x = (*(in+2)).x;
  (*(out+2)).y = (*(in+2)).y;
  (*(out+2)).z = (*(in+2)).z;
  (*(out+2)).w = (*(in+2)).w;

  (*(out+3)).x = -1.0*(*(in+3)).x;
  (*(out+3)).y = -1.0*(*(in+3)).y;
  (*(out+3)).z = -1.0*(*(in+3)).z;
  (*(out+3)).w = -1.0*(*(in+3)).w;
  (*(out+4)).x = -1.0*(*(in+4)).x;
  (*(out+4)).y = -1.0*(*(in+4)).y;

  (*(out+4)).z = -1.0*(*(in+4)).z;
  (*(out+4)).w = -1.0*(*(in+4)).w;
  (*(out+5)).x = -1.0*(*(in+5)).x;
  (*(out+5)).y = -1.0*(*(in+5)).y;
  (*(out+5)).z = -1.0*(*(in+5)).z;
  (*(out+5)).w = -1.0*(*(in+5)).w;
}




// older version, all in one function
template<class RealT>
__device__ void dev_GammatV(int mu, typename dev_spinorT<RealT>::type * in){//multipliziert Gamma(mu)*V effizientes ausnutzen der Nullen 
 RealT tempre,tempim;
 /* ORDER: t, z, y, x*/
 switch (mu){
 
 case 0:
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+3)).x;
     (*(in+0)).y = -1.0*(*(in+3)).y;
     (*(in+3)).x = -1.0*tempre;
     (*(in+3)).y = -1.0*tempim;     
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+3)).z;
     (*(in+0)).w = -1.0*(*(in+3)).w;
     (*(in+3)).z = -1.0*tempre;
     (*(in+3)).w = -1.0*tempim; 
 
 
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+4)).x;
     (*(in+1)).y = -1.0*(*(in+4)).y;
     (*(in+4)).x = -1.0*tempre;
     (*(in+4)).y = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).z;
     (*(in+1)).w = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).x;
     (*(in+2)).y = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
   
   
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).z;
     (*(in+2)).w = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;

 break;
 
 case 1:
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+3)).y;
     (*(in+0)).y = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+3)).w;
     (*(in+0)).w = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;    
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+4)).y;
     (*(in+1)).y = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;     
     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).w;
     (*(in+1)).w = (*(in+4)).z;
     (*(in+4)).z  = tempim;
     (*(in+4)).w  = -1.0*tempre;     
     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).y;
     (*(in+2)).y = (*(in+5)).x;
     (*(in+5)).x = tempim;
     (*(in+5)).y = -1.0*tempre;    
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).w;
     (*(in+2)).w = (*(in+5)).z;
     (*(in+5)).z = tempim;
     (*(in+5)).w = -1.0*tempre;


 break;
 
 case 2:
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+4)).z;
     (*(in+0)).y = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+5)).x;
     (*(in+0)).w = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+5)).z;
     (*(in+1)).y = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).x;
     (*(in+1)).w = (*(in+3)).y;
     (*(in+3)).x = tempre;
     (*(in+3)).y = tempim;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).z;
     (*(in+2)).y = (*(in+3)).w;
     (*(in+3)).z = tempre;
     (*(in+3)).w = tempim;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).x;
     (*(in+2)).w = (*(in+4)).y;
     (*(in+4)).x = tempre;
     (*(in+4)).y = tempim;

 break; 
 
 case 3:


     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+4)).w;
     (*(in+0)).y = -1.0*(*(in+4)).z;
     (*(in+4)).z  = -1.0*tempim;
     (*(in+4)).w  = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+5)).y;
     (*(in+0)).w = -1.0*(*(in+5)).x;
     (*(in+5)).x = -1.0*tempim;
     (*(in+5)).y = tempre;     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+5)).w;
     (*(in+1)).y = -1.0*(*(in+5)).z;
     (*(in+5)).z = -1.0*tempim;
     (*(in+5)).w = tempre;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).y;
     (*(in+1)).w = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).w;
     (*(in+2)).y = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).y;
     (*(in+2)).w = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;
     
     
 break;
 
 
 case 4:
  
          (*(in+3)).x = -1.0*(*(in+3)).x;
          (*(in+3)).y = -1.0*(*(in+3)).y;
          (*(in+3)).z = -1.0*(*(in+3)).z;
          (*(in+3)).w = -1.0*(*(in+3)).w;
          (*(in+4)).x = -1.0*(*(in+4)).x;
          (*(in+4)).y = -1.0*(*(in+4)).y; 

          (*(in+4)).z = -1.0*(*(in+4)).z;
          (*(in+4)).w = -1.0*(*(in+4)).w;
          (*(in+5)).x = -1.0*(*(in+5)).x;
          (*(in+5)).y = -1.0*(*(in+5)).y;
          (*(in+5)).z = -1.0*(*(in+5)).z;
          (*(in+5)).w = -1.0*(*(in+5)).w;  
 break;
 }
}



