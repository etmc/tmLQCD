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
 * File: Hopping_Matrix.cuh
 *
 * CUDA Hopping_Matrix and associated functions
 *
 * 
 *
 **************************************************************************/






//-kappa(r - gamma_mu)
__device__ void dev_kappaP1_plus(dev_spinor * out, dev_spinor * in, float kappa){

     (*(out+0)).x -= kappa*( (*(in+0)).x - (*(in+4)).w);
     (*(out+0)).y -= kappa*( (*(in+0)).y + (*(in+4)).z);
     (*(out+0)).z -= kappa*( (*(in+0)).z - (*(in+5)).y);
     (*(out+0)).w -= kappa*( (*(in+0)).w + (*(in+5)).x);    

     (*(out+1)).x -= kappa*((*(in+1)).x - (*(in+5)).w);
     (*(out+1)).y -= kappa*((*(in+1)).y + (*(in+5)).z); 
     (*(out+1)).z -= kappa*((*(in+1)).z - (*(in+3)).y);
     (*(out+1)).w -= kappa*((*(in+1)).w + (*(in+3)).x); 
     
     (*(out+2)).x -= kappa*((*(in+2)).x - (*(in+3)).w);
     (*(out+2)).y -= kappa*((*(in+2)).y + (*(in+3)).z);
     (*(out+2)).z -= kappa*((*(in+2)).z - (*(in+4)).y);
     (*(out+2)).w -= kappa*((*(in+2)).w + (*(in+4)).x);     
     
     (*(out+3)).x -= kappa*((*(in+3)).x + (*(in+1)).w);
     (*(out+3)).y -= kappa*((*(in+3)).y - (*(in+1)).z);     
     (*(out+3)).z -= kappa*((*(in+3)).z + (*(in+2)).y);
     (*(out+3)).w -= kappa*((*(in+3)).w - (*(in+2)).x);       
     
     (*(out+4)).z -= kappa*( (*(in+4)).z + (*(in+0)).y);
     (*(out+4)).w -= kappa*( (*(in+4)).w - (*(in+0)).x);    
     (*(out+4)).x -= kappa*((*(in+4)).x + (*(in+2)).w);
     (*(out+4)).y -= kappa*((*(in+4)).y - (*(in+2)).z);     

     (*(out+5)).x -= kappa*( (*(in+5)).x + (*(in+0)).w);
     (*(out+5)).y -= kappa*( (*(in+5)).y - (*(in+0)).z);     
     (*(out+5)).z -= kappa*((*(in+5)).z + (*(in+1)).y);
     (*(out+5)).w -= kappa*((*(in+5)).w - (*(in+1)).x);     
     
}


//-kappa(r + gamma_mu)
__device__ void dev_kappaP1_minus(dev_spinor * out, dev_spinor * in, float kappa){

     (*(out+0)).x -= kappa*( (*(in+0)).x + (*(in+4)).w);
     (*(out+0)).y -= kappa*( (*(in+0)).y - (*(in+4)).z);
     (*(out+0)).z -= kappa*( (*(in+0)).z + (*(in+5)).y);
     (*(out+0)).w -= kappa*( (*(in+0)).w - (*(in+5)).x);    

     (*(out+1)).x -= kappa*((*(in+1)).x + (*(in+5)).w);
     (*(out+1)).y -= kappa*((*(in+1)).y - (*(in+5)).z); 
     (*(out+1)).z -= kappa*((*(in+1)).z + (*(in+3)).y);
     (*(out+1)).w -= kappa*((*(in+1)).w - (*(in+3)).x); 
     
     (*(out+2)).x -= kappa*((*(in+2)).x + (*(in+3)).w);
     (*(out+2)).y -= kappa*((*(in+2)).y - (*(in+3)).z);
     (*(out+2)).z -= kappa*((*(in+2)).z + (*(in+4)).y);
     (*(out+2)).w -= kappa*((*(in+2)).w - (*(in+4)).x);     
     
     (*(out+3)).x -= kappa*((*(in+3)).x - (*(in+1)).w);
     (*(out+3)).y -= kappa*((*(in+3)).y + (*(in+1)).z);     
     (*(out+3)).z -= kappa*((*(in+3)).z - (*(in+2)).y);
     (*(out+3)).w -= kappa*((*(in+3)).w + (*(in+2)).x);       
     
     (*(out+4)).z -= kappa*( (*(in+4)).z - (*(in+0)).y);
     (*(out+4)).w -= kappa*( (*(in+4)).w + (*(in+0)).x);    
     (*(out+4)).x -= kappa*((*(in+4)).x - (*(in+2)).w);
     (*(out+4)).y -= kappa*((*(in+4)).y + (*(in+2)).z);     

     (*(out+5)).x -= kappa*( (*(in+5)).x - (*(in+0)).w);
     (*(out+5)).y -= kappa*( (*(in+5)).y + (*(in+0)).z);     
     (*(out+5)).z -= kappa*((*(in+5)).z - (*(in+1)).y);
     (*(out+5)).w -= kappa*((*(in+5)).w + (*(in+1)).x);     
     
}





//-kappa(r - gamma_mu)
__device__ void dev_kappaP2_plus(dev_spinor * out, dev_spinor * in, float kappa){


     (*(out+0)).x -= kappa*( (*(in+0)).x + (*(in+4)).z);
     (*(out+0)).y -= kappa*( (*(in+0)).y + (*(in+4)).w);
     (*(out+4)).z -= kappa*( (*(in+4)).z + (*(in+0)).x);
     (*(out+4)).w -= kappa*( (*(in+4)).w + (*(in+0)).y);    
     
 
     (*(out+0)).z -= kappa*( (*(in+0)).z + (*(in+5)).x);
     (*(out+0)).w -= kappa*( (*(in+0)).w + (*(in+5)).y);
     (*(out+5)).x -= kappa*( (*(in+5)).x + (*(in+0)).z);
     (*(out+5)).y -= kappa*( (*(in+5)).y + (*(in+0)).w);     
     
     
     (*(out+1)).x -= kappa*( (*(in+1)).x + (*(in+5)).z);
     (*(out+1)).y -= kappa*( (*(in+1)).y + (*(in+5)).w);
     (*(out+5)).z -= kappa*( (*(in+5)).z + (*(in+1)).x);
     (*(out+5)).w -= kappa*( (*(in+5)).w + (*(in+1)).y);     
     
     
     (*(out+1)).z -= kappa*( (*(in+1)).z - (*(in+3)).x);
     (*(out+1)).w -= kappa*( (*(in+1)).w - (*(in+3)).y);
     (*(out+3)).x -= kappa*( (*(in+3)).x - (*(in+1)).z);
     (*(out+3)).y -= kappa*( (*(in+3)).y - (*(in+1)).w);     
     
     
     (*(out+2)).x -= kappa*( (*(in+2)).x - (*(in+3)).z);
     (*(out+2)).y -= kappa*( (*(in+2)).y - (*(in+3)).w);
     (*(out+3)).z -= kappa*( (*(in+3)).z - (*(in+2)).x);
     (*(out+3)).w -= kappa*( (*(in+3)).w - (*(in+2)).y);     
     
     
     (*(out+2)).z -= kappa*( (*(in+2)).z - (*(in+4)).x);
     (*(out+2)).w -= kappa*( (*(in+2)).w - (*(in+4)).y);
     (*(out+4)).x -= kappa*( (*(in+4)).x - (*(in+2)).z);
     (*(out+4)).y -= kappa*( (*(in+4)).y - (*(in+2)).w);
   
     
}


//-kappa(r + gamma_mu)  kappa reell !!!!
__device__ void dev_kappaP2_minus(dev_spinor * out, dev_spinor * in, float kappa){


     (*(out+0)).x -= kappa*( (*(in+0)).x - (*(in+4)).z);
     (*(out+0)).y -= kappa*( (*(in+0)).y - (*(in+4)).w);
     (*(out+4)).z -= kappa*( (*(in+4)).z - (*(in+0)).x);
     (*(out+4)).w -= kappa*( (*(in+4)).w - (*(in+0)).y);    
     
 
     (*(out+0)).z -= kappa*( (*(in+0)).z - (*(in+5)).x);
     (*(out+0)).w -= kappa*( (*(in+0)).w - (*(in+5)).y);
     (*(out+5)).x -= kappa*( (*(in+5)).x - (*(in+0)).z);
     (*(out+5)).y -= kappa*( (*(in+5)).y - (*(in+0)).w);     
     
     
     (*(out+1)).x -= kappa*( (*(in+1)).x - (*(in+5)).z);
     (*(out+1)).y -= kappa*( (*(in+1)).y - (*(in+5)).w);
     (*(out+5)).z -= kappa*( (*(in+5)).z - (*(in+1)).x);
     (*(out+5)).w -= kappa*( (*(in+5)).w - (*(in+1)).y);     
     
     
     (*(out+1)).z -= kappa*( (*(in+1)).z + (*(in+3)).x);
     (*(out+1)).w -= kappa*( (*(in+1)).w + (*(in+3)).y);
     (*(out+3)).x -= kappa*( (*(in+3)).x + (*(in+1)).z);
     (*(out+3)).y -= kappa*( (*(in+3)).y + (*(in+1)).w);     
     
     
     (*(out+2)).x -= kappa*( (*(in+2)).x + (*(in+3)).z);
     (*(out+2)).y -= kappa*( (*(in+2)).y + (*(in+3)).w);
     (*(out+3)).z -= kappa*( (*(in+3)).z + (*(in+2)).x);
     (*(out+3)).w -= kappa*( (*(in+3)).w + (*(in+2)).y);     
     
     
     (*(out+2)).z -= kappa*( (*(in+2)).z + (*(in+4)).x);
     (*(out+2)).w -= kappa*( (*(in+2)).w + (*(in+4)).y);
     (*(out+4)).x -= kappa*( (*(in+4)).x + (*(in+2)).z);
     (*(out+4)).y -= kappa*( (*(in+4)).y + (*(in+2)).w);
   
     
}



//-kappa(r - gamma_mu) kappa reell !!!!
__device__ void dev_kappaP3_plus(dev_spinor * out, dev_spinor * in, float kappa){

     (*(out+0)).x -= kappa*( (*(in+0)).x - (*(in+3)).y);
     (*(out+0)).y -= kappa*( (*(in+0)).y + (*(in+3)).x);
     (*(out+3)).x -= kappa*( (*(in+3)).x + (*(in+0)).y);
     (*(out+3)).y -= kappa*( (*(in+3)).y - (*(in+0)).x);    
     

     (*(out+0)).z -= kappa*( (*(in+0)).z - (*(in+3)).w);
     (*(out+0)).w -= kappa*( (*(in+0)).w + (*(in+3)).z);
     (*(out+3)).z -= kappa*( (*(in+3)).z + (*(in+0)).w);
     (*(out+3)).w -= kappa*( (*(in+3)).w - (*(in+0)).z);    
     
     
     (*(out+1)).x -= kappa*( (*(in+1)).x - (*(in+4)).y);
     (*(out+1)).y -= kappa*( (*(in+1)).y + (*(in+4)).x);
     (*(out+4)).x -= kappa*( (*(in+4)).x + (*(in+1)).y);
     (*(out+4)).y -= kappa*( (*(in+4)).y - (*(in+1)).x);     
     
     
     (*(out+1)).z -= kappa*( (*(in+1)).z + (*(in+4)).w);
     (*(out+1)).w -= kappa*( (*(in+1)).w - (*(in+4)).z);
     (*(out+4)).z -= kappa*( (*(in+4)).z - (*(in+1)).w);
     (*(out+4)).w -= kappa*( (*(in+4)).w + (*(in+1)).z);     
     
       
     (*(out+2)).x -= kappa*( (*(in+2)).x + (*(in+5)).y);
     (*(out+2)).y -= kappa*( (*(in+2)).y - (*(in+5)).x);
     (*(out+5)).x -= kappa*( (*(in+5)).x - (*(in+2)).y);
     (*(out+5)).y -= kappa*( (*(in+5)).y + (*(in+2)).x);    
     

     (*(out+2)).z -= kappa*( (*(in+2)).z + (*(in+5)).w);
     (*(out+2)).w -= kappa*( (*(in+2)).w - (*(in+5)).z);
     (*(out+5)).z -= kappa*( (*(in+5)).z - (*(in+2)).w);
     (*(out+5)).w -= kappa*( (*(in+5)).w + (*(in+2)).z);
  
}


//-kappa(r + gamma_mu) kappa reell !!!
__device__ void dev_kappaP3_minus(dev_spinor * out, dev_spinor * in, float kappa){

     (*(out+0)).x -= kappa*( (*(in+0)).x + (*(in+3)).y);
     (*(out+0)).y -= kappa*( (*(in+0)).y - (*(in+3)).x);
     (*(out+3)).x -= kappa*( (*(in+3)).x - (*(in+0)).y);
     (*(out+3)).y -= kappa*( (*(in+3)).y + (*(in+0)).x);    
     

     (*(out+0)).z -= kappa*( (*(in+0)).z + (*(in+3)).w);
     (*(out+0)).w -= kappa*( (*(in+0)).w - (*(in+3)).z);
     (*(out+3)).z -= kappa*( (*(in+3)).z - (*(in+0)).w);
     (*(out+3)).w -= kappa*( (*(in+3)).w + (*(in+0)).z);    
     
     
     (*(out+1)).x -= kappa*( (*(in+1)).x + (*(in+4)).y);
     (*(out+1)).y -= kappa*( (*(in+1)).y - (*(in+4)).x);
     (*(out+4)).x -= kappa*( (*(in+4)).x - (*(in+1)).y);
     (*(out+4)).y -= kappa*( (*(in+4)).y + (*(in+1)).x);     
     
     
     (*(out+1)).z -= kappa*( (*(in+1)).z - (*(in+4)).w);
     (*(out+1)).w -= kappa*( (*(in+1)).w + (*(in+4)).z);
     (*(out+4)).z -= kappa*( (*(in+4)).z + (*(in+1)).w);
     (*(out+4)).w -= kappa*( (*(in+4)).w - (*(in+1)).z);     
     
       
     (*(out+2)).x -= kappa*( (*(in+2)).x - (*(in+5)).y);
     (*(out+2)).y -= kappa*( (*(in+2)).y + (*(in+5)).x);
     (*(out+5)).x -= kappa*( (*(in+5)).x + (*(in+2)).y);
     (*(out+5)).y -= kappa*( (*(in+5)).y - (*(in+2)).x);    
     

     (*(out+2)).z -= kappa*( (*(in+2)).z - (*(in+5)).w);
     (*(out+2)).w -= kappa*( (*(in+2)).w + (*(in+5)).z);
     (*(out+5)).z -= kappa*( (*(in+5)).z + (*(in+2)).w);
     (*(out+5)).w -= kappa*( (*(in+5)).w - (*(in+2)).z);
  
}







//-kappa(r - gamma_mu)
__device__ void dev_kappaP0_plus(dev_spinor * out, dev_spinor * in, dev_complex kappa){


     (*(out+0)).x -= (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+0)).y -= (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im;
     (*(out+0)).x -= (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+0)).y -= (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;
     
     (*(out+3)).x -= (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+3)).y -= (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;   
     (*(out+3)).x -= (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+3)).y -= (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im; 



     (*(out+0)).z -= (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+0)).w -= (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;
     (*(out+0)).z -= (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+0)).w -= (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;     
          
     (*(out+3)).z -= (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+3)).w -= (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;
     (*(out+3)).z -= (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+3)).w -= (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;

 
  
     (*(out+1)).x -= (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+1)).y -= (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im;
     (*(out+1)).x -= (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+1)).y -= (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;     
     
     (*(out+4)).x -= (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+4)).y -= (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;    
     (*(out+4)).x -= (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+4)).y -= (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im; 
 
 
 
     (*(out+1)).z -= (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+1)).w -= (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;
     (*(out+1)).z -= (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+1)).w -= (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;     
     
     (*(out+4)).z -= (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+4)).w -= (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;    
     (*(out+4)).z -= (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+4)).w -= (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;      
 

     
     (*(out+2)).x -= (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+2)).y -= (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im;
     (*(out+2)).x -= (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+2)).y -= (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;
                
     (*(out+5)).x -= (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+5)).y -= (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;   
     (*(out+5)).x -= (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+5)).y -= (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im; 
   
   
        
     (*(out+2)).z -= (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+2)).w -= (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im;
     (*(out+2)).z -= (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+2)).w -= (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;
               
     (*(out+5)).z -= (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+5)).w -= (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;   
     (*(out+5)).z -= (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+5)).w -= (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im; 
  
}






//-kappa(r + gamma_mu)
__device__ void dev_kappaP0_minus(dev_spinor * out, dev_spinor * in, dev_complex kappa){


     (*(out+0)).x -= (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+0)).y -= (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im;
     (*(out+0)).x += (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+0)).y += (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;
     
     (*(out+3)).x -= (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+3)).y -= (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;   
     (*(out+3)).x += (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+3)).y += (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im; 



     (*(out+0)).z -= (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+0)).w -= (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;
     (*(out+0)).z += (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+0)).w += (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;     
          
     (*(out+3)).z -= (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+3)).w -= (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;
     (*(out+3)).z += (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+3)).w += (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;

 
 
     (*(out+1)).x -= (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+1)).y -= (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im;
     (*(out+1)).x += (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+1)).y += (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;     
     
     (*(out+4)).x -= (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+4)).y -= (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;    
     (*(out+4)).x += (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+4)).y += (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im; 
 
 
 
     (*(out+1)).z -= (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+1)).w -= (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;
     (*(out+1)).z += (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+1)).w += (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;     
     
     (*(out+4)).z -= (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+4)).w -= (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;    
     (*(out+4)).z += (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+4)).w += (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;      
 

     
     (*(out+2)).x -= (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+2)).y -= (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im;
     (*(out+2)).x += (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+2)).y += (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;
                
     (*(out+5)).x -= (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+5)).y -= (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;   
     (*(out+5)).x += (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+5)).y += (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im; 
   
   
       
     (*(out+2)).z -= (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+2)).w -= (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im;
     (*(out+2)).z += (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+2)).w += (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;
               
     (*(out+5)).z -= (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+5)).w -= (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;   
     (*(out+5)).z += (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+5)).w += (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im; 
  
}







#ifdef RELATIVISTIC_BASIS
//  here comes P0+- for the relativistic basis
//  in this basis we have:
//
//  gamma0 =
//  -1  0  0  0 
//   0 -1  0  0
//   0  0  1  0
//   0  0  0  1
//
//



//-kappa(r - gamma_mu)
__device__ void dev_kappaP0_plus_relativistic(dev_spinor * out, dev_spinor * in, dev_complex kappa){


     (*(out+0)).x -= 2.0*( (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im );
     (*(out+0)).y -= 2.0*( (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im );
     (*(out+0)).z -= 2.0*( (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im );
     (*(out+0)).w -= 2.0*( (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im );
     
     (*(out+1)).x -= 2.0*( (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im );
     (*(out+1)).y -= 2.0*( (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im );     
     (*(out+1)).z -= 2.0*( (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im );
     (*(out+1)).w -= 2.0*( (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im );     
     
     (*(out+2)).x -= 2.0*( (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im );
     (*(out+2)).y -= 2.0*( (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im );
     (*(out+2)).z -= 2.0*( (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im );
     (*(out+2)).w -= 2.0*( (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im );

}






//-kappa(r + gamma_mu)
__device__ void dev_kappaP0_minus_relativistic(dev_spinor * out, dev_spinor * in, dev_complex kappa){

    
     (*(out+3)).x -= 2.0*( (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im );
     (*(out+3)).y -= 2.0*( (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im );     
     (*(out+3)).z -= 2.0*( (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im );
     (*(out+3)).w -= 2.0*( (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im );


     (*(out+4)).x -= 2.0*( (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im );
     (*(out+4)).y -= 2.0*( (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im );     
     (*(out+4)).z -= 2.0*( (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im );
     (*(out+4)).w -= 2.0*( (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im );    

           
     (*(out+5)).x -= 2.0*( (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im );
     (*(out+5)).y -= 2.0*( (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im );            
     (*(out+5)).z -= 2.0*( (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im );
     (*(out+5)).w -= 2.0*( (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im );   

  
}

//RELATIVISTIC_BASIS
#endif 





//applies the Hopping Part Even-Odd !
//the gauge field is the complete gaugefield!
//the gauge field at the local point is reconstructed by 2*pos+eo where pos is the eo-position
//from 0..VOLUME/2-1, eo = 0 or 1
//the positions in the gauge fields are passed in "gfindex_site" for gf's that are attached at
//the actual positions and in "gfindex_nextsite" for gf's that start at a position of the 
//other eo-sublattice.
//for the hopping positions of the eo-spinor field we use on of the two dedicated eo-nn fields
//the boundary conditions are implemented as in Hopping_Matrix.c
//mult with complex conjugate k0,k1,k2,k3 in positive direction because
// psi(x+mu) != exp(i theta_mu) psi(x)  
// we start from site index start and go to start+size
// by this we can split up bulk and rand in mpi


__device__ void dev_hopping_kernel(const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * ssum, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo, int pos){

  int hoppos;
    dev_spinor shelp1[6];


  #ifdef GPU_3DBLOCK
    dev_su3_pad gfsmem;  
  #else
    dev_su3_pad gfsmem;  
  #endif


 
  //dev_zero_spinor_local(&(ssum[0])); // zero sum 
  ssum[0].x=0.0f; ssum[0].y=0.0f; ssum[0].z=0.0f; ssum[0].w=0.0f;
  ssum[1].x=0.0f; ssum[1].y=0.0f; ssum[1].z=0.0f; ssum[1].w=0.0f;
  ssum[2].x=0.0f; ssum[2].y=0.0f; ssum[2].z=0.0f; ssum[2].w=0.0f;
  ssum[3].x=0.0f; ssum[3].y=0.0f; ssum[3].z=0.0f; ssum[3].w=0.0f;  
  ssum[4].x=0.0f; ssum[4].y=0.0f; ssum[4].z=0.0f; ssum[4].w=0.0f;
  ssum[5].x=0.0f; ssum[5].y=0.0f; ssum[5].z=0.0f; ssum[5].w=0.0f;    
  
 #ifdef TEMPORALGAUGE
  int spatialvol = dev_LX*dev_LY*dev_LZ;
 #endif
  

 //the volume for the gaugefield is 2x volume of spinors!!!
 #ifdef MPI
   int gaugevol = 2*dev_VOLUMEPLUSRAND;
 #else
   int gaugevol = 2*dev_VOLUME; 
 #endif 

//hopping term                
//l==0,t
            //positive direction
            hoppos = nn_evenodd[8*pos];
             //hoppos = tex1Dfetch(nn_tex,8*pos);
            //color
            
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_site[pos]) < (dev_T-1)*spatialvol) || (dev_rank < dev_nproc-1) ) {
                //if ((gfindex_site[pos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_site[pos]/spatialvol) != (dev_T-1) ) {
              #endif
              
              #ifdef USETEXTURE 
                shelp1[0] = tex1Dfetch(spin_tex0,hoppos);
                shelp1[1] = tex1Dfetch(spin_tex1,hoppos);
                shelp1[2] = tex1Dfetch(spin_tex2,hoppos);
                #ifdef RELATIVISTIC_BASIS
                  shelp1[3].x = 0.0f; shelp1[3].y = 0.0f; shelp1[3].z = 0.0f; shelp1[3].w = 0.0f;
		  shelp1[4].x = 0.0f; shelp1[4].y = 0.0f; shelp1[4].z = 0.0f; shelp1[4].w = 0.0f;
		  shelp1[5].x = 0.0f; shelp1[5].y = 0.0f; shelp1[5].z = 0.0f; shelp1[5].w = 0.0f;
		#else
		  shelp1[3] = tex1Dfetch(spin_tex3,hoppos);
                  shelp1[4] = tex1Dfetch(spin_tex4,hoppos);
                  shelp1[5] = tex1Dfetch(spin_tex5,hoppos);
		#endif
              #else
                shelp1[0] = sin[hoppos+0*DEVOFF];
                shelp1[1] = sin[hoppos+1*DEVOFF];
                shelp1[2] = sin[hoppos+2*DEVOFF];
                #ifdef RELATIVISTIC_BASIS
                  shelp1[3].x = 0.0f; shelp1[3].y = 0.0f; shelp1[3].z = 0.0f; shelp1[3].w = 0.0f;
		  shelp1[4].x = 0.0f; shelp1[4].y = 0.0f; shelp1[4].z = 0.0f; shelp1[4].w = 0.0f;
		  shelp1[5].x = 0.0f; shelp1[5].y = 0.0f; shelp1[5].z = 0.0f; shelp1[5].w = 0.0f;
                #else
		  shelp1[3] = sin[hoppos+3*DEVOFF];
                  shelp1[4] = sin[hoppos+4*DEVOFF];
                  shelp1[5] = sin[hoppos+5*DEVOFF];
                #endif
	      #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem.m));
                #else
                dev_reconstructgf_2vtexref(gf,gfindex_site[pos], 0, gaugevol ,&(gfsmem.m));
                #endif
                
                #ifdef RELATIVISTIC_BASIS
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex_rel_up(gfsmem.m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV_rel_up(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #else
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #endif
              }
            #else
              #ifdef GF_8
              dev_reconstructgf_8texref(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem.m));
              #else
              dev_reconstructgf_2vtexref(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem.m));
              #endif
              #ifdef RELATIVISTIC_BASIS
                #ifdef USETEXTURE
                  dev_su3MtV_spintex_rel_up(gfsmem.m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV_rel_up(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                #endif
              #else
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                #endif
              #endif
            #endif
            
            #ifdef RELATIVISTIC_BASIS
              dev_kappaP0_plus_relativistic(&(ssum[0]), &(shelp1[0]), dev_cconj(dev_k0));
            #else
            //-kappa(r - gamma_mu)
              dev_kappaP0_plus(&(ssum[0]), &(shelp1[0]), dev_cconj(dev_k0));
            #endif
	    
//l==0,t
            //negative direction
            hoppos = nn_evenodd[8*pos+4]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+4);
            //color
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_nextsite[hoppos]) < (dev_T-1)*spatialvol) || (dev_rank > 0) ) {
                //if ((gfindex_nextsite[hoppos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_nextsite[hoppos]/spatialvol) != (dev_T-1) ) {
              #endif
              
               #ifdef USETEXTURE
                #ifdef RELATIVISTIC_BASIS
                  shelp1[0].x = 0.0f; shelp1[0].y = 0.0f; shelp1[0].z = 0.0f; shelp1[0].w = 0.0f;
		  shelp1[1].x = 0.0f; shelp1[1].y = 0.0f; shelp1[1].z = 0.0f; shelp1[1].w = 0.0f;
		  shelp1[2].x = 0.0f; shelp1[2].y = 0.0f; shelp1[2].z = 0.0f; shelp1[2].w = 0.0f;
                #else
                  shelp1[0] = tex1Dfetch(spin_tex0,hoppos);
                  shelp1[1] = tex1Dfetch(spin_tex1,hoppos);
                  shelp1[2] = tex1Dfetch(spin_tex2,hoppos);
                #endif
		shelp1[3] = tex1Dfetch(spin_tex3,hoppos);
                shelp1[4] = tex1Dfetch(spin_tex4,hoppos);
                shelp1[5] = tex1Dfetch(spin_tex5,hoppos);
               #else
                #ifdef RELATIVISTIC_BASIS
                  shelp1[0].x = 0.0f; shelp1[0].y = 0.0f; shelp1[0].z = 0.0f; shelp1[0].w = 0.0f;
		  shelp1[1].x = 0.0f; shelp1[1].y = 0.0f; shelp1[1].z = 0.0f; shelp1[1].w = 0.0f;
		  shelp1[2].x = 0.0f; shelp1[2].y = 0.0f; shelp1[2].z = 0.0f; shelp1[2].w = 0.0f;
                #else
                  shelp1[0] = sin[hoppos+0*DEVOFF];
                  shelp1[1] = sin[hoppos+1*DEVOFF];
                  shelp1[2] = sin[hoppos+2*DEVOFF];
                #endif
		shelp1[3] = sin[hoppos+3*DEVOFF];
                shelp1[4] = sin[hoppos+4*DEVOFF];
                shelp1[5] = sin[hoppos+5*DEVOFF];
               #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem.m));
                #else
                dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem.m));
                #endif
                #ifdef RELATIVISTIC_BASIS
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex_rel_down(gfsmem.m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV_rel_down(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #else
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #endif
              }
            #else            
              #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem.m));
              #else
              dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem.m));
              #endif
              #ifdef RELATIVISTIC_BASIS
                #ifdef USETEXTURE
                  dev_su3MtV_spintex_rel_down(gfsmem.m, hoppos, &(shelp1[0]));  
                #else
                  dev_su3MtV_rel_down(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                #endif 
              #else
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));  
                #else
                  dev_su3MtV(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                #endif 
              #endif
            #endif
            
            #ifdef RELATIVISTIC_BASIS
              dev_kappaP0_minus_relativistic(&(ssum[0]), &(shelp1[0]), dev_k0);
            #else
              //-kappa(r + gamma_mu)
              dev_kappaP0_minus(&(ssum[0]), &(shelp1[0]), dev_k0);
            #endif




//l==3,z 
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,gfindex_site[pos], 3, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref(gf, gfindex_site[pos], 3, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet  
              //dev_kappaP3_plus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
	      
              dev_su3MtV_kappaP3_plus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k3);	      
	    #else
              dev_su3MtV_kappaP3_plus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k3);
            #endif
            

//l==3,z               
            
            //negative direction
            hoppos = nn_evenodd[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 3, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 3, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet 
              //dev_kappaP3_minus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
	      
	      dev_su3MtV_kappaP3_minus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k3);	      
	    #else
	      dev_su3MtV_kappaP3_minus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k3);
            #endif
            





//l==2,y 
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,gfindex_site[pos], 2, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref(gf,gfindex_site[pos], 2, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet
              //dev_kappaP2_plus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
	      
	      dev_su3MtV_kappaP2_plus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k2);	      
	    #else
	      dev_su3MtV_kappaP2_plus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k2);
            #endif
            


//l==2,y        

            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 2, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 2, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet
              //dev_kappaP2_minus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
	      
	      dev_su3MtV_kappaP2_minus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k2);	      
	    #else
	      dev_su3MtV_kappaP2_minus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k2);
            #endif
            
            



//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,gfindex_site[pos], 1, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref(gf,gfindex_site[pos], 1, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet
              //dev_kappaP1_plus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
	      
	      dev_su3MtV_kappaP1_plus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k1);	      
	    #else
	      dev_su3MtV_kappaP1_plus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k1);
            #endif
            
            


//l==1,x 
            
            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 1, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 1, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet
              //dev_kappaP1_minus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
	      
	      dev_su3MtV_kappaP1_minus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k1);	      
	    #else
	      dev_su3MtV_kappaP1_minus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k1);
            #endif
            
           
}














__global__ void dev_Hopping_Matrix(const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo, int start, int size){

  int pos;
    dev_spinor ssum[6];


  #ifdef GPU_3DBLOCK
    pos = start  
          + (threadIdx.z + blockDim.z*(threadIdx.y + blockDim.y*(threadIdx.x))) 
          + blockDim.z*blockDim.y*blockDim.x*blockIdx.x;     
  #else
    pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;  
  #endif

  
  
  if(pos < (start + size)){
 
     //apply hopping kernel
     dev_hopping_kernel(gf, sin, &(ssum[0]), gfindex_site, gfindex_nextsite, nn_evenodd, eo, pos);

       
        //copy to output spinor
        //dev_write_spinor(&(ssum[0]),&(sout[pos])); 
	sout[pos+0*DEVOFF].x = ssum[0].x; sout[pos+0*DEVOFF].y = ssum[0].y; 
	sout[pos+0*DEVOFF].z = ssum[0].z; sout[pos+0*DEVOFF].w = ssum[0].w;	
	sout[pos+1*DEVOFF].x = ssum[1].x; sout[pos+1*DEVOFF].y = ssum[1].y; 
	sout[pos+1*DEVOFF].z = ssum[1].z; sout[pos+1*DEVOFF].w = ssum[1].w;
	sout[pos+2*DEVOFF].x = ssum[2].x; sout[pos+2*DEVOFF].y = ssum[2].y; 
	sout[pos+2*DEVOFF].z = ssum[2].z; sout[pos+2*DEVOFF].w = ssum[2].w;
	sout[pos+3*DEVOFF].x = ssum[3].x; sout[pos+3*DEVOFF].y = ssum[3].y; 
	sout[pos+3*DEVOFF].z = ssum[3].z; sout[pos+3*DEVOFF].w = ssum[3].w;
	sout[pos+4*DEVOFF].x = ssum[4].x; sout[pos+4*DEVOFF].y = ssum[4].y; 
	sout[pos+4*DEVOFF].z = ssum[4].z; sout[pos+4*DEVOFF].w = ssum[4].w;	
	sout[pos+5*DEVOFF].x = ssum[5].x; sout[pos+5*DEVOFF].y = ssum[5].y; 
	sout[pos+5*DEVOFF].z = ssum[5].z; sout[pos+5*DEVOFF].w = ssum[5].w;	
	
  }
}





//applies Hopping Matrix and multiplies output spinor with (1 +- i mubar gamma5)
__global__ void dev_Hopping_Matrix_ext(const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout,dev_spinor * sout2,  float sign, 
				       const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, 
				       const int eo, int start, int size){

  int pos;
    dev_spinor ssum[6], shelp1[6];


  #ifdef GPU_3DBLOCK
    pos = start  
          + (threadIdx.z + blockDim.z*(threadIdx.y + blockDim.y*(threadIdx.x))) 
          + blockDim.z*blockDim.y*blockDim.x*blockIdx.x;     
  #else 
    pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;  
  #endif

  
  
  if(pos < (start + size)){
 
 
     //apply hopping kernel
     dev_hopping_kernel(gf, sin, &(ssum[0]), gfindex_site, gfindex_nextsite, nn_evenodd, eo, pos);

         
	sout[pos+0*DEVOFF].x = ssum[0].x; sout[pos+0*DEVOFF].y = ssum[0].y; 
	sout[pos+0*DEVOFF].z = ssum[0].z; sout[pos+0*DEVOFF].w = ssum[0].w;	
	sout[pos+1*DEVOFF].x = ssum[1].x; sout[pos+1*DEVOFF].y = ssum[1].y; 
	sout[pos+1*DEVOFF].z = ssum[1].z; sout[pos+1*DEVOFF].w = ssum[1].w;
	sout[pos+2*DEVOFF].x = ssum[2].x; sout[pos+2*DEVOFF].y = ssum[2].y; 
	sout[pos+2*DEVOFF].z = ssum[2].z; sout[pos+2*DEVOFF].w = ssum[2].w;
	sout[pos+3*DEVOFF].x = ssum[3].x; sout[pos+3*DEVOFF].y = ssum[3].y; 
	sout[pos+3*DEVOFF].z = ssum[3].z; sout[pos+3*DEVOFF].w = ssum[3].w;
	sout[pos+4*DEVOFF].x = ssum[4].x; sout[pos+4*DEVOFF].y = ssum[4].y; 
	sout[pos+4*DEVOFF].z = ssum[4].z; sout[pos+4*DEVOFF].w = ssum[4].w;	
	sout[pos+5*DEVOFF].x = ssum[5].x; sout[pos+5*DEVOFF].y = ssum[5].y; 
	sout[pos+5*DEVOFF].z = ssum[5].z; sout[pos+5*DEVOFF].w = ssum[5].w;	


	dev_complex pm_imu = dev_initcomplex(0.0, sign * mubar);
	#ifdef RELATIVISTIC_BASIS
	  dev_skalarmult_gamma5_spinor_rel(&(shelp1[0]), pm_imu, &(ssum[0]) );
	#else
	   dev_skalarmult_gamma5_spinor(&(shelp1[0]), pm_imu, &(ssum[0]) );
	#endif            
        dev_add_spinor_assign(&(shelp1[0]), &(ssum[0]));					// slocal  =  slocal + sin  =  pm_imu * (gamma5) * sin + sin

        //copy to output spinor
        //dev_write_spinor(&(ssum[0]),&(sout[pos])); 
	sout2[pos+0*DEVOFF].x = shelp1[0].x; sout2[pos+0*DEVOFF].y = shelp1[0].y; 
	sout2[pos+0*DEVOFF].z = shelp1[0].z; sout2[pos+0*DEVOFF].w = shelp1[0].w;	
	sout2[pos+1*DEVOFF].x = shelp1[1].x; sout2[pos+1*DEVOFF].y = shelp1[1].y; 
	sout2[pos+1*DEVOFF].z = shelp1[1].z; sout2[pos+1*DEVOFF].w = shelp1[1].w;
	sout2[pos+2*DEVOFF].x = shelp1[2].x; sout2[pos+2*DEVOFF].y = shelp1[2].y; 
	sout2[pos+2*DEVOFF].z = shelp1[2].z; sout2[pos+2*DEVOFF].w = shelp1[2].w;
	sout2[pos+3*DEVOFF].x = shelp1[3].x; sout2[pos+3*DEVOFF].y = shelp1[3].y; 
	sout2[pos+3*DEVOFF].z = shelp1[3].z; sout2[pos+3*DEVOFF].w = shelp1[3].w;
	sout2[pos+4*DEVOFF].x = shelp1[4].x; sout2[pos+4*DEVOFF].y = shelp1[4].y; 
	sout2[pos+4*DEVOFF].z = shelp1[4].z; sout2[pos+4*DEVOFF].w = shelp1[4].w;	
	sout2[pos+5*DEVOFF].x = shelp1[5].x; sout2[pos+5*DEVOFF].y = shelp1[5].y; 
	sout2[pos+5*DEVOFF].z = shelp1[5].z; sout2[pos+5*DEVOFF].w = shelp1[5].w;	
	

  }
}






// sout = gamma5 (1+sign imubar gamma5) sin2 - H sin
__global__ void dev_Hopping_Matrix_ext2(const dev_su3_2v * gf, const dev_spinor * sin, 
                                        dev_spinor * sout, float sign,dev_spinor * sin2,   
				        const int * gfindex_site, const int* gfindex_nextsite, 
                                        const int * nn_evenodd, const int eo, int start, int size){

  int pos;
    dev_spinor shelp1[6], ssum[6];


  #ifdef GPU_3DBLOCK
    pos = start  
          + (threadIdx.z + blockDim.z*(threadIdx.y + blockDim.y*(threadIdx.x))) 
          + blockDim.z*blockDim.y*blockDim.x*blockIdx.x;     
  #else
    pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;  
  #endif

  
  
  if(pos < (start + size)){
 
 
//apply hopping kernel
    dev_hopping_kernel(gf, sin, &(ssum[0]), gfindex_site, gfindex_nextsite, nn_evenodd, eo, pos);

         
     dev_complex pm_imu = dev_initcomplex(0.0, sign*twokappamu); // i mutilde
     dev_spinor slocal[6];
     #ifdef RELATIVISTIC_BASIS
       dev_skalarmult_gamma5_globalspinor_rel(&(slocal[0]), pm_imu, &(sin2[pos]) );
     #else
       dev_skalarmult_gamma5_globalspinor(&(slocal[0]),pm_imu,&(sin2[pos]));
     #endif
     dev_add_globalspinor_assign(&(slocal[0]), &(sin2[pos]));
     dev_sub_spinor_assign(&(slocal[0]), &(ssum[0]));
     #ifdef RELATIVISTIC_BASIS
       dev_Gamma5_assigntoglobal_rel(&(sout[pos]), &(slocal[0]));
     #else
       dev_Gamma5_assigntoglobal(&(sout[pos]), &(slocal[0]));
     #endif


  }
}



//applies Hopping Matrix and multiplies output spinor with 1.0/(1.0 + twokappamu*twokappamu)*(1 +- i mubar gamma5)
__global__ void dev_Hopping_Matrix_ext3(const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout, float sign, 
				       const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, 
				       const int eo, int start, int size){

  int pos;
    dev_spinor ssum[6], shelp1[6];


  #ifdef GPU_3DBLOCK
    pos = start  
          + (threadIdx.z + blockDim.z*(threadIdx.y + blockDim.y*(threadIdx.x))) 
          + blockDim.z*blockDim.y*blockDim.x*blockIdx.x;     
  #else 
    pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;  
  #endif

  
  
  if(pos < (start + size)){

     //apply hopping kernel
     dev_hopping_kernel(gf, sin, &(ssum[0]), gfindex_site, gfindex_nextsite, nn_evenodd, eo, pos);

     // 1/(1+mubar^2) (1 +- imubar gamma5) H*in
        float one_plus_musquare_inv = 1.0/(1.0 + twokappamu*twokappamu);
	dev_complex pm_imu = dev_initcomplex(0.0, -1.0*sign * twokappamu);
	#ifdef RELATIVISTIC_BASIS
	  dev_skalarmult_gamma5_spinor_rel(&(shelp1[0]), pm_imu, &(ssum[0]) );
	#else
	   dev_skalarmult_gamma5_spinor(&(shelp1[0]), pm_imu, &(ssum[0]) );
	#endif  
	dev_add_spinor_assign(&(shelp1[0]), &(ssum[0]));

	dev_realmult_spinor_assign(&(ssum[0]), one_plus_musquare_inv, &(shelp1[0]));					// slocal  =  slocal + sin  =  pm_imu * (gamma5) * sin + sin

	sout[pos+0*DEVOFF].x = ssum[0].x; sout[pos+0*DEVOFF].y = ssum[0].y; 
	sout[pos+0*DEVOFF].z = ssum[0].z; sout[pos+0*DEVOFF].w = ssum[0].w;	
	sout[pos+1*DEVOFF].x = ssum[1].x; sout[pos+1*DEVOFF].y = ssum[1].y; 
	sout[pos+1*DEVOFF].z = ssum[1].z; sout[pos+1*DEVOFF].w = ssum[1].w;
	sout[pos+2*DEVOFF].x = ssum[2].x; sout[pos+2*DEVOFF].y = ssum[2].y; 
	sout[pos+2*DEVOFF].z = ssum[2].z; sout[pos+2*DEVOFF].w = ssum[2].w;
	sout[pos+3*DEVOFF].x = ssum[3].x; sout[pos+3*DEVOFF].y = ssum[3].y; 
	sout[pos+3*DEVOFF].z = ssum[3].z; sout[pos+3*DEVOFF].w = ssum[3].w;
	sout[pos+4*DEVOFF].x = ssum[4].x; sout[pos+4*DEVOFF].y = ssum[4].y; 
	sout[pos+4*DEVOFF].z = ssum[4].z; sout[pos+4*DEVOFF].w = ssum[4].w;	
	sout[pos+5*DEVOFF].x = ssum[5].x; sout[pos+5*DEVOFF].y = ssum[5].y; 
	sout[pos+5*DEVOFF].z = ssum[5].z; sout[pos+5*DEVOFF].w = ssum[5].w;	

     //dev_realmult_spinor_assigntoglobal(&(sout[pos]), one_plus_musquare_inv, &(shelp1[0]));
  }
}












__global__ void dev_Hopping_Matrix_updn(const dev_su3_2v * gf, const dev_spinor * sin_up, const dev_spinor * sin_dn, dev_spinor * sout_up, dev_spinor * sout_dn , const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo, int start, int size){

  int pos,hoppos;
    dev_spinor shelp1_up[6], ssum_up[6];
    dev_spinor shelp1_dn[6], ssum_dn[6];

  #ifdef GPU_3DBLOCK
    dev_su3 gfsmem;  
    pos = start  
          + (threadIdx.z + blockDim.z*(threadIdx.y + blockDim.y*(threadIdx.x))) 
          + blockDim.z*blockDim.y*blockDim.x*blockIdx.x;     
  #else
    dev_su3 gfsmem;  
    pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;  
  #endif

  
  
  if(pos < (start + size)){
 
  //dev_zero_spinor_local(&(ssum[0])); // zero sum 
  ssum_up[0].x=0.0f; ssum_up[0].y=0.0f; ssum_up[0].z=0.0f; ssum_up[0].w=0.0f;
  ssum_up[1].x=0.0f; ssum_up[1].y=0.0f; ssum_up[1].z=0.0f; ssum_up[1].w=0.0f;
  ssum_up[2].x=0.0f; ssum_up[2].y=0.0f; ssum_up[2].z=0.0f; ssum_up[2].w=0.0f;
  ssum_up[3].x=0.0f; ssum_up[3].y=0.0f; ssum_up[3].z=0.0f; ssum_up[3].w=0.0f;  
  ssum_up[4].x=0.0f; ssum_up[4].y=0.0f; ssum_up[4].z=0.0f; ssum_up[4].w=0.0f;
  ssum_up[5].x=0.0f; ssum_up[5].y=0.0f; ssum_up[5].z=0.0f; ssum_up[5].w=0.0f;    

  ssum_dn[0].x=0.0f; ssum_dn[0].y=0.0f; ssum_dn[0].z=0.0f; ssum_dn[0].w=0.0f;
  ssum_dn[1].x=0.0f; ssum_dn[1].y=0.0f; ssum_dn[1].z=0.0f; ssum_dn[1].w=0.0f;
  ssum_dn[2].x=0.0f; ssum_dn[2].y=0.0f; ssum_dn[2].z=0.0f; ssum_dn[2].w=0.0f;
  ssum_dn[3].x=0.0f; ssum_dn[3].y=0.0f; ssum_dn[3].z=0.0f; ssum_dn[3].w=0.0f;  
  ssum_dn[4].x=0.0f; ssum_dn[4].y=0.0f; ssum_dn[4].z=0.0f; ssum_dn[4].w=0.0f;
  ssum_dn[5].x=0.0f; ssum_dn[5].y=0.0f; ssum_dn[5].z=0.0f; ssum_dn[5].w=0.0f;    
  
 #ifdef TEMPORALGAUGE
  int spatialvol = dev_LX*dev_LY*dev_LZ;
 #endif
  

 //the volume for the gaugefield is 2x volume of spinors!!!
 #ifdef MPI
   int gaugevol = 2*dev_VOLUMEPLUSRAND;
 #else
   int gaugevol = 2*dev_VOLUME; 
 #endif 

//hopping term                
//l==0,t
            //positive direction
            hoppos = nn_evenodd[8*pos];
             //hoppos = tex1Dfetch(nn_tex,8*pos);
            //color
            
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_site[pos]) < (dev_T-1)*spatialvol) || (dev_rank < dev_nproc-1) ) {
                //if ((gfindex_site[pos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_site[pos]/spatialvol) != (dev_T-1) ) {
              #endif
              
              #ifdef USETEXTURE 
                shelp1_up[0] = tex1Dfetch(spin_tex0,hoppos);
                shelp1_up[1] = tex1Dfetch(spin_tex1,hoppos);
                shelp1_up[2] = tex1Dfetch(spin_tex2,hoppos);
                shelp1_dn[0] = tex1Dfetch(spin_tex_dn0,hoppos);
                shelp1_dn[1] = tex1Dfetch(spin_tex_dn1,hoppos);
                shelp1_dn[2] = tex1Dfetch(spin_tex_dn2,hoppos);		
                #ifdef RELATIVISTIC_BASIS
                  shelp1_up[3].x = 0.0f; shelp1_up[3].y = 0.0f; shelp1_up[3].z = 0.0f; shelp1_up[3].w = 0.0f;
		  shelp1_up[4].x = 0.0f; shelp1_up[4].y = 0.0f; shelp1_up[4].z = 0.0f; shelp1_up[4].w = 0.0f;
		  shelp1_up[5].x = 0.0f; shelp1_up[5].y = 0.0f; shelp1_up[5].z = 0.0f; shelp1_up[5].w = 0.0f;
                  shelp1_dn[3].x = 0.0f; shelp1_dn[3].y = 0.0f; shelp1_dn[3].z = 0.0f; shelp1_dn[3].w = 0.0f;
		  shelp1_dn[4].x = 0.0f; shelp1_dn[4].y = 0.0f; shelp1_dn[4].z = 0.0f; shelp1_dn[4].w = 0.0f;
		  shelp1_dn[5].x = 0.0f; shelp1_dn[5].y = 0.0f; shelp1_dn[5].z = 0.0f; shelp1_dn[5].w = 0.0f;		  
		#else
		  shelp1_up[3] = tex1Dfetch(spin_tex3,hoppos);
                  shelp1_up[4] = tex1Dfetch(spin_tex4,hoppos);
                  shelp1_up[5] = tex1Dfetch(spin_tex5,hoppos);
		  shelp1_dn[3] = tex1Dfetch(spin_tex_dn3,hoppos);
                  shelp1_dn[4] = tex1Dfetch(spin_tex_dn4,hoppos);
                  shelp1_dn[5] = tex1Dfetch(spin_tex_dn5,hoppos);		  
		#endif
              #else
                shelp1_up[0] = sin_up[hoppos+0*DEVOFF];
                shelp1_up[1] = sin_up[hoppos+1*DEVOFF];
                shelp1_up[2] = sin_up[hoppos+2*DEVOFF];
                shelp1_dn[0] = sin_dn[hoppos+0*DEVOFF];
                shelp1_dn[1] = sin_dn[hoppos+1*DEVOFF];
                shelp1_dn[2] = sin_dn[hoppos+2*DEVOFF];		
                #ifdef RELATIVISTIC_BASIS
                  shelp1_up[3].x = 0.0f; shelp1_up[3].y = 0.0f; shelp1_up[3].z = 0.0f; shelp1_up[3].w = 0.0f;
		  shelp1_up[4].x = 0.0f; shelp1_up[4].y = 0.0f; shelp1_up[4].z = 0.0f; shelp1_up[4].w = 0.0f;
		  shelp1_up[5].x = 0.0f; shelp1_up[5].y = 0.0f; shelp1_up[5].z = 0.0f; shelp1_up[5].w = 0.0f;
                  shelp1_dn[3].x = 0.0f; shelp1_dn[3].y = 0.0f; shelp1_dn[3].z = 0.0f; shelp1_dn[3].w = 0.0f;
		  shelp1_dn[4].x = 0.0f; shelp1_dn[4].y = 0.0f; shelp1_dn[4].z = 0.0f; shelp1_dn[4].w = 0.0f;
		  shelp1_dn[5].x = 0.0f; shelp1_dn[5].y = 0.0f; shelp1_dn[5].z = 0.0f; shelp1_dn[5].w = 0.0f;
		  #else
		  shelp1_up[3] = sin_up[hoppos+3*DEVOFF];
                  shelp1_up[4] = sin_up[hoppos+4*DEVOFF];
                  shelp1_up[5] = sin_up[hoppos+5*DEVOFF];
		  shelp1_dn[3] = sin_dn[hoppos+3*DEVOFF];
                  shelp1_dn[4] = sin_dn[hoppos+4*DEVOFF];
                  shelp1_dn[5] = sin_dn[hoppos+5*DEVOFF];		  
                #endif
	      #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem));
                #else
                dev_reconstructgf_2vtexref(gf,gfindex_site[pos], 0, gaugevol ,&(gfsmem));
                #endif
                
                #ifdef RELATIVISTIC_BASIS
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex_rel_up(gfsmem, hoppos, &(shelp1_up[0]));
                    dev_su3MtV_spintex2_rel_up(gfsmem, hoppos, &(shelp1_dn[0]));		    
                  #else
                    dev_su3MtV_rel_up(gfsmem, &(sin_up[hoppos]), &(shelp1_up[0]));
                    dev_su3MtV_rel_up(gfsmem, &(sin_dn[hoppos]), &(shelp1_dn[0]));		    
                  #endif
                #else
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex(gfsmem, hoppos, &(shelp1_up[0]));
                    dev_su3MtV_spintex2(gfsmem, hoppos, &(shelp1_dn[0]));		    
                  #else
                    dev_su3MtV(gfsmem, &(sin_up[hoppos]), &(shelp1_up[0]));
                    dev_su3MtV(gfsmem, &(sin_dn[hoppos]), &(shelp1_dn[0]));		    
                  #endif
                #endif
              }
            #else
              #ifdef GF_8
              dev_reconstructgf_8texref(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem));
              #else
              dev_reconstructgf_2vtexref(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem));
              #endif
              #ifdef RELATIVISTIC_BASIS
                #ifdef USETEXTURE
                  dev_su3MtV_spintex_rel_up(gfsmem, hoppos, &(shelp1_up[0]));
                  dev_su3MtV_spintex2_rel_up(gfsmem, hoppos, &(shelp1_dn[0]));		  
                #else
                  dev_su3MtV_rel_up(gfsmem, &(sin_up[hoppos]), &(shelp1_up[0]));
                  dev_su3MtV_rel_up(gfsmem, &(sin_dn[hoppos]), &(shelp1_dn[0]));		  
                #endif
              #else
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem, hoppos, &(shelp1_up[0]));
                  dev_su3MtV_spintex2(gfsmem, hoppos, &(shelp1_dn[0]));		  
                #else
                  dev_su3MtV(gfsmem, &(sin_up[hoppos]), &(shelp1_up[0]));
                  dev_su3MtV(gfsmem, &(sin_dn[hoppos]), &(shelp1_dn[0]));		  
                #endif
              #endif
            #endif
            
            #ifdef RELATIVISTIC_BASIS
              dev_kappaP0_plus_relativistic(&(ssum_up[0]), &(shelp1_up[0]), dev_cconj(dev_k0));
              dev_kappaP0_plus_relativistic(&(ssum_dn[0]), &(shelp1_dn[0]), dev_cconj(dev_k0));	      
            #else
            //-kappa(r - gamma_mu)
              dev_kappaP0_plus(&(ssum_up[0]), &(shelp1_up[0]), dev_cconj(dev_k0));
              dev_kappaP0_plus(&(ssum_dn[0]), &(shelp1_dn[0]), dev_cconj(dev_k0));	      
            #endif
	    
//l==0,t
            //negative direction
            hoppos = nn_evenodd[8*pos+4]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+4);
            //color
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_nextsite[hoppos]) < (dev_T-1)*spatialvol) || (dev_rank > 0) ) {
                //if ((gfindex_nextsite[hoppos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_nextsite[hoppos]/spatialvol) != (dev_T-1) ) {
              #endif
              
               #ifdef USETEXTURE
                #ifdef RELATIVISTIC_BASIS
                  shelp1_up[0].x = 0.0f; shelp1_up[0].y = 0.0f; shelp1_up[0].z = 0.0f; shelp1_up[0].w = 0.0f;
		  shelp1_up[1].x = 0.0f; shelp1_up[1].y = 0.0f; shelp1_up[1].z = 0.0f; shelp1_up[1].w = 0.0f;
		  shelp1_up[2].x = 0.0f; shelp1_up[2].y = 0.0f; shelp1_up[2].z = 0.0f; shelp1_up[2].w = 0.0f;
		  shelp1_dn[0].x = 0.0f; shelp1_dn[0].y = 0.0f; shelp1_dn[0].z = 0.0f; shelp1_dn[0].w = 0.0f;
		  shelp1_dn[1].x = 0.0f; shelp1_dn[1].y = 0.0f; shelp1_dn[1].z = 0.0f; shelp1_dn[1].w = 0.0f;
		  shelp1_dn[2].x = 0.0f; shelp1_dn[2].y = 0.0f; shelp1_dn[2].z = 0.0f; shelp1_dn[2].w = 0.0f;
                #else
                  shelp1_up[0] = tex1Dfetch(spin_tex0,hoppos);
                  shelp1_up[1] = tex1Dfetch(spin_tex1,hoppos);
                  shelp1_up[2] = tex1Dfetch(spin_tex2,hoppos);
                  shelp1_dn[0] = tex1Dfetch(spin_tex_dn0,hoppos);
                  shelp1_dn[1] = tex1Dfetch(spin_tex_dn1,hoppos);
                  shelp1_dn[2] = tex1Dfetch(spin_tex_dn2,hoppos);		  
                #endif
		shelp1_up[3] = tex1Dfetch(spin_tex3,hoppos);
                shelp1_up[4] = tex1Dfetch(spin_tex4,hoppos);
                shelp1_up[5] = tex1Dfetch(spin_tex5,hoppos);
		shelp1_dn[3] = tex1Dfetch(spin_tex_dn3,hoppos);
                shelp1_dn[4] = tex1Dfetch(spin_tex_dn4,hoppos);
                shelp1_dn[5] = tex1Dfetch(spin_tex_dn5,hoppos);		
               #else
                #ifdef RELATIVISTIC_BASIS
                  shelp1_up[0].x = 0.0f; shelp1_up[0].y = 0.0f; shelp1_up[0].z = 0.0f; shelp1_up[0].w = 0.0f;
		  shelp1_up[1].x = 0.0f; shelp1_up[1].y = 0.0f; shelp1_up[1].z = 0.0f; shelp1_up[1].w = 0.0f;
		  shelp1_up[2].x = 0.0f; shelp1_up[2].y = 0.0f; shelp1_up[2].z = 0.0f; shelp1_up[2].w = 0.0f;
                  shelp1_dn[0].x = 0.0f; shelp1_dn[0].y = 0.0f; shelp1_dn[0].z = 0.0f; shelp1_dn[0].w = 0.0f;
		  shelp1_dn[1].x = 0.0f; shelp1_dn[1].y = 0.0f; shelp1_dn[1].z = 0.0f; shelp1_dn[1].w = 0.0f;
		  shelp1_dn[2].x = 0.0f; shelp1_dn[2].y = 0.0f; shelp1_dn[2].z = 0.0f; shelp1_dn[2].w = 0.0f;
		  #else
                  shelp1_up[0] = sin_up[hoppos+0*DEVOFF];
                  shelp1_up[1] = sin_up[hoppos+1*DEVOFF];
                  shelp1_up[2] = sin_up[hoppos+2*DEVOFF];
                  shelp1_dn[0] = sin_dn[hoppos+0*DEVOFF];
                  shelp1_dn[1] = sin_dn[hoppos+1*DEVOFF];
                  shelp1_dn[2] = sin_dn[hoppos+2*DEVOFF];		  
                #endif
		shelp1_up[3] = sin_up[hoppos+3*DEVOFF];
                shelp1_up[4] = sin_up[hoppos+4*DEVOFF];
                shelp1_up[5] = sin_up[hoppos+5*DEVOFF];
		shelp1_dn[3] = sin_dn[hoppos+3*DEVOFF];
                shelp1_dn[4] = sin_dn[hoppos+4*DEVOFF];
                shelp1_dn[5] = sin_dn[hoppos+5*DEVOFF];		
               #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem));
                #else
                dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem));
                #endif
                #ifdef RELATIVISTIC_BASIS
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex_rel_down(gfsmem, hoppos, &(shelp1_up[0]));
                    dev_su3MtV_spintex2_rel_down(gfsmem, hoppos, &(shelp1_dn[0]));		    
                  #else
                    dev_su3MtV_rel_down(gfsmem, &(sin_up[hoppos]), &(shelp1_up[0]));
                    dev_su3MtV_rel_down(gfsmem, &(sin_dn[hoppos]), &(shelp1_dn[0]));		    
                  #endif
                #else
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex(gfsmem, hoppos, &(shelp1_up[0]));
                    dev_su3MtV_spintex2(gfsmem, hoppos, &(shelp1_dn[0]));		    
                  #else
                    dev_su3MtV(gfsmem, &(sin_up[hoppos]), &(shelp1_up[0]));
                    dev_su3MtV(gfsmem, &(sin_dn[hoppos]), &(shelp1_dn[0]));	    
                  #endif
                #endif
              }
            #else            
              #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem));
              #else
              dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem));
              #endif
              #ifdef RELATIVISTIC_BASIS
                #ifdef USETEXTURE
                  dev_su3MtV_spintex_rel_down(gfsmem, hoppos, &(shelp1_up[0])); 
                  dev_su3MtV_spintex2_rel_down(gfsmem, hoppos, &(shelp1_dn[0])); 		  
                #else
                  dev_su3MtV_rel_down(gfsmem, &(sin_up[hoppos]), &(shelp1_up[0]));
                  dev_su3MtV_rel_down(gfsmem, &(sin_dn[hoppos]), &(shelp1_dn[0]));		  
                #endif 
              #else
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem, hoppos, &(shelp1_up[0])); 
                  dev_su3MtV_spintex2(gfsmem, hoppos, &(shelp1_dn[0])); 		  
                #else
                  dev_su3MtV(gfsmem, &(sin_up[hoppos]), &(shelp1_up[0]));
                  dev_su3MtV(gfsmem, &(sin_dn[hoppos]), &(shelp1_dn[0]));		  
                #endif 
              #endif
            #endif
            
            #ifdef RELATIVISTIC_BASIS
              dev_kappaP0_minus_relativistic(&(ssum_up[0]), &(shelp1_up[0]), dev_k0);
              dev_kappaP0_minus_relativistic(&(ssum_dn[0]), &(shelp1_dn[0]), dev_k0);	      
            #else
              //-kappa(r + gamma_mu)
              dev_kappaP0_minus(&(ssum_up[0]), &(shelp1_up[0]), dev_k0);
              dev_kappaP0_minus(&(ssum_dn[0]), &(shelp1_dn[0]), dev_k0);	      
            #endif




//l==3,z 
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,gfindex_site[pos], 3, gaugevol ,&(gfsmem));
            #else
            dev_reconstructgf_2vtexref(gf, gfindex_site[pos], 3, gaugevol ,&(gfsmem));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet  
              //dev_kappaP3_plus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
	      
              //dev_su3MtV_kappaP3_plus_spintex(gfsmem,hoppos, &(ssum_up[0]), dev_k3);
              //dev_su3MtV_kappaP3_plus_spintex2(gfsmem,hoppos, &(ssum_dn[0]), dev_k3);
	      dev_su3MtV_kappaP3_plus_spintex_ud(gfsmem,hoppos, &(ssum_up[0]), &(ssum_dn[0]), dev_k3);
	      
	    #else
              dev_su3MtV_kappaP3_plus(gfsmem,&(sin_up[hoppos]), &(ssum_up[0]), dev_k3);
              dev_su3MtV_kappaP3_plus(gfsmem,&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k3);	      
            #endif
            

//l==3,z               
            
            //negative direction
            hoppos = nn_evenodd[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 3, gaugevol ,&(gfsmem));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 3, gaugevol ,&(gfsmem));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet 
              //dev_kappaP3_minus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
	      
	      //dev_su3MtV_kappaP3_minus_spintex(gfsmem,hoppos, &(ssum_up[0]), dev_k3);
	      //dev_su3MtV_kappaP3_minus_spintex2(gfsmem,hoppos, &(ssum_dn[0]), dev_k3);	   
	      dev_su3MtV_kappaP3_minus_spintex_ud(gfsmem,hoppos,&(ssum_up[0]), &(ssum_dn[0]), dev_k3);
	    #else
	      dev_su3MtV_kappaP3_minus(gfsmem,&(sin_up[hoppos]), &(ssum_up[0]), dev_k3);
	      dev_su3MtV_kappaP3_minus(gfsmem,&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k3);	      
            #endif
            





//l==2,y 
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,gfindex_site[pos], 2, gaugevol ,&(gfsmem));
            #else
            dev_reconstructgf_2vtexref(gf,gfindex_site[pos], 2, gaugevol ,&(gfsmem));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet
              //dev_kappaP2_plus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
	      
	      dev_su3MtV_kappaP2_plus_spintex(gfsmem,hoppos, &(ssum_up[0]), dev_k2);
	      dev_su3MtV_kappaP2_plus_spintex2(gfsmem,hoppos, &(ssum_dn[0]), dev_k2);	      
	    #else
	      dev_su3MtV_kappaP2_plus(gfsmem,&(sin_up[hoppos]), &(ssum_up[0]), dev_k2);
	      dev_su3MtV_kappaP2_plus(gfsmem,&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k2);	      
            #endif
            


//l==2,y        

            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 2, gaugevol ,&(gfsmem));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 2, gaugevol ,&(gfsmem));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet
              //dev_kappaP2_minus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
	      
	      dev_su3MtV_kappaP2_minus_spintex(gfsmem,hoppos, &(ssum_up[0]), dev_k2);
	      dev_su3MtV_kappaP2_minus_spintex2(gfsmem,hoppos, &(ssum_dn[0]), dev_k2);	      
	    #else
	      dev_su3MtV_kappaP2_minus(gfsmem,&(sin_up[hoppos]), &(ssum_up[0]), dev_k2);
	      dev_su3MtV_kappaP2_minus(gfsmem,&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k2);	      
            #endif
            
            



//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,gfindex_site[pos], 1, gaugevol ,&(gfsmem));
            #else
            dev_reconstructgf_2vtexref(gf,gfindex_site[pos], 1, gaugevol ,&(gfsmem));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet
              //dev_kappaP1_plus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
	      
	      dev_su3MtV_kappaP1_plus_spintex(gfsmem,hoppos, &(ssum_up[0]), dev_k1);
	      dev_su3MtV_kappaP1_plus_spintex2(gfsmem,hoppos, &(ssum_dn[0]), dev_k1);	      
	    #else
	      dev_su3MtV_kappaP1_plus(gfsmem,&(sin_up[hoppos]), &(ssum_up[0]), dev_k1);
	      dev_su3MtV_kappaP1_plus(gfsmem,&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k1);	      
            #endif
            
            


//l==1,x 
            
            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 1, gaugevol ,&(gfsmem));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 1, gaugevol ,&(gfsmem));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet
              //dev_kappaP1_minus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
	      
	      dev_su3MtV_kappaP1_minus_spintex(gfsmem,hoppos, &(ssum_up[0]), dev_k1);
	      dev_su3MtV_kappaP1_minus_spintex2(gfsmem,hoppos, &(ssum_dn[0]), dev_k1);	      
	    #else
	      dev_su3MtV_kappaP1_minus(gfsmem,&(sin_up[hoppos]), &(ssum_up[0]), dev_k1);
	      dev_su3MtV_kappaP1_minus(gfsmem,&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k1);	      
            #endif
            
            
 
        //copy to output spinor
        //dev_write_spinor(&(ssum[0]),&(sout[pos])); 
	sout_up[pos+0*DEVOFF].x = ssum_up[0].x; sout_up[pos+0*DEVOFF].y = ssum_up[0].y; 
	sout_up[pos+0*DEVOFF].z = ssum_up[0].z; sout_up[pos+0*DEVOFF].w = ssum_up[0].w;	
	sout_up[pos+1*DEVOFF].x = ssum_up[1].x; sout_up[pos+1*DEVOFF].y = ssum_up[1].y; 
	sout_up[pos+1*DEVOFF].z = ssum_up[1].z; sout_up[pos+1*DEVOFF].w = ssum_up[1].w;
	sout_up[pos+2*DEVOFF].x = ssum_up[2].x; sout_up[pos+2*DEVOFF].y = ssum_up[2].y; 
	sout_up[pos+2*DEVOFF].z = ssum_up[2].z; sout_up[pos+2*DEVOFF].w = ssum_up[2].w;
	sout_up[pos+3*DEVOFF].x = ssum_up[3].x; sout_up[pos+3*DEVOFF].y = ssum_up[3].y; 
	sout_up[pos+3*DEVOFF].z = ssum_up[3].z; sout_up[pos+3*DEVOFF].w = ssum_up[3].w;
	sout_up[pos+4*DEVOFF].x = ssum_up[4].x; sout_up[pos+4*DEVOFF].y = ssum_up[4].y; 
	sout_up[pos+4*DEVOFF].z = ssum_up[4].z; sout_up[pos+4*DEVOFF].w = ssum_up[4].w;	
	sout_up[pos+5*DEVOFF].x = ssum_up[5].x; sout_up[pos+5*DEVOFF].y = ssum_up[5].y; 
	sout_up[pos+5*DEVOFF].z = ssum_up[5].z; sout_up[pos+5*DEVOFF].w = ssum_up[5].w;	

	sout_dn[pos+0*DEVOFF].x = ssum_dn[0].x; sout_dn[pos+0*DEVOFF].y = ssum_dn[0].y; 
	sout_dn[pos+0*DEVOFF].z = ssum_dn[0].z; sout_dn[pos+0*DEVOFF].w = ssum_dn[0].w;	
	sout_dn[pos+1*DEVOFF].x = ssum_dn[1].x; sout_dn[pos+1*DEVOFF].y = ssum_dn[1].y; 
	sout_dn[pos+1*DEVOFF].z = ssum_dn[1].z; sout_dn[pos+1*DEVOFF].w = ssum_dn[1].w;
	sout_dn[pos+2*DEVOFF].x = ssum_dn[2].x; sout_dn[pos+2*DEVOFF].y = ssum_dn[2].y; 
	sout_dn[pos+2*DEVOFF].z = ssum_dn[2].z; sout_dn[pos+2*DEVOFF].w = ssum_dn[2].w;
	sout_dn[pos+3*DEVOFF].x = ssum_dn[3].x; sout_dn[pos+3*DEVOFF].y = ssum_dn[3].y; 
	sout_dn[pos+3*DEVOFF].z = ssum_dn[3].z; sout_dn[pos+3*DEVOFF].w = ssum_dn[3].w;
	sout_dn[pos+4*DEVOFF].x = ssum_dn[4].x; sout_dn[pos+4*DEVOFF].y = ssum_dn[4].y; 
	sout_dn[pos+4*DEVOFF].z = ssum_dn[4].z; sout_dn[pos+4*DEVOFF].w = ssum_dn[4].w;	
	sout_dn[pos+5*DEVOFF].x = ssum_dn[5].x; sout_dn[pos+5*DEVOFF].y = ssum_dn[5].y; 
	sout_dn[pos+5*DEVOFF].z = ssum_dn[5].z; sout_dn[pos+5*DEVOFF].w = ssum_dn[5].w;		
	
  }
}














#ifdef HALF

//applies the Hopping Part Even-Odd  HALF PRECISION !
//else aequivalent to the above version 
__global__ void dev_Hopping_Matrix_half(const dev_su3_2v_half * gf, const dev_spinor_half * sin, const float* sin_norm, dev_spinor_half * sout, float* sout_norm, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo){

  int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];
    __shared__ dev_su3_pad gfsmem[BLOCK];


  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  
  //the volume for the gaugefield is 2x volume of spinors!!!
 //the volume for the gaugefield is 2x volume of spinors!!!
 #ifdef MPI
   int gaugevol = 2*dev_VOLUMEPLUSRAND;
 #else
   int gaugevol = 2*dev_VOLUME; 
 #endif
 
  if(pos < dev_VOLUME){
  

  dev_zero_spinor(&(ssum[0])); // zero sum        
 #ifdef TEMPORALGAUGE
  int spatialvol = dev_LX*dev_LY*dev_LZ;
 #endif


//hopping term                
//l==0,t
            //positive direction
            hoppos = nn_evenodd[8*pos];
             //hoppos = tex1Dfetch(nn_tex,8*pos);
            //color
            
            #ifdef TEMPORALGAUGE
              float norm;
	      int i;
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_site[pos]) < (dev_T-1)*spatialvol) || (dev_rank < dev_nproc-1) ) {
                //if ((gfindex_site[pos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_site[pos]/spatialvol) != (dev_T-1) ) {
              #endif
              
              #ifdef USETEXTURE
                norm = tex1Dfetch(spinnormhalf_tex, hoppos);
                shelp1[0] = tex1Dfetch(spinhalf_tex0,hoppos);
                shelp1[1] = tex1Dfetch(spinhalf_tex1,hoppos);
                shelp1[2] = tex1Dfetch(spinhalf_tex2,hoppos);
                #ifdef RELATIVISTIC_BASIS
                  shelp1[3].x = 0.0f; shelp1[3].y = 0.0f; shelp1[3].z = 0.0f; shelp1[3].w = 0.0f;
		  shelp1[4].x = 0.0f; shelp1[4].y = 0.0f; shelp1[4].z = 0.0f; shelp1[4].w = 0.0f;
		  shelp1[5].x = 0.0f; shelp1[5].y = 0.0f; shelp1[5].z = 0.0f; shelp1[5].w = 0.0f;
		#else
		  shelp1[3] = tex1Dfetch(spinhalf_tex3,hoppos);
                  shelp1[4] = tex1Dfetch(spinhalf_tex4,hoppos);
                  shelp1[5] = tex1Dfetch(spinhalf_tex5,hoppos);
		#endif
                //normalize
                #pragma unroll 6
                for(i=0; i<6; i++){
                  shelp1[i].x = norm*shelp1[i].x;
                  shelp1[i].y = norm*shelp1[i].y;
                  shelp1[i].z = norm*shelp1[i].z;
                  shelp1[i].w = norm*shelp1[i].w;
                }
              #else
                norm = sin_norm[hoppos];
                //read and normalize
                #ifdef RELATIVISTIC_BASIS
                  #pragma unroll 3
                  for(i=0; i<3; i++){
                    shelp1[i].x = norm*sh2fl(sin[hoppos+i*DEVOFF].x);
                    shelp1[i].y = norm*sh2fl(sin[hoppos+i*DEVOFF].y);
                    shelp1[i].z = norm*sh2fl(sin[hoppos+i*DEVOFF].z);
                    shelp1[i].w = norm*sh2fl(sin[hoppos+i*DEVOFF].w);
		  }
		  #pragma unroll 3
                  for(i=3; i<6; i++){
                    shelp1[i].x = 0.0f;
                    shelp1[i].y = 0.0f;
                    shelp1[i].z = 0.0f;
                    shelp1[i].w = 0.0f;
		  }
                #else
                  #pragma unroll 6
                  for(i=0; i<6; i++){
                    shelp1[i].x = norm*sh2fl(sin[hoppos+i*DEVOFF].x);
                    shelp1[i].y = norm*sh2fl(sin[hoppos+i*DEVOFF].y);
                    shelp1[i].z = norm*sh2fl(sin[hoppos+i*DEVOFF].z);
                    shelp1[i].w = norm*sh2fl(sin[hoppos+i*DEVOFF].w);
		  }
		#endif
              #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref_half(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem[ix].m));
                #else
                dev_reconstructgf_2vtexref_half(gf,gfindex_site[pos], 0, gaugevol ,&(gfsmem[ix].m));
                #endif
                #ifdef RELATIVISTIC_BASIS
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex_rel_up(gfsmem[ix].m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV_half_rel_up(gfsmem[ix].m, &(sin[hoppos]), &(sin_norm[hoppos]), &(shelp1[0]));
                  #endif
                #else
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV_half(gfsmem[ix].m, &(sin[hoppos]), &(sin_norm[hoppos]), &(shelp1[0]));
                  #endif
                #endif
              }
            #else
              #ifdef GF_8
              dev_reconstructgf_8texref_half(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem[ix].m));
              #else
              dev_reconstructgf_2vtexref_half(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem[ix].m));
              #endif
              #ifdef RELATIVISTIC_BASIS
                #ifdef USETEXTURE
                  dev_su3MtV_spintex_rel_up(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV_half_rel_up(gfsmem[ix].m, &(sin[hoppos]), &(sin_norm[hoppos]), &(shelp1[0]));
                #endif
              #else
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV_half(gfsmem[ix].m, &(sin[hoppos]), &(sin_norm[hoppos]), &(shelp1[0]));
                #endif
              #endif
            #endif
            
            #ifdef RELATIVISTIC_BASIS
              dev_kappaP0_plus_relativistic(&(ssum[0]), &(shelp1[0]), dev_cconj(dev_k0));
            #else
              //-kappa(r - gamma_mu)
              dev_kappaP0_plus(&(ssum[0]), &(shelp1[0]), dev_cconj(dev_k0));
            #endif
            
//l==0,t
            //negative direction
            hoppos = nn_evenodd[8*pos+4]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+4);
            //color
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_nextsite[hoppos]) < (dev_T-1)*spatialvol) || (dev_rank > 0) ) {
                //if ((gfindex_nextsite[hoppos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_nextsite[hoppos]/spatialvol) != (dev_T-1) ) {
              #endif
              
              #ifdef USETEXTURE
                norm = tex1Dfetch(spinnormhalf_tex, hoppos);
                #ifdef RELATIVISTIC_BASIS
                  shelp1[0].x = 0.0f; shelp1[0].y = 0.0f; shelp1[0].z = 0.0f; shelp1[0].w = 0.0f;
		  shelp1[1].x = 0.0f; shelp1[1].y = 0.0f; shelp1[1].z = 0.0f; shelp1[1].w = 0.0f;
		  shelp1[2].x = 0.0f; shelp1[2].y = 0.0f; shelp1[2].z = 0.0f; shelp1[2].w = 0.0f;
                #else
		  shelp1[0] = tex1Dfetch(spinhalf_tex0,hoppos);
                  shelp1[1] = tex1Dfetch(spinhalf_tex1,hoppos);
                  shelp1[2] = tex1Dfetch(spinhalf_tex2,hoppos);
                #endif
		shelp1[3] = tex1Dfetch(spinhalf_tex3,hoppos);
                shelp1[4] = tex1Dfetch(spinhalf_tex4,hoppos);
                shelp1[5] = tex1Dfetch(spinhalf_tex5,hoppos);

                //normalize
                #pragma unroll 6
                for(i=0; i<6; i++){
                  shelp1[i].x = norm*shelp1[i].x;
                  shelp1[i].y = norm*shelp1[i].y;
                  shelp1[i].z = norm*shelp1[i].z;
                  shelp1[i].w = norm*shelp1[i].w;
                }
              #else
                norm = sin_norm[hoppos];
                //read and normalize
                #ifdef RELATIVISTIC_BASIS
                #pragma unroll 3
                  for(i=0; i<3; i++){
                    shelp1[i].x = 0.0f;
                    shelp1[i].y = 0.0f;
                    shelp1[i].z = 0.0f;
                    shelp1[i].w = 0.0f;
                  } 
                  #pragma unroll 3
                  for(i=3; i<6; i++){
                    shelp1[i].x = norm*sh2fl(sin[hoppos+i*DEVOFF].x);
                    shelp1[i].y = norm*sh2fl(sin[hoppos+i*DEVOFF].y);
                    shelp1[i].z = norm*sh2fl(sin[hoppos+i*DEVOFF].z);
                    shelp1[i].w = norm*sh2fl(sin[hoppos+i*DEVOFF].w);
                  } 
                #else
                #pragma unroll 6
                for(i=0; i<6; i++){
                  shelp1[i].x = norm*sh2fl(sin[hoppos+i*DEVOFF].x);
                  shelp1[i].y = norm*sh2fl(sin[hoppos+i*DEVOFF].y);
                  shelp1[i].z = norm*sh2fl(sin[hoppos+i*DEVOFF].z);
                  shelp1[i].w = norm*sh2fl(sin[hoppos+i*DEVOFF].w);
                }
                #endif
              #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref_dagger_half(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem[ix].m));
                #else
                dev_reconstructgf_2vtexref_dagger_half(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem[ix].m));
                #endif
                #ifdef RELATIVISTIC_BASIS
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex_rel_down(gfsmem[ix].m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV_half_rel_down(gfsmem[ix].m, &(sin[hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
                  #endif
                #else
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV_half(gfsmem[ix].m, &(sin[hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
                  #endif
                #endif
              }
            #else 
            
              #ifdef GF_8
              dev_reconstructgf_8texref_dagger_half(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem[ix].m));
              #else
              dev_reconstructgf_2vtexref_dagger_half(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem[ix].m));
              #endif
              #ifdef RELATIVISTIC_BASIS
                #ifdef USETEXTURE
                  dev_su3MtV_spintex_rel_down(gfsmem[ix].m, hoppos, &(shelp1[0]));  
                #else
                  dev_su3MtV_half_rel_down(gfsmem[ix].m, &(sin[hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
                #endif 
              #else
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));  
                #else
                  dev_su3MtV_half(gfsmem[ix].m, &(sin[hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
                #endif 
              #endif
            
            #endif
            
            #ifdef RELATIVISTIC_BASIS
              //-kappa(r + gamma_mu)
              dev_kappaP0_minus_relativistic(&(ssum[0]), &(shelp1[0]), dev_k0);
            #else
              //-kappa(r + gamma_mu)
              dev_kappaP0_minus(&(ssum[0]), &(shelp1[0]), dev_k0);
            #endif




//l==3,z 
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_half(gf,gfindex_site[pos], 3, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_half(gf, gfindex_site[pos], 3, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)    
            dev_kappaP3_plus(&(ssum[0]), &(shelp1[0]), dev_k3.re);

//l==3,z               
            
            //negative direction
            hoppos = nn_evenodd[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger_half(gf,gfindex_nextsite[hoppos], 3, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger_half(gf, gfindex_nextsite[hoppos], 3, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_kappaP3_minus(&(ssum[0]), &(shelp1[0]), dev_k3.re);





//l==2,y 
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_half(gf,gfindex_site[pos], 2, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_half(gf,gfindex_site[pos], 2, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
	    dev_kappaP2_plus(&(ssum[0]), &(shelp1[0]), dev_k2.re);


//l==2,y        

            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger_half(gf,gfindex_nextsite[hoppos], 2, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger_half(gf,gfindex_nextsite[hoppos], 2, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[hoppos]), &(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_kappaP2_minus(&(ssum[0]), &(shelp1[0]), dev_k2.re);




//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_half(gf,gfindex_site[pos], 1, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_half(gf, gfindex_site[pos], 1, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_kappaP1_plus(&(ssum[0]), &(shelp1[0]), dev_k1.re);



//l==1,x 
            
            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger_half(gf,gfindex_nextsite[hoppos], 1, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger_half(gf,gfindex_nextsite[hoppos], 1, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_kappaP1_minus(&(ssum[0]), &(shelp1[0]), dev_k1.re);

 
        //write to output spinor and write the norm
        dev_write_spinor_half(&(ssum[0]),&(sout[pos]), &(sout_norm[pos])); 
  }
}


#endif









