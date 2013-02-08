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
__device__ void dev_kappaP1_plus(dev_spinor * out, dev_spinor * in, REAL kappa){

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
__device__ void dev_kappaP1_minus(dev_spinor * out, dev_spinor * in, REAL kappa){

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
__device__ void dev_kappaP2_plus(dev_spinor * out, dev_spinor * in, REAL kappa){


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
__device__ void dev_kappaP2_minus(dev_spinor * out, dev_spinor * in, REAL kappa){


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
__device__ void dev_kappaP3_plus(dev_spinor * out, dev_spinor * in, REAL kappa){

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
__device__ void dev_kappaP3_minus(dev_spinor * out, dev_spinor * in, REAL kappa){

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



#define NEW2
#ifdef NEW2


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
__global__ void dev_Hopping_Matrix(const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo, int start, int size){

  int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];


  #ifdef GPU_3DBLOCK
    dev_su3_pad gfsmem;  
    pos = start  
          + (threadIdx.z + blockDim.z*(threadIdx.y + blockDim.y*(threadIdx.x))) 
          + blockDim.z*blockDim.y*blockDim.x*blockIdx.x;
    int ix = threadIdx.z + blockDim.z*(threadIdx.y + blockDim.y*(threadIdx.x));       
  #else
    dev_su3_pad gfsmem;  
    pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;  
    int ix = threadIdx.x;
  #endif

  
  
  if(pos < (start + size)){
 
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
	      
              dev_su3MtV_kappaP3_plus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k3.re);	      
	    #else
              dev_su3MtV_kappaP3_plus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k3.re);
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
	      
	      dev_su3MtV_kappaP3_minus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k3.re);	      
	    #else
	      dev_su3MtV_kappaP3_minus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k3.re);
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
	      
	      dev_su3MtV_kappaP2_plus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k2.re);	      
	    #else
	      dev_su3MtV_kappaP2_plus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k2.re);
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
	      
	      dev_su3MtV_kappaP2_minus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k2.re);	      
	    #else
	      dev_su3MtV_kappaP2_minus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k2.re);
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
	      
	      dev_su3MtV_kappaP1_plus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k1.re);	      
	    #else
	      dev_su3MtV_kappaP1_plus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k1.re);
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
	      
	      dev_su3MtV_kappaP1_minus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k1.re);	      
	    #else
	      dev_su3MtV_kappaP1_minus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k1.re);
            #endif
            
            
 
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


#else

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
__global__ void dev_Hopping_Matrix(const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo, int start, int size){

  int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];


  #ifdef GPU_3DBLOCK
    __shared__ dev_su3_pad gfsmem[BLOCK*BLOCKSUB*BLOCKSUB];  
    pos = start  
          + (threadIdx.z + blockDim.z*(threadIdx.y + blockDim.y*(threadIdx.x))) 
          + blockDim.z*blockDim.y*blockDim.x*blockIdx.x;
    int ix = threadIdx.z + blockDim.z*(threadIdx.y + blockDim.y*(threadIdx.x));       
  #else
    __shared__ dev_su3_pad gfsmem[BLOCK];  
    pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;  
    int ix = threadIdx.x;
  #endif

  
  
  if(pos < (start + size)){
 
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
                dev_reconstructgf_8texref(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem[ix].m));
                #else
                dev_reconstructgf_2vtexref(gf,gfindex_site[pos], 0, gaugevol ,&(gfsmem[ix].m));
                #endif
                
                #ifdef RELATIVISTIC_BASIS
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex_rel_up(gfsmem[ix].m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV_rel_up(gfsmem[ix].m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #else
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV(gfsmem[ix].m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #endif
              }
            #else
              #ifdef GF_8
              dev_reconstructgf_8texref(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem[ix].m));
              #else
              dev_reconstructgf_2vtexref(gf, gfindex_site[pos], 0, gaugevol ,&(gfsmem[ix].m));
              #endif
              #ifdef RELATIVISTIC_BASIS
                #ifdef USETEXTURE
                  dev_su3MtV_spintex_rel_up(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV_rel_up(gfsmem[ix].m, &(sin[hoppos]), &(shelp1[0]));
                #endif
              #else
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV(gfsmem[ix].m, &(sin[hoppos]), &(shelp1[0]));
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
                dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem[ix].m));
                #else
                dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem[ix].m));
                #endif
                #ifdef RELATIVISTIC_BASIS
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex_rel_down(gfsmem[ix].m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV_rel_down(gfsmem[ix].m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #else
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV(gfsmem[ix].m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #endif
              }
            #else            
              #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem[ix].m));
              #else
              dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 0, gaugevol ,&(gfsmem[ix].m));
              #endif
              #ifdef RELATIVISTIC_BASIS
                #ifdef USETEXTURE
                  dev_su3MtV_spintex_rel_down(gfsmem[ix].m, hoppos, &(shelp1[0]));  
                #else
                  dev_su3MtV_rel_down(gfsmem[ix].m, &(sin[hoppos]), &(shelp1[0]));
                #endif 
              #else
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));  
                #else
                  dev_su3MtV(gfsmem[ix].m, &(sin[hoppos]), &(shelp1[0]));
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
            dev_reconstructgf_8texref(gf,gfindex_site[pos], 3, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref(gf, gfindex_site[pos], 3, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet  
              //dev_kappaP3_plus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
	      
	      dev_su3MtV_kappaP3_plus_spintex(gfsmem[ix].m,hoppos, &(ssum[0]), dev_k3.re);
	    #else
              dev_su3MtV_kappaP3_plus(gfsmem[ix].m,&(sin[hoppos]), &(ssum[0]), dev_k3.re);
            #endif
            

//l==3,z               
            
            //negative direction
            hoppos = nn_evenodd[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 3, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 3, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet 
              //dev_kappaP3_minus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
	      
	      dev_su3MtV_kappaP3_minus_spintex(gfsmem[ix].m,hoppos, &(ssum[0]), dev_k3.re);
	    #else
	      dev_su3MtV_kappaP3_minus(gfsmem[ix].m,&(sin[hoppos]), &(ssum[0]), dev_k3.re);
            #endif
            





//l==2,y 
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,gfindex_site[pos], 2, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref(gf,gfindex_site[pos], 2, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet
              //dev_kappaP2_plus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
	      
	      dev_su3MtV_kappaP2_plus_spintex(gfsmem[ix].m,hoppos, &(ssum[0]), dev_k2.re);	      
	    #else
	      dev_su3MtV_kappaP2_plus(gfsmem[ix].m,&(sin[hoppos]), &(ssum[0]), dev_k2.re);
            #endif
            


//l==2,y        

            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 2, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 2, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet
              //dev_kappaP2_minus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
	      
	      dev_su3MtV_kappaP2_minus_spintex(gfsmem[ix].m,hoppos, &(ssum[0]), dev_k2.re);	      
	    #else
	      dev_su3MtV_kappaP2_minus(gfsmem[ix].m,&(sin[hoppos]), &(ssum[0]), dev_k2.re);
            #endif
            
            



//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,gfindex_site[pos], 1, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref(gf,gfindex_site[pos], 1, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet
              //dev_kappaP1_plus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
	      
	      dev_su3MtV_kappaP1_plus_spintex(gfsmem[ix].m,hoppos, &(ssum[0]), dev_k1.re);	      
	    #else
	      dev_su3MtV_kappaP1_plus(gfsmem[ix].m,&(sin[hoppos]), &(ssum[0]), dev_k1.re);
            #endif
            
            


//l==1,x 
            
            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,gfindex_nextsite[hoppos], 1, gaugevol ,&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,gfindex_nextsite[hoppos], 1, gaugevol ,&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet
              //dev_kappaP1_minus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
	      
	      dev_su3MtV_kappaP1_minus_spintex(gfsmem[ix].m,hoppos, &(ssum[0]), dev_k1.re);	      
	    #else
	      dev_su3MtV_kappaP1_minus(gfsmem[ix].m,&(sin[hoppos]), &(ssum[0]), dev_k1.re);
            #endif
            
            
 
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



#endif




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









