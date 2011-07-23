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
template<class RealT>
__device__ void dev_kappaP1_plus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa){

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
template<class RealT>
__device__ void dev_kappaP1_minus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa){

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
template<class RealT>
__device__ void dev_kappaP2_plus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa){


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
template<class RealT>
__device__ void dev_kappaP2_minus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa){


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
template<class RealT>
__device__ void dev_kappaP3_plus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa){

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
template<class RealT>
__device__ void dev_kappaP3_minus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa){

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
template<class RealT>
__device__ void dev_kappaP0_plus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, dev_complexM(RealT) kappa){


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






//-kappa(r - gamma_mu)
template<class RealT>
__device__ void dev_kappaP0_minus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, dev_complexM(RealT) kappa){


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
template<class RealT>
__global__ void dev_Hopping_Matrix(const dev_su3_2vM(RealT) * gf, const dev_spinorM(RealT) * sin, dev_spinorM(RealT) * sout, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo){

  int pos,hoppos;
    dev_spinorM(RealT)             shelp1[6], ssum[6];
    __shared__ dev_su3_padM(RealT) gfsmem[BLOCK];



  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  
  
  if(pos < dev_VOLUME){
  

  dev_zero_spinor<RealT>(&(ssum[0])); // zero sum        
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
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_site[pos]) < (dev_T-1)*spatialvol) || (dev_rank < dev_nproc-1) ) {
                //if ((gfindex_site[pos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_site[pos]/spatialvol) != (dev_T-1) ) {
              #endif
              
              #ifdef USETEXTURE
                shelp1[0] = tex1Dfetch(spin_tex,6*hoppos);
                shelp1[1] = tex1Dfetch(spin_tex,6*hoppos+1);
                shelp1[2] = tex1Dfetch(spin_tex,6*hoppos+2);
                shelp1[3] = tex1Dfetch(spin_tex,6*hoppos+3);
                shelp1[4] = tex1Dfetch(spin_tex,6*hoppos+4);
                shelp1[5] = tex1Dfetch(spin_tex,6*hoppos+5);
              #else
                shelp1[0] = sin[6*hoppos];
                shelp1[1] = sin[6*hoppos+1];
                shelp1[2] = sin[6*hoppos+2];
                shelp1[3] = sin[6*hoppos+3];
                shelp1[4] = sin[6*hoppos+4];
                shelp1[5] = sin[6*hoppos+5];
              #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref <RealT>(gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
                #else
                dev_reconstructgf_2vtexref<RealT>(gf,4*(gfindex_site[pos]),&(gfsmem[ix].m));
                #endif
                #ifdef USETEXTURE
                  dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV        <RealT>(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
                #endif
              }
            #else
              #ifdef GF_8
              dev_reconstructgf_8texref <RealT>(gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
              #else
              dev_reconstructgf_2vtexref<RealT>(gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
              #endif
              #ifdef USETEXTURE
                dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
              #else
                dev_su3MtV        <RealT>(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
              #endif
            #endif
            
            //-kappa(r - gamma_mu)
            #ifdef GF_8
             dev_kappaP0_plus(&(ssum[0]), &(shelp1[0]), dev_cconj(dev_complexT<RealT>(dev_k0)));
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));
              dev_Gamma0<RealT>(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k0),&(shelp1[0]), &(ssum[0]));
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
                shelp1[0] = tex1Dfetch(spin_tex,6*hoppos);
                shelp1[1] = tex1Dfetch(spin_tex,6*hoppos+1);
                shelp1[2] = tex1Dfetch(spin_tex,6*hoppos+2);
                shelp1[3] = tex1Dfetch(spin_tex,6*hoppos+3);
                shelp1[4] = tex1Dfetch(spin_tex,6*hoppos+4);
                shelp1[5] = tex1Dfetch(spin_tex,6*hoppos+5);
               #else
                shelp1[0] = sin[6*hoppos];
                shelp1[1] = sin[6*hoppos+1];
                shelp1[2] = sin[6*hoppos+2];
                shelp1[3] = sin[6*hoppos+3];
                shelp1[4] = sin[6*hoppos+4];
                shelp1[5] = sin[6*hoppos+5];
               #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref_dagger <RealT>(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
                #else
                dev_reconstructgf_2vtexref_dagger<RealT>(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
                #endif
                #ifdef USETEXTURE
                  dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV        <RealT>(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
                #endif 
              }
            #else            
              #ifdef GF_8
              dev_reconstructgf_8texref_dagger <RealT>(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
              #else
              dev_reconstructgf_2vtexref_dagger<RealT>(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
              #endif
              #ifdef USETEXTURE
                dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));  
              #else
                dev_su3MtV        <RealT>(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
              #endif 
            #endif
            
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP0_minus(&(ssum[0]), &(shelp1[0]), dev_complexT<RealT>(dev_k0));
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));
              dev_Gamma0<RealT>(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));
            #endif




//l==3,z 
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref <RealT>(gf,4*(gfindex_site[pos])+(3),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref<RealT>(gf, 4*(gfindex_site[pos])+(3),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)    
            #ifdef GF_8
            dev_kappaP3_plus(&(ssum[0]), &(shelp1[0]), RealT(dev_k3.re));
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
              dev_Gamma3<RealT>(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k3),&(shelp1[0]), &(ssum[0]));
	    #endif
//l==3,z               
            
            //negative direction
            hoppos = nn_evenodd[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger <RealT>(gf,4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger<RealT>(gf,4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP3_minus(&(ssum[0]), &(shelp1[0]), RealT(dev_k3.re));
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
              dev_Gamma3<RealT>(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
            #endif




//l==2,y 
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref <RealT>(gf,4*(gfindex_site[pos])+(2),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref<RealT>(gf,4*(gfindex_site[pos])+(2),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP2_plus(&(ssum[0]), &(shelp1[0]), RealT(dev_k2.re));
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));
              dev_Gamma2<RealT>(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k2),&(shelp1[0]), &(ssum[0]));
            #endif

//l==2,y        

            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger <RealT>(gf,4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger<RealT>(gf,4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP2_minus(&(ssum[0]), &(shelp1[0]), RealT(dev_k2.re));
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));
              dev_Gamma2<RealT>(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));
	    #endif



//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref <RealT>(gf,4*(gfindex_site[pos])+(1),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref<RealT>(gf,4*(gfindex_site[pos])+(1),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP1_plus(&(ssum[0]), &(shelp1[0]), RealT(dev_k1.re));
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));
              dev_Gamma1<RealT>(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k1),&(shelp1[0]), &(ssum[0]));
	    #endif


//l==1,x 
            
            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger <RealT>(gf,4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger<RealT>(gf,4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP1_minus(&(ssum[0]), &(shelp1[0]), RealT(dev_k1.re));
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));
              dev_Gamma1<RealT>(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));      
            #endif
 
        //copy to output spinor
        dev_copy_spinor<RealT,RealT>(&(ssum[0]),&(sout[6*pos])); 
  }
}



/*

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
__global__ void dev_Hopping_Matrix(dev_su3_2v * gf, dev_spinor * sin, dev_spinor * sout, int * gfindex_site,int* gfindex_nextsite, int * nn_evenodd, const int eo){

  int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];
    __shared__ dev_su3 gfsmem[BLOCK];
    

  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  if(pos < dev_VOLUME){

  dev_zero_spinor(&(ssum[0])); // zero sum        
//hopping term                
//l==0,t
            //positive direction
            hoppos = nn_evenodd[8*pos];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(4*(gfindex_site[pos]),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos]),&(gfsmem[ix]));
            #endif
            
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k0),&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = nn_evenodd[8*pos+4]; 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(4*gfindex_nextsite[hoppos],&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos],&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));     
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));


//l==3,z               
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(4*(gfindex_site[pos])+(3),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos])+(3),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k3),&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = nn_evenodd[8*pos+7]; 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
         
         
//l==2,y        
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(4*(gfindex_site[pos])+(2),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos])+(2),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k2),&(shelp1[0]), &(ssum[0]));
            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));


//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(4*(gfindex_site[pos])+(1),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos])+(1),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k1),&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));      
               
 
        //copy to output spinor
        dev_copy_spinor(&(ssum[0]),&(sout[6*pos])); 
  }
}



*/





#ifdef HALF

//applies the Hopping Part Even-Odd  HALF PRECISION !
//else aequivalent to the above version 
__global__ void dev_Hopping_Matrix_half(const dev_su3_2v_half * gf, const dev_spinor_half * sin, const float* sin_norm, dev_spinor_half * sout, float* sout_norm, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo){

  typedef REAL RealT;
  int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];
    __shared__ dev_su3_pad gfsmem[BLOCK];


  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  
  
  if(pos < dev_VOLUME){
  

  dev_zero_spinor<RealT>(&(ssum[0])); // zero sum        
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
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_site[pos]) < (dev_T-1)*spatialvol) || (dev_rank < dev_nproc-1) ) {
                //if ((gfindex_site[pos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_site[pos]/spatialvol) != (dev_T-1) ) {
              #endif
              
              #ifdef USETEXTURE
                double norm = tex1Dfetch(spinnormhalf_tex, hoppos);
                shelp1[0] = tex1Dfetch(spinhalf_tex,6*hoppos);
                shelp1[1] = tex1Dfetch(spinhalf_tex,6*hoppos+1);
                shelp1[2] = tex1Dfetch(spinhalf_tex,6*hoppos+2);
                shelp1[3] = tex1Dfetch(spinhalf_tex,6*hoppos+3);
                shelp1[4] = tex1Dfetch(spinhalf_tex,6*hoppos+4);
                shelp1[5] = tex1Dfetch(spinhalf_tex,6*hoppos+5);
                //normalize
                #pragma unroll 6
                for(int i=0; i<6; i++){
                  shelp1[i].x = norm*shelp1[i].x;
                  shelp1[i].y = norm*shelp1[i].y;
                  shelp1[i].z = norm*shelp1[i].z;
                  shelp1[i].w = norm*shelp1[i].w;
                }
              #else
                norm = sin_norm[hoppos];
                //read and normalize
                #pragma unroll 6
                for(i=0; i<6; i++){
                  shelp1[i].x = norm*sh2fl(sin[6*hoppos+i].x);
                  shelp1[i].y = norm*sh2fl(sin[6*hoppos+i].y);
                  shelp1[i].z = norm*sh2fl(sin[6*hoppos+i].z);
                  shelp1[i].w = norm*sh2fl(sin[6*hoppos+i].w);
                }
              #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref_half (gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
                #else
                dev_reconstructgf_2vtexref_half(gf,4*(gfindex_site[pos]),&(gfsmem[ix].m));
                #endif
                #ifdef USETEXTURE
                  dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV_half(gfsmem[ix].m, &(sin[6*hoppos]), &(sin_norm[hoppos]), &(shelp1[0]));
                #endif
              }
            #else
              #ifdef GF_8
              dev_reconstructgf_8texref_half (gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
              #else
              dev_reconstructgf_2vtexref_half(gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
              #endif
              #ifdef USETEXTURE
                dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
              #else
                dev_su3MtV_half(gfsmem[ix].m, &(sin[6*hoppos]), &(sin_norm[hoppos]), &(shelp1[0]));
              #endif
            #endif
            
            //-kappa(r - gamma_mu)
            #ifdef GF_8
             dev_kappaP0_plus(&(ssum[0]), &(shelp1[0]), dev_cconj(dev_complexT<RealT>(dev_k0)));
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));
              dev_Gamma0<RealT>(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k0),&(shelp1[0]), &(ssum[0]));
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
                double norm = tex1Dfetch(spinnormhalf_tex, hoppos);
                shelp1[0] = tex1Dfetch(spinhalf_tex,6*hoppos);
                shelp1[1] = tex1Dfetch(spinhalf_tex,6*hoppos+1);
                shelp1[2] = tex1Dfetch(spinhalf_tex,6*hoppos+2);
                shelp1[3] = tex1Dfetch(spinhalf_tex,6*hoppos+3);
                shelp1[4] = tex1Dfetch(spinhalf_tex,6*hoppos+4);
                shelp1[5] = tex1Dfetch(spinhalf_tex,6*hoppos+5);
                //normalize
                #pragma unroll 6
                for(int i=0; i<6; i++){
                  shelp1[i].x = norm*shelp1[i].x;
                  shelp1[i].y = norm*shelp1[i].y;
                  shelp1[i].z = norm*shelp1[i].z;
                  shelp1[i].w = norm*shelp1[i].w;
                }
              #else
                norm = sin_norm[hoppos];
                //read and normalize
                #pragma unroll 6
                for(i=0; i<6; i++){
                  shelp1[i].x = norm*sh2fl(sin[6*hoppos+i].x);
                  shelp1[i].y = norm*sh2fl(sin[6*hoppos+i].y);
                  shelp1[i].z = norm*sh2fl(sin[6*hoppos+i].z);
                  shelp1[i].w = norm*sh2fl(sin[6*hoppos+i].w);
                }
              #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref_dagger_half (gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
                #else
                dev_reconstructgf_2vtexref_dagger_half(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
                #endif
                #ifdef USETEXTURE
                  dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV_half(gfsmem[ix].m, &(sin[6*hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
                #endif 
              }
            #else            
              #ifdef GF_8
              dev_reconstructgf_8texref_dagger_half (gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
              #else
              dev_reconstructgf_2vtexref_dagger_half(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
              #endif
              #ifdef USETEXTURE
                dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));  
              #else
                dev_su3MtV_half(gfsmem[ix].m, &(sin[6*hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
              #endif 
            #endif
            
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP0_minus(&(ssum[0]), &(shelp1[0]), dev_complexT<RealT>(dev_k0));
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));
              dev_Gamma0<RealT>(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));
            #endif




//l==3,z 
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_half (gf,4*(gfindex_site[pos])+(3),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_half(gf, 4*(gfindex_site[pos])+(3),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[6*hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)    
            #ifdef GF_8
            dev_kappaP3_plus(&(ssum[0]), &(shelp1[0]), RealT(dev_k3.re));
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
              dev_Gamma3<RealT>(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k3),&(shelp1[0]), &(ssum[0]));
            #endif
//l==3,z               
            
            //negative direction
            hoppos = nn_evenodd[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger_half (gf,4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger_half(gf,4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[6*hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP3_minus(&(ssum[0]), &(shelp1[0]), RealT(dev_k3.re));
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
              dev_Gamma3<RealT>(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
            #endif




//l==2,y 
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_half (gf,4*(gfindex_site[pos])+(2),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_half(gf,4*(gfindex_site[pos])+(2),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[6*hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP2_plus(&(ssum[0]), &(shelp1[0]), RealT(dev_k2.re));
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));
              dev_Gamma2<RealT>(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k2),&(shelp1[0]), &(ssum[0]));
            #endif

//l==2,y        

            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger_half (gf,4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger_half(gf,4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[6*hoppos]), &(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP2_minus(&(ssum[0]), &(shelp1[0]), RealT(dev_k2.re));
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));
              dev_Gamma2<RealT>(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));
            #endif



//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_half (gf,4*(gfindex_site[pos])+(1),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_half(gf,4*(gfindex_site[pos])+(1),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[6*hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP1_plus(&(ssum[0]), &(shelp1[0]), RealT(dev_k1.re));
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));
              dev_Gamma1<RealT>(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k1),&(shelp1[0]), &(ssum[0]));
            #endif


//l==1,x 
            
            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger_half (gf,4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger_half(gf,4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV_half(gfsmem[ix].m, &(sin[6*hoppos]),&(sin_norm[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP1_minus(&(ssum[0]), &(shelp1[0]), RealT(dev_k1.re));
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));
              dev_Gamma1<RealT>(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));      
            #endif
 
        //write to output spinor and write the norm
        dev_write_spinor_half(&(ssum[0]),&(sout[6*pos]), &(sout_norm[pos])); 
  }
}


#endif









