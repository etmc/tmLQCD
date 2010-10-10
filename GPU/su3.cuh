



// u = v * w
__device__ void dev_su3_ti_su3(dev_su3* u, dev_su3 * v, dev_su3 * w){
  dev_complex help1, help2;
  dev_complex zero = dev_initcomplex(0.0f,0.0f);
  int i,j,k;
  for(i=0; i<3;i++){
    for(j=0; j<3; j++){
    
      help2 = zero;
      for(k=0; k<3; k++){
          help1 = dev_cmult((*v)[i][k],(*w)[k][j]);
          help2 = dev_cadd(help1, help2);
        }
        (*u)[i][j] = help2;    
    }
  }
}


__device__ float dev_su3Retrace(dev_su3 * M){
  float erg;
  int i;
  erg = 0.0;
  for(i=0; i<3; i++){
    erg = erg + (*M)[i][i].re;
  }
return erg;
}


__device__ void dev_su3trace(float2 * erg, dev_su3 * M){
  int i;
  float2 tr;
  tr.x = 0.0; tr.y = 0.0;
  for(i=0; i<3; i++){
    tr.x += (*M)[i][i].re;
    tr.y += (*M)[i][i].im;
  }
  (*erg).x = tr.x;
  (*erg).y = tr.y;
}



__device__ void dev_unit_su3(dev_su3 * M){
  (*M)[0][0].re = 1.0;
  (*M)[0][0].im = 0.0;
  (*M)[0][1].re = 0.0;
  (*M)[0][1].im = 0.0;
  (*M)[0][2].re = 0.0;
  (*M)[0][2].im = 0.0;
  
  (*M)[1][0].re = 0.0;
  (*M)[1][0].im = 0.0;
  (*M)[1][1].re = 1.0;
  (*M)[1][1].im = 0.0;
  (*M)[1][2].re = 0.0;
  (*M)[1][2].im = 0.0;
  
  (*M)[2][0].re = 0.0;
  (*M)[2][0].im = 0.0;
  (*M)[2][1].re = 0.0;
  (*M)[2][1].im = 0.0;
  (*M)[2][2].re = 1.0;
  (*M)[2][2].im = 0.0;
}




// to := from
__device__ void dev_su3_assign( dev_su3 * to, dev_su3 * from){
  int i,j;
  for(i=0; i<3; i++){
    for(j=0; j<3; j++){
      (*to)[i][j] = (*from)[i][j];
    }
  }
return;
}





