



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


__device__ double dev_su3Retrace(dev_su3 * M){
  double erg;
  int i;
  erg = 0.0;
  for(i=0; i<3; i++){
    erg = erg + (*M)[i][i].re;
  }
return erg;
}



