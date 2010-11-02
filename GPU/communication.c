
// this was taken from MPI.cuh
// and preliminarily put here for individual compilation



#ifdef HAVE_CONFIG_H
  #include<config.h>
#endif


#ifdef MPI
  
  #include <cublas.h>
  #include <mpi.h>
  // #include "cudaglobal.h"
  
  
  // a wrapper function for cublasSdot() (with the same interface)
  // provides the MPI communication via MPI_Allreduce()
  
  float cublasSdot_wrapper2(int size, float * A, int incx, float * B, int incy) {
  
    float result;
    float buffer;
    
    buffer = cublasSdot(size, (float *) A, incx, (float *) B, incy);
    MPI_Allreduce(&buffer, &result, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
    return(result);
  
  }
  
  
#endif
