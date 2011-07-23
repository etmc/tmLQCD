//prepare    cublas[typeId el. {s,S,d,D,c,C,z,Z}][function]([parameterList])  for template usage by defining
//overloaded cublas                              [Function]([parameterList])

#include "cublas.h"
#include "cudaglobal.h"


//#ifdef OLD_CUBLAS

  float  cublasDot(int n,const float* x,int incx,const float* y,int incy)
  { return cublasSdot(n,x,incx,y,incy); }
  double cublasDot(int n,const double* x,int incx,const double* y,int incy)
  { return cublasDdot(n,x,incx,y,incy); }
 
  void cublasAxpy(int n,float alpha,const float* x,int incx,float* y,int incy)
  { cublasSaxpy(n,alpha,x,incx,y,incy); }
  void cublasAxpy(int n,double alpha,const double* x,int incx,double* y,int incy)
  { cublasDaxpy(n,alpha,x,incx,y,incy); }

  void cublasScal(int n,float alpha,float* x,int incx)
  { cublasSscal(n,alpha,x,incx); }
  void cublasScal(int n,double alpha,double* x,int incx)
  { cublasDscal(n,alpha,x,incx); }

  void cublasCopy(int n,const float* x,int incx,float* y,int incy)
  { cublasScopy(n,x,incx,y,incy); }
  void cublasCopy(int n,const double* x,int incx,double* y,int incy)
  { cublasDcopy(n,x,incx,y,incy); }

/*#else
  
  template<class RealT> inline cublasStatus_t cublasDot        (cublasHandle_t handle,int n,const RealT* x,int incx,const RealT* y,int incy,RealT* result) 
  { return RealT.cublasWrapperError(); } //produces an error when called with wrong template type
  template<           > inline cublasStatus_t cublasDot<float >(cublasHandle_t handle,int n,const RealT* x,int incx,const RealT* y,int incy,RealT* result) 
  { return cublasSdot(handle,n,x,incx,y,incy,result); }
  template<           > inline cublasStatus_t cublasDot<double>(cublasHandle_t handle,int n,const RealT* x,int incx,const RealT* y,int incy,RealT* result) 
  { return cublasDdot(handle,n,x,incx,y,incy,result); }
#endif*/

