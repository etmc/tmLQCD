#include "cudadefs.h"


#ifndef __CUDADEFS_H
 #define  __CUDADEFS_H

#ifndef __cplusplus
 #error "GPU code needs C++ due to templates"
#endif

#ifdef OLD_CUBLAS
 #define cublasStatus_t cublasStatus
#else
 #define cublasStatus cublasStatus_t 
#endif
/* GPU Stuff */
template<class RealT> struct dev_complexT
{
  RealT re;
  RealT im;

  template<class TargetRealT>operator dev_complexT<TargetRealT>() const //enables conversions between dev_complexT instances with different template arguments
  {
    dev_complexT<TargetRealT> tmp;
    tmp.re=TargetRealT(re); tmp.im=TargetRealT(im);

    return tmp;
  }
};
//template<>operator dev_complexT::dev_complexT<RealT>() const {return *this}
#define dev_complexM(RealT) dev_complexT<RealT>
#define dev_complex  dev_complexT<REAL >
#define dev_complexD dev_complexT<REALD>
//typedef dev_complexT<REAL > dev_complex ;
//typedef dev_complexT<REALD> dev_complexD;


/* non-scalar types x: usage xT<RealT>::type */

template<class RealT>struct REAL4T       //usefull to select REAL4-Type according to RealT
{ 
  struct type  {  RealT  w,x,y,z;  };
};
template<>           struct REAL4T<REAL> //template specialisation to select the nvidia-based type, which may offer some optimisation 
{  typedef REAL4  type;  };
template<>           struct REAL4T<REALD> 
{  typedef REAL4D type;  };
#define REAL4M(RealT) typename REAL4T<RealT>::type


/* Device Gauge Fields */
// Typedef  dev_su3 [3][3];         /* su(3)-Matrix 3x3 komplexe Einträge DEVICE */

template<class RealT> struct dev_su3T //no template typedef in c++ yet; structure will be usefull in function templates
{  typedef dev_complexT<RealT> type[3][3];  };/* su(3)-Matrix 3x3 komplexe Einträge DEVICE */
#define dev_su3M(RealT) typename dev_su3T<RealT>::type
#define dev_su3  dev_su3T<REAL >::type
#define dev_su3D dev_su3T<REALD>::type
//typedef dev_su3T<REAL >::type dev_su3;
//typedef dev_su3T<REALD>::type dev_su3D;

template<typename RealT> struct dev_su3_padT
{
  struct type //only for conistency
  {
    typename dev_su3T<RealT>::type m;
    RealT pad;
  };
};
#define dev_su3_padM(RealT) typename dev_su3_padT<RealT>::type
#define dev_su3_pad  dev_su3_padT<REAL >::type
#define dev_su3_padD dev_su3_padT<REALD>::type
//typedef dev_su3_padT<REAL >::type dev_su3_pad ;
//typedef dev_su3_padT<REALD>::type dev_su3_padD;

//#define su3_2vT REAL4T /* 2 Zeilen der su(3)-Matrix, 6 komplexe Einträge HOST 3*4*VOLUME in array -> texture */
template<class RealT> struct su3_2vT:REAL4T<RealT> {};
#define su3_2vM(RealT) typename su3_2vT<RealT>::type
#define su3_2v  su3_2vT<REAL >::type
#define su3_2vD su3_2vT<REALD>::type
//typedef su3_2vT<REAL >::type su3_2v ;
//typedef su3_2vT<REALD>::type su3_2vD;

//#define dev_su3_2vT REAL4T /* 2 Zeilen der su(3)-Matrix 3*2 komplexe Einträge DEVICE 3*4*VOLUME in array -> texture*/
template<class RealT> struct dev_su3_2vT:REAL4T<RealT> {};
#define dev_su3_2vM(RealT) typename dev_su3_2vT<RealT>::type
#define dev_su3_2v  dev_su3_2vT<REAL >::type
#define dev_su3_2vD dev_su3_2vT<REALD>::type
//typedef dev_su3_2vT<REAL >::type dev_su3_2v ;
//typedef dev_su3_2vT<REALD>::type dev_su3_2vD;
              
//#define dev_su3_8T REAL4T /* 8 numbers to reconstruct the gauge field as described in M. Clark */
template<class RealT> struct dev_su3_8T:REAL4T<RealT> {};
#define dev_su3_8M(RealT) typename dev_su3_8T<RealT>::type
#define dev_su3_8  dev_su3_8T<REAL>::type
#define dev_su3_8D dev_su3_8T<REAL>::type
//typedef dev_su3_8T<REAL >::type dev_su3_8;
//typedef dev_su3_8T<REALD>::type dev_su3_8D;


/* Device Spinor Fields */
//#define dev_spinorT REAL4T
template<class RealT> struct dev_spinorT:REAL4T<RealT> {};
#define dev_spinorM(RealT) typename dev_spinorT<RealT>::type
#define dev_spinor  dev_spinorT<REAL >::type
#define dev_spinorD dev_spinorT<REALD>::type
//typedef REAL4T<REAL >::type dev_spinor;
//typedef REAL4T<REALD>::type dev_spinorD;

template<class RealT> struct dev_spinor_smemT
{
  struct type
  {
    dev_spinorT<RealT> spin;
    RealT dummy; // used to fit memory usage to GPU architecture? - then we probably need template specialisation - otherwise delete this comment
  };
};
#define dev_spinor_smemM(RealT) typename dev_spinor_smemT<RealT>::type
#define dev_spinor_smem  dev_spinor_smemT<REAL >::type
#define dev_spinor_smemD dev_spinor_smemT<REALD>::type
//typedef dev_spinor_smemT<REAL >::type dev_spinor_smem ;
//typedef dev_spinor_smemT<REALD>::type dev_spinor_smemD;

template<class RealT> struct         dev_propmatrixT {  typedef dev_complexT<RealT> type[12][12];  };
#define dev_propmatrixM(RealT) typename dev_propmatrixT<RealT>::type
#define dev_propmatrix  dev_propmatrixT<REAL >::type
#define dev_propmatrixD dev_propmatrixT<REALD>::type
//typedef dev_propmatrixT<REAL >::type dev_propmatrix ;
//typedef dev_propmatrixT<REALD>::type dev_propmatrixD;

template<class RealT> struct dev_fbyfT {  typedef dev_complexT<RealT> type[4][4];  };
#define dev_fbyfM(RealT) typename dev_fbyfT<RealT>::type
#define dev_fbyf  dev_fbyfT<REAL >::type
#define dev_fbyfD dev_fbyfT<REALD>::type
//typedef dev_fbyfT<REAL >::type dev_fbyf ;
//typedef dev_fbyfT<REALD>::type dev_fbyfD;


#ifdef HALF
 typedef short4 dev_spinor_half;
 typedef short4 dev_su3_2v_half;
 typedef short4 dev_su3_8_half;
#endif


#endif
/* END GPU Stuff */
