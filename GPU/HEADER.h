
#ifndef _MIXED_SOLVE_H_




///////////////////
// own functions //
///////////////////


// eo, nd 

void to_device (dev_spinor * device, spinor * host, dev_spinor * auxiliary, int size);

void to_host (spinor * host, dev_spinor * device, dev_spinor * auxiliary, int size);

__global__ void he_cg_init_nd_additional (float param_mubar, float param_epsbar);

__global__ void dev_mul_one_pm_imubar_gamma5 (dev_spinor * sin, dev_spinor * sout, REAL sign);

void init_mixedsolve_eo_nd(su3** gf);

void finalize_mixedsolve_eo_nd(void);

void matrix_multiplication32 (dev_spinor * , dev_spinor * , dev_spinor * , dev_spinor * , int, int, int, int, int, int, int, int);

void flopcount(unsigned long long int& total, int add);

extern "C" void benchmark_eo_nd (spinor * const Q_up, spinor * const Q_dn, int N);

int cg_eo_nd (dev_su3_2v * gf,
              dev_spinor * P_up, dev_spinor * P_dn,
              dev_spinor * Q_up, dev_spinor * Q_dn,
              int max_iter,
              int check_abs , int check_rel,
              double eps_abs, double eps_rel       );

extern "C" int mixedsolve_eo_nd (spinor * P_up, spinor * P_dn,
                                 spinor * Q_up, spinor * Q_dn,
                                 int max_iter, double eps_sq, int rel_prec);

void set_global_sizes();

float cublasDot_wrapper(int size, float * A, int incx, float * B, int incy);


// eo, nd, MPI

void convert2double_spin_mpi (dev_spinor* spin, spinor* h2d, int start, int end);
void convert2REAL4_spin_mpi (spinor* spin, dev_spinor* h2d, int start, int end);
void to_device_mpi (dev_spinor * device, spinor * host, dev_spinor * auxiliary, int size, int start, int end);
void to_host_mpi (spinor * host, dev_spinor * device, dev_spinor * auxiliary, int size, int start, int end);
void xchange_field_wrapper (dev_spinor * dev_spin, int ieo);
void Hopping_Matrix_wrapper (int ieo, dev_spinor * out, dev_spinor * in);
void su3to2vf4_mpi(su3** gf, dev_su3_2v* h2d_gf);
void su3to8_mpi(su3** gf, dev_su3_8* h2d_gf);
void init_iseven();
void init_nnspinor_eo_mpi();
void init_idxgauge_mpi();
void init_gpu_indexfields();
void free_gpu_indexfields();
__global__ void he_cg_init_nd_additional_mpi (int param_VOLUMEPLUSRAND, int param_RAND, int rank, int nproc);
void init_mixedsolve_eo_nd_mpi(su3** gf);
void finalize_mixedsolve_eo_nd_mpi(void);

__global__ void dev_Hopping_Matrix_mpi (const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout,
                                        int * dev_iup, int * dev_idn, int * dev_eo2lexic, int * dev_lexic2eosub,
                                        int ieo);

void matrix_multiplication32_mpi (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                                  dev_spinor * spinin_up , dev_spinor * spinin_dn ,
                                  int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                                  int gridsize3, int blocksize3, int gridsize4, int blocksize4);

int cg_eo_nd_mpi (dev_su3_2v * gf,
                  dev_spinor * P_up, dev_spinor * P_dn,
                  dev_spinor * Q_up, dev_spinor * Q_dn,
                  int max_iter,
                  int check_abs , int check_rel,
                  double eps_abs, double eps_rel       );

extern "C" int mixedsolve_eo_nd_mpi (spinor * P_up, spinor * P_dn,
                                     spinor * Q_up, spinor * Q_dn,
                                     int max_iter, double eps_sq, int rel_prec);




// ASYNC

__global__ void dev_Hopping_Matrix_ASYNC (const dev_su3_2v * gf, 
                                          const dev_spinor * sin, dev_spinor * sout,
                                          const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd,
                                          const int eo,
                                          int start, int size);

void HOPPING_ASYNC (dev_su3_2v * gf, 
                    dev_spinor * spinin, dev_spinor * spinout,
                    int * gfindex_site, int * gfindex_nextsite, int * nn_evenodd,
                    int ieo,
                    int gridsize, int blocksize);

void matrix_multiplication32_mpi_ASYNC (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                                        dev_spinor * spinin_up , dev_spinor * spinin_dn ,
                                        int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                                        int gridsize3, int blocksize3, int gridsize4, int blocksize4);









/////////////////////////
// Florian's functions //
/////////////////////////

// eo, non-eo, non-nd
template<class RealT>struct MixedsolveParameter;
template<class RealT>class MixedsolveOperator // interface class
{
public:
  virtual ~MixedsolveOperator();

  virtual void gpuInit  (dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim);
  virtual void gpu      (dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim)=0;  //linkage error due to pure virtual functions with gcc 4.3.2
  virtual void gpuDeinit(dev_spinorM(RealT)* spininout,dev_spinorM(RealT)* spinTmp,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim,const RealT scaleparam);

  virtual void checkInit  (spinor* const spinin,spinor* const spinTmp,spinor* const spinout,int volume);
  virtual void check      (spinor* const conjungateBasisPSpinin,spinor* const spinout,const int volume)=0;
  virtual void checkDeinit(spinor* const spinin,spinor* const spinTmp,spinor* const spinout,int volume);
};

// cublasWrapper
float  cublasDot(int n,const float* x,int incx,const float* y,int incy);
double cublasDot(int n,const double* x,int incx,const double* y,int incy);
void cublasAxpy(int n,float alpha,const float* x,int incx,float* y,int incy);
void cublasAxpy(int n,double alpha,const double* x,int incx,double* y,int incy);
void cublasScal(int n,float alpha,float* x,int incx);
void cublasScal(int n,double alpha,double* x,int incx);
void cublasCopy(int n,const float* x,int incx,float* y,int incy);
void cublasCopy(int n,const double* x,int incx,double* y,int incy);


//__device__ inline dev_complex dev_cconj (dev_complex c);
template<class RealT>__device__ inline dev_complexT<RealT> dev_cconj (dev_complexT<RealT> c);
//__device__ inline void dev_ccopy(dev_complex* von, dev_complex* nach);
template<class RealTVon,class RealTNach>__device__ inline void dev_ccopy(dev_complexT<RealTVon>* von, dev_complexT<RealTNach>* nach);
//__device__ inline REAL dev_cabssquare (dev_complex c);
template<class RealT>__device__ inline RealT dev_cabssquare (dev_complexT<RealT> c);
//__device__ inline REAL dev_cabsolute (dev_complex c);
template<class RealT>__device__ inline RealT dev_cabsolute (dev_complexT<RealT> c);
//__device__ inline  dev_complex dev_crealmult(dev_complex c1, REAL real);
template<class RealT>__device__ inline  dev_complexT<RealT> dev_crealmult(dev_complexT<RealT> c1, RealT real);
//__device__ inline dev_complex dev_cmult (dev_complex c1, dev_complex c2);
template<class RealT>__device__ inline dev_complexT<RealT> dev_cmult (dev_complexT<RealT> c1, dev_complexT<RealT> c2);
//__device__ inline dev_complex dev_cadd (dev_complex c1, dev_complex c2);
template<class RealT>__device__ inline dev_complexT<RealT> dev_cadd (dev_complexT<RealT> c1, dev_complexT<RealT> c2);
//__device__ inline dev_complex dev_cdiv(dev_complex c1, dev_complex c2);
template<class RealT>__device__ inline dev_complexT<RealT> dev_cdiv(dev_complexT<RealT> c1, dev_complexT<RealT> c2);
//__device__ inline dev_complex dev_csub(dev_complex c1, dev_complex c2);
template<class RealT>__device__ inline dev_complexT<RealT> dev_csub(dev_complexT<RealT> c1, dev_complexT<RealT> c2);
//__device__ inline dev_complex dev_initcomplex(REAL re, REAL im);
template<class RealT>__device__ inline dev_complexT<RealT> dev_initcomplex(RealT re, RealT im);
//__device__ inline void dev_copy_spinor(dev_spinor *i1, dev_spinor *i2);
template<class RealT1,class RealT2>__device__ inline void dev_copy_spinor(typename dev_spinorT<RealT1>::type *i1, typename dev_spinorT<RealT2>::type *i2);
//__device__ inline void dev_zero_spinor(dev_spinor *sin);
template<class RealT>__device__ inline void dev_zero_spinor(typename dev_spinorT<RealT>::type *sin);
//__device__ inline void dev_skalarmult_add_assign_spinor(dev_spinor *in, REAL lambda,dev_spinor * in2, dev_spinor * out);
template<class RealT>__device__ inline void dev_skalarmult_add_assign_spinor(typename dev_spinorT<RealT>::type *in, RealT lambda, typename dev_spinorT<RealT>::type * in2, typename dev_spinorT<RealT>::type * out);
//__device__ inline void dev_complexmult_add_assign_spinor(dev_spinor * in, dev_complex lambda,dev_spinor * in2, dev_spinor * out);
template<class RealT>__device__ inline void dev_complexmult_add_assign_spinor( typename dev_spinorT<RealT>::type* in, dev_complexT<RealT> lambda, typename dev_spinorT<RealT>::type* in2, typename dev_spinorT<RealT>::type* out);
//__device__ inline void dev_complexcgmult_add_assign_spinor(dev_spinor * in, dev_complex lambda,dev_spinor * in2, dev_spinor * out);
template<class RealT>__device__ inline void dev_complexcgmult_add_assign_spinor( typename dev_spinorT<RealT>::type * in, dev_complexT<RealT> lambda, typename dev_spinorT<RealT>::type* in2, typename dev_spinorT<RealT>::type* out);
//__device__ void inline dev_skalarmult_spinor(dev_spinor * in, dev_complex lambda, dev_spinor * out);
template<class RealT>__device__ void inline dev_skalarmult_spinora( typename dev_spinorT<RealT>::type* in, dev_complexT<RealT> lambda, typename dev_spinorT<RealT>::type* out);
//__device__ void inline dev_skalarmult_gamma5_spinor(dev_spinor * out, dev_complex lambda, dev_spinor * in);
template<class RealT>__device__ void inline dev_skalarmult_gamma5_spinor(typename dev_spinorT<RealT>::type* out, dev_complexT<RealT> lambda, typename dev_spinorT<RealT>::type* in);
//__device__ void inline dev_realmult_spinor(dev_spinor * in, REAL lambda);
template<class RealT>__device__ void inline dev_realmult_spinor(typename dev_spinorT<RealT>::type* in, RealT lambda);
//__device__ void inline dev_realmult_spinor_assign(dev_spinor* out, REAL lambda, dev_spinor* in);
template<class RealT>__device__ void inline dev_realmult_spinor_assign(typename dev_spinorT<RealT>::type* out, RealT lambda, typename dev_spinorT<RealT>::type* in);
//__device__ void dev_assign_realmult_add_spinor(dev_spinor* out, REAL lambda, dev_spinor* in1,  dev_spinor* in2);
template<class RealT>__device__ void dev_assign_realmult_add_spinor( typename dev_spinorT<RealT>::type* out, RealT lambda, typename dev_spinorT<RealT>::type* in1, typename dev_spinorT<RealT>::type* in2);
//__device__ inline void dev_add_spinor_assign(dev_spinor * i1, dev_spinor * i2);
template<class RealT>__device__ inline void dev_add_spinor_assign(typename dev_spinorT<RealT>::type * i1, typename dev_spinorT<RealT>::type * i2);
//__device__ inline void dev_sub_spinor_assign(dev_spinor * i1, dev_spinor * i2);
template<class RealT>__device__ inline void dev_sub_spinor_assign(typename dev_spinorT<RealT>::type * i1, typename dev_spinorT<RealT>::type * i2);
//__device__ void dev_su3MtV_spintex(dev_su3 M, int pos, dev_spinor * out);
template<class RealT> __device__ void dev_su3MtV_spintex(dev_su3M(RealT) M, int pos, dev_spinorM(RealT) * out);
//__device__ void dev_su3MtV(dev_su3 M, const dev_spinor * s, dev_spinor * out);
template<class RealT>__device__ void dev_su3MtV(typename dev_su3T<RealT>::type M, const typename dev_spinorT<RealT>::type * s, typename dev_spinorT<RealT>::type * out);
//__device__ void dev_su3MdaggertV(dev_su3 M, dev_spinor * s, dev_spinor * out);
template<class RealT>__device__ void dev_su3MdaggertV(typename dev_su3T<RealT>::type M, typename dev_spinorT<RealT>::type * s, typename dev_spinorT<RealT>::type * out);
//__device__ void dev_Gamma0(dev_spinor * in);
template<class RealT>__device__ void dev_Gamma0(typename dev_spinorT<RealT>::type * in);
//__device__ void dev_Gamma3(dev_spinor * in);
template<class RealT>__device__ void dev_Gamma3(typename dev_spinorT<RealT>::type * in);
//__device__ void dev_Gamma2(dev_spinor * in);
template<class RealT>__device__ void dev_Gamma2(typename dev_spinorT<RealT>::type * in);
//__device__ void dev_Gamma1(dev_spinor * in);
template<class RealT>__device__ void dev_Gamma1(typename dev_spinorT<RealT>::type * in);
//__device__ void dev_Gamma5(dev_spinor * in);
template<class RealT>__device__ void dev_Gamma5(typename dev_spinorT<RealT>::type * in);
//__device__ void dev_Gamma5_assign(dev_spinor* out, dev_spinor* in);
template<class RealT>__device__ void dev_Gamma5_assign(typename dev_spinorT<RealT>::type* out, typename dev_spinorT<RealT>::type* in);
//__device__ void dev_GammatV(int mu, dev_spinor * in);
template<class RealT>__device__ void dev_GammatV(int mu, typename dev_spinorT<RealT>::type * in);
//__device__ void dev_reconstructgf_2vtexref (const dev_su3_2v* field, int pos, dev_su3* gf);
//__device__ void dev_reconstructgf_2vtexref_dagger (const dev_su3_2v* field, int pos, dev_su3* gf);
//__device__ void dev_reconstructgf_8texref (const dev_su3_2v * field, int pos, dev_su3* gf);
//__device__ void dev_reconstructgf_8texref_dagger (const dev_su3_2v* field,int pos, dev_su3* gf);
template<class RealT> __device__ void dev_reconstructgf_2vtexref        (const typename dev_su3_2vT<RealT>::type* field, int pos, typename dev_su3T<RealT>::type* gf);
template<class RealT> __device__ void dev_reconstructgf_2vtexref_dagger (const typename dev_su3_2vT<RealT>::type* field, int pos, typename dev_su3T<RealT>::type* gf);
template<class RealT> __device__ void dev_reconstructgf_8texref         (const typename dev_su3_2vT<RealT>::type* field, int pos, typename dev_su3T<RealT>::type* gf);
template<class RealT> __device__ void dev_reconstructgf_8texref_dagger  (const typename dev_su3_2vT<RealT>::type* field, int pos, typename dev_su3T<RealT>::type* gf);
template<class RealT> __global__ void dev_gamma5(typename dev_spinorT<RealT>::type * sin, typename dev_spinorT<RealT>::type * sout);
__global__ void dev_swapmu();
//__global__ void dev_mul_one_pm_imu_inv(dev_spinor* sin, dev_spinor* sout, const REAL sign);
template<class RealT>__global__ void dev_mul_one_pm_imu_inv(dev_spinorM(RealT)* sin, dev_spinorM(RealT)* sout, const RealT sign);
//__global__ void dev_mul_one_pm_imu_sub_mul_gamma5(dev_spinor* sin1, dev_spinor* sin2, dev_spinor* sout, const REAL sign);
template<class RealT>__global__ void dev_mul_one_pm_imu_sub_mul_gamma5(dev_spinorM(RealT)* sin1, dev_spinorM(RealT)* sin2, dev_spinorM(RealT)* sout, const RealT sign);
template<class RealT>__device__ void dev_kappaP1_plus (dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa);
template<class RealT>__device__ void dev_kappaP1_minus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa);
template<class RealT>__device__ void dev_kappaP2_plus (dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa);
template<class RealT>__device__ void dev_kappaP2_minus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa);
template<class RealT>__device__ void dev_kappaP3_plus (dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa);
template<class RealT>__device__ void dev_kappaP3_minus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, RealT kappa);
template<class RealT>__device__ void dev_kappaP0_plus (dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, dev_complexM(RealT) kappa);
template<class RealT>__device__ void dev_kappaP0_minus(dev_spinorM(RealT) * out, dev_spinorM(RealT) * in, dev_complexM(RealT) kappa);
//__global__ void dev_Hopping_Matrix(const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo);
template<class RealT>__global__ void dev_Hopping_Matrix(const dev_su3_2vM(RealT) * gf, const dev_spinorM(RealT) * sin, dev_spinorM(RealT) * sout, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo);
//extern "C" void dev_Qtm_pm_psi(dev_spinor* spinin, dev_spinor* spinout, int gridsize, int blocksize, int gridsize2, int blocksize2);
template<class RealT>void dev_Qtm_pm_psi(dev_spinorM(RealT)* spinin, dev_spinorM(RealT)* spinout, int gridsize, int blocksize, int gridsize2, int blocksize2, MixedsolveParameter<RealT>& mixedsolveParameter);
//__global__ void dev_tm_dirac_kappa(dev_su3_2v * gf, dev_spinor * sin, dev_spinor * sout, int * dev_nn);
template<class RealT>
__global__ void dev_tm_dirac_kappa
(
  typename dev_su3_2vT<RealT>::type * gf,
  typename dev_spinorT<RealT>::type * sin,
  typename dev_spinorT<RealT>::type * sout,
  int * dev_nn
);
//extern "C" void dev_tm_dirac_dagger_kappa(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize);
template<class RealT>void dev_tm_dirac_dagger_kappa(typename dev_su3_2vT<RealT>::type * gf,typename dev_spinorT<RealT>::type* spinin,typename dev_spinorT<RealT>::type* spinout, int *grid, int * nn_grid, RealT* output,RealT* erg, int xsize, int ysize);
__device__ inline REAL dev_skalarprod_spinor(dev_spinor * s1, dev_spinor * s2);
__device__ inline REAL dev_squarenorm_spinor(dev_spinor * s1);
__device__ inline REAL dev_squarenorm_spinor_tex(int pos);
__global__ void dev_skalarprod_spinor_field2(dev_spinor* s1, dev_spinor* s2, REAL* erg);
__global__ void dev_squarenorm_spinor_field(dev_spinor* s1, REAL* erg);
__global__ void dev_skalarprod_spinor_field(dev_spinor* s1, dev_spinor* s2, REAL* erg);
template<class RealT>__global__ void dev_zero_spinor_field(typename dev_spinorT<RealT>::type* s1);
//__global__ void dev_copy_spinor_field(dev_spinor* s1, dev_spinor* s2);
template<class RealT1,class RealT2>__global__ void dev_copy_spinor_field(dev_spinorM(RealT1)* s1, dev_spinorM(RealT2)* s2);
//__global__ void dev_skalarmult_add_assign_spinor_field(dev_spinor* s1, REAL lambda, dev_spinor* s2, dev_spinor* so);
template<class RealT>__global__ void dev_skalarmult_add_assign_spinor_field(dev_spinorM(RealT)* s1, RealT lambda, dev_spinorM(RealT)* s2, dev_spinorM(RealT)* so);
//__global__ void dev_skalarmult_spinor_field(dev_spinor* s1, REAL lambda, dev_spinor* so);
template<class RealT>__global__ void dev_skalarmult_spinor_field(dev_spinorM(RealT)* s1, RealT lambda, dev_spinorM(RealT)* so);
//__global__ void dev_complexmult_spinor_field(dev_spinor* s1, dev_complex lambda, dev_spinor* so);
template<class RealT>__global__ void dev_complexmult_spinor_field(dev_spinorM(RealT)* s1, dev_complexM(RealT) lambda, dev_spinorM(RealT)* so);
__global__ void he_cg_init (int* grid, REAL param_kappa, REAL param_mu, dev_complex k0, dev_complex k1, dev_complex k2, dev_complex k3);
extern "C" int find_devices();
extern "C" int bind_texture_spin(dev_spinor* s, int i);
extern "C" int unbind_texture_spin(int i);
extern "C" int bind_texture_gf(dev_su3_2v * gf);
extern "C" int unbind_texture_gf();
extern "C" int bind_texture_nn(int* nn);
extern "C" int unbind_texture_nn();
extern "C" void test_operator(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize);
//extern "C" int dev_cg(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, int rescalekappa);
template<class RealT,template<class MixedsolveOperatorSRealT>class MixedsolveOperatorT>int dev_cg   (dev_su3_2vM(RealT)* gf, dev_spinorM(RealT)* spinin, dev_spinorM(RealT)* spinout, dev_spinorM(RealT)* spin0, dev_spinorM(RealT)* spin1, dev_spinorM(RealT)* spin2, dev_spinorM(RealT)* spin3, dev_spinorM(RealT)* spin4, int* grid, int* nn_grid, int rescalekappa,MixedsolveOperatorT<RealT>& mixedsolveOperator, REALD initial_sourcesquarenorm, bool rel_prec, double finalEps/*, bool& reachedFinalPrecision*/);
//extern "C" int dev_cg_eo(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, REAL epsfinal);
template<class RealT                                                                  >int dev_cg_eo(dev_su3_2vM(RealT)* gf, dev_spinorM(RealT)* spinin, dev_spinorM(RealT)* spinout, dev_spinorM(RealT)* spin0, dev_spinorM(RealT)* spin1, dev_spinorM(RealT)* spin2, dev_spinorM(RealT)* spin3, dev_spinorM(RealT)* spin4, int* grid, int* nn_grid, RealT epsfinal, MixedsolveParameter<RealT>& mixedsolveParameter);
void initnn();
void initnn_eo();
void shownn_eo();
void show_su3(su3 gf1);
void show_dev_su3(dev_su3 gf1);
void lptovec(int k);
void shownn();
//void su3to2vf4(su3** gf, dev_su3_2v* h2d_gf);
template<class RealT> void su3to2vf4(su3** gf, typename dev_su3_2vT<RealT>::type* h2d_gf);
//void su3to8(su3** gf, dev_su3_8* h2d_gf);
template<class RealT> void su3to8(su3** gf, typename dev_su3_8T<RealT>::type* h2d_gf);
void reconstructgf_2v (dev_su3* gf);
template<class RealT>__global__ void dev_check_gauge_reconstruction_8(typename dev_su3_2vT<RealT>::type* gf, int pos, typename dev_su3T<RealT>::type * outgf1, typename dev_su3T<RealT>::type* outgf2);
//void check_gauge_reconstruction_8(su3 ** gf1, dev_su3_2v * gf2, int ind1, int mu);
template<class RealT>void check_gauge_reconstruction_8(su3 ** gf1, dev_su3_2vM(RealT) * gf2, int ind1, int mu, MixedsolveParameter<RealT>& mixedsolveParameter);
void reconstructgf_8 (dev_su3_8 * h2d_gf, dev_su3* gf);
//void showcompare_gf(int t, int x, int y, int z, int mu);
template<class RealT>void showcompare_gf(int t, int x, int y, int z, int mu, MixedsolveParameter<RealT>& mixedsolveParameter);
//void convert2double_spin(dev_spinor* spin, spinor* h2d);
template<class RealT> void convert2double_spin(typename dev_spinorT<RealT>::type* spin, spinor* h2d);
//void convert2REAL4_spin(spinor* spin, dev_spinor* h2d);
template<class RealT> void convert2REAL4_spin(spinor* spin, typename dev_spinorT<RealT>::type* h2d);
//void init_mixedsolve(su3** gf);
template<class RealT>MixedsolveParameter<RealT>* init_mixedsolve(su3** gf);
//void init_mixedsolve_eo(su3** gf);
template<class RealT>MixedsolveParameter<RealT>* init_mixedsolve_eo(su3** gf);
//void finalize_mixedsolve();
template<class RealT>void finalize_mixedsolve(MixedsolveParameter<RealT>* mixedsolveParameterP);
//extern "C" int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, double eps, const int rel_prec,const int N);
template<class RealT,template<class MixedsolveOperatorSRealT>class MixedsolveOperatorT>int mixed_solveT                                 (spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N, MixedsolveOperatorT<RealT>& mixedsolveOperator);
extern "C"                                                                             int mixed_solve                                  (spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N);
extern "C"                                                                             int mixed_solveD                                 (spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N);
extern "C"                                                                             int mixed_solve_DiracDaggerDirac                 (spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N);
extern "C"                                                                             int mixed_solve_DiracDaggerDiracD                (spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N);
extern "C"                                                                             int mixed_solve_DiracDaggerDiracDiracDaggerDirac (spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N);
extern "C"                                                                             int mixed_solve_DiracDaggerDiracDiracDaggerDiracD(spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N);

void dummy (dev_spinor* a, dev_spinor* b);
//void benchmark(spinor * const Q);
template<class RealT>void benchmark(spinor * const Q,MixedsolveParameter<RealT>& mixedsolveParameter);
//extern "C" int mixed_solve_eo (spinor * const P, spinor * const Q, const int max_iter, double eps, const int rel_prec, const int N);
template<class RealT>int mixed_solve_eoT (spinor * const P, spinor * const Q, const int max_iter, double eps, const int rel_prec, const int N);
extern "C"           int mixed_solve_eo  (spinor * const P, spinor * const Q, const int max_iter, double eps, const int rel_prec, const int N);
extern "C"           int mixed_solve_eoD (spinor * const P, spinor * const Q, const int max_iter, double eps, const int rel_prec, const int N);



#define _MIXED_SOLVE_H_

#endif




