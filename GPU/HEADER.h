
#ifndef _MIXED_SOLVE_H_




///////////////////
// own functions //
///////////////////


// eo, nd 

void to_device (dev_spinor * device, spinor * host, dev_spinor * auxiliary, int size);

void to_host (spinor * host, dev_spinor * device, dev_spinor * auxiliary, int size);

__global__ void he_cg_init_nd_additional (double param_mubar, double param_epsbar);

__global__ void dev_mul_one_pm_imubar_gamma5 (dev_spinor * sin, dev_spinor * sout, float sign);

void init_mixedsolve_eo_nd(su3** gf, int use_eo);

void finalize_mixedsolve_eo_nd(void);
void update_gpu_gf(su3** gf);


void flopcount(unsigned long long int& total, int add);

extern "C" void benchmark_eo_nd (spinor * const Q_up, spinor * const Q_dn, int N, int use_eo);

int dev_cg_eo_nd (dev_su3_2v * gf,
              dev_spinor * P_up, dev_spinor * P_dn,
              dev_spinor * Q_up, dev_spinor * Q_dn,
	      float shift,
              int max_iter,
              int check_abs , int check_rel,
              double eps_abs, double eps_rel       );

extern "C" int mixedsolve_eo_nd (spinor * P_up, spinor * P_dn,
                                 spinor * Q_up, spinor * Q_dn, double shift,
                                 int max_iter, double eps_sq, int rel_prec, int use_eo, matrix_mult_nd f);

void set_global_sizes(int use_eo);

float cublasSdot_wrapper(int size, float * A, int incx, float * B, int incy);


// eo, nd, MPI

void convert2double_spin_mpi (dev_spinor* spin, spinor* h2d, int start, int end);
void convert2REAL4_spin_mpi (spinor* spin, dev_spinor* h2d, int start, int end);
void to_device_mpi (dev_spinor * device, spinor * host, dev_spinor * auxiliary, int size, int start, int end);
void to_host_mpi (spinor * host, dev_spinor * device, dev_spinor * auxiliary, int size, int start, int end);
void xchange_field_wrapper (dev_spinor * dev_spin, int ieo);
void xchange_field_wrapper_d (dev_spinor_d * dev_spin, int ieo);
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


int cg_eo_nd_mpi (dev_su3_2v * gf,
                  dev_spinor * P_up, dev_spinor * P_dn,
                  dev_spinor * Q_up, dev_spinor * Q_dn,
                  int max_iter,
                  int check_abs , int check_rel,
                  double eps_abs, double eps_rel       );

extern "C" int mixedsolve_eo_nd_mpi (spinor * P_up, spinor * P_dn,
                                     spinor * Q_up, spinor * Q_dn, double shift,
                                     int max_iter, double eps_sq, int rel_prec);

extern "C" int dev_cg_mms_tm_nd(spinor ** const Pup, spinor ** const Pdn, 
		 spinor * const Qup, spinor * const Qdn, 
		 solver_pm_t * solver_pm);





__global__ void dev_Hopping_Matrix_mpi (const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout,
                                        int * dev_iup, int * dev_idn, int * dev_eo2lexic, int * dev_lexic2eosub,
                                        int ieo);
__global__ void dev_Hopping_Matrix_ASYNC (const dev_su3_2v * gf, 
                                          const dev_spinor * sin, dev_spinor * sout,
                                          const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd,
                                          const int eo,
                                          int start, int size);

void HOPPING_ASYNC (dev_su3_2v * gf, 
                    dev_spinor * spinin, dev_spinor * spinout,
                    int * gfindex_site, int * gfindex_nextsite, int * nn_evenodd,
                    int ieo,
                    int gridsize, dim3 blocksize);
void HOPPING_ASYNC_D (dev_su3_2v_d * gf, 
                    dev_spinor_d * spinin, dev_spinor_d * spinout,
                    int * gfindex_site, int * gfindex_nextsite, int * nn_evenodd,
                    int ieo,
                    int gridsize, int blocksize);
void HOPPING_ASYNC_UPDN_D (dev_su3_2v_d * gf, 
                    dev_spinor_d * spinin_up, dev_spinor_d* spinin_dn, 
		    dev_spinor_d * spinout_up, dev_spinor_d* spinout_dn,
                    int * gfindex_site, int * gfindex_nextsite, int * nn_evenodd,
                    int ieo,
                    int gridsize, int blocksize);




void dev_Qtm_pm_ndpsi_mpi (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                                  dev_spinor * spinin_up , dev_spinor * spinin_dn , 
                                  int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                                  int gridsize3, int blocksize3, int gridsize4, int blocksize4);



void dev_Qtm_pm_ndpsi_mpi_ASYNC (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                                        dev_spinor * spinin_up , dev_spinor * spinin_dn ,
                                        int gridsize1, dim3 blocksize1, int gridsize2, int blocksize2,
                                        int gridsize3, int blocksize3, int gridsize4, int blocksize4);



void dev_Qtm_pm_psi(dev_spinor* spinin, dev_spinor* spinout);
void dev_Qtm_pm_psi_d(dev_spinor_d* spinin, dev_spinor_d* spinout);

void dev_Qtm_pm_ndpsi(dev_spinor * , dev_spinor * , dev_spinor * , dev_spinor * );
void dev_Qtm_pm_ndpsi_d (dev_spinor_d * spinout_up, dev_spinor_d * spinout_dn,
                         dev_spinor_d * spinin_up , dev_spinor_d * spinin_dn);

void dev_Q_pm_ndpsi(dev_spinor * const l_strange, dev_spinor * const l_charm, dev_spinor * const k_strange, dev_spinor * const k_charm);
void dev_Q_pm_ndpsi_d(dev_spinor_d * const l_strange, dev_spinor_d * const l_charm, dev_spinor_d * const k_strange, dev_spinor_d * const k_charm);

void dev_Qtm_minus_psi_d(dev_spinor_d* spinin, dev_spinor_d* spinout,dev_spinor_d* spin_eo1_d, int gridsize, int blocksize, int gridsize2, int blocksize2, int* dev_eoidx_even, int* dev_eoidx_odd,  int* dev_nn_eo, int* dev_nn_oe);				


__global__ void dev_tm_dirac_kappa(dev_su3_2v * gf, dev_spinor * sin, dev_spinor * sout, int * dev_nn);
void dev_tm_dirac_dagger_kappa(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, int *grid, int * nn_grid, float* output,float* erg, int xsize, int ysize);






/////////////////////////
// Florian's functions //
/////////////////////////

// eo, non-eo, non-nd

__device__ inline dev_complex dev_cconj (dev_complex c);
__device__ inline void dev_ccopy(dev_complex* von, dev_complex* nach);
__device__ inline float dev_cabssquare (dev_complex c);
__device__ inline float dev_cabsolute (dev_complex c);
__device__ inline  dev_complex dev_crealmult(dev_complex c1, float real);
__device__ inline dev_complex dev_cmult (dev_complex c1, dev_complex c2);
__device__ inline dev_complex dev_cadd (dev_complex c1, dev_complex c2);
__device__ inline dev_complex dev_cdiv(dev_complex c1, dev_complex c2);
__device__ inline dev_complex dev_csub(dev_complex c1, dev_complex c2);
__device__ inline dev_complex dev_initcomplex(float re, float im);
__device__ inline void dev_copy_spinor(dev_spinor *i1, dev_spinor *i2);
__device__ inline void dev_zero_spinor(dev_spinor *sin);
__device__ inline void dev_skalarmult_add_assign_spinor(dev_spinor *in, float lambda,dev_spinor * in2, dev_spinor * out);
__device__ inline void dev_complexmult_add_assign_spinor(dev_spinor * in, dev_complex lambda,dev_spinor * in2, dev_spinor * out);
__device__ inline void dev_complexcgmult_add_assign_spinor(dev_spinor * in, dev_complex lambda,dev_spinor * in2, dev_spinor * out);
__device__ void inline dev_skalarmult_spinor(dev_spinor * in, dev_complex lambda, dev_spinor * out);
__device__ void inline dev_skalarmult_gamma5_spinor(dev_spinor * out, dev_complex lambda, dev_spinor * in);
__device__ void inline dev_realmult_spinor(dev_spinor * in, float lambda);
__device__ void inline dev_realmult_spinor_assign(dev_spinor* out,  float lambda, dev_spinor* in);
__device__ void dev_assign_realmult_add_spinor(dev_spinor* out, float lambda, dev_spinor* in1,  dev_spinor* in2);
__device__ inline void dev_add_spinor_assign(dev_spinor * i1, dev_spinor * i2);
__device__ inline void dev_sub_spinor_assign(dev_spinor * i1, dev_spinor * i2);
__device__ void dev_su3MtV_spintex(dev_su3 M, int pos, dev_spinor * out);
__device__ void dev_su3MtV(dev_su3 M, const dev_spinor * s, dev_spinor * out);
__device__ void dev_su3MdaggertV(dev_su3 M, dev_spinor * s, dev_spinor * out);
__device__ void dev_Gamma0(dev_spinor * in);
__device__ void dev_Gamma3(dev_spinor * in);
__device__ void dev_Gamma2(dev_spinor * in);
__device__ void dev_Gamma1(dev_spinor * in);
__device__ void dev_Gamma5(dev_spinor * in);
__device__ void dev_Gamma5_assign(dev_spinor* out, dev_spinor* in);
__device__ void dev_GammatV(int mu, dev_spinor * in);
__device__ void dev_reconstructgf_2vtexref (const dev_su3_2v* field, int pos, dev_su3* gf);
__device__ void dev_reconstructgf_2vtexref_dagger (const dev_su3_2v* field, int pos, dev_su3* gf);
__device__ void dev_reconstructgf_8texref (const dev_su3_2v * field, int pos, dev_su3* gf);
__device__ void dev_reconstructgf_8texref_dagger (const dev_su3_2v* field,int pos, dev_su3* gf);
__global__ void dev_gamma5(dev_spinor * sin, dev_spinor * sout);
__global__ void dev_swapmu();
__global__ void dev_mul_one_pm_imu_inv(dev_spinor* sin, dev_spinor* sout, const float sign);
__global__ void dev_mul_one_pm_imu_sub_mul_gamma5(dev_spinor* sin1, dev_spinor* sin2, dev_spinor* sout, const float sign);
__device__ void dev_kappaP1_plus(dev_spinor * out, dev_spinor * in, float kappa);
__device__ void dev_kappaP1_minus(dev_spinor * out, dev_spinor * in, float kappa);
__device__ void dev_kappaP2_plus(dev_spinor * out, dev_spinor * in, float kappa);
__device__ void dev_kappaP2_minus(dev_spinor * out, dev_spinor * in, float kappa);
__device__ void dev_kappaP3_plus(dev_spinor * out, dev_spinor * in, float kappa);
__device__ void dev_kappaP3_minus(dev_spinor * out, dev_spinor * in, float kappa);
__device__ void dev_kappaP0_plus(dev_spinor * out, dev_spinor * in, dev_complex kappa);
__device__ void dev_kappaP0_minus(dev_spinor * out, dev_spinor * in, dev_complex kappa);
__global__ void dev_Hopping_Matrix(const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo);
__device__ inline float dev_skalarprod_spinor(dev_spinor * s1, dev_spinor * s2);
__device__ inline float dev_squarenorm_spinor(dev_spinor * s1);
__device__ inline float dev_squarenorm_spinor_tex(int pos);
__global__ void dev_skalarprod_spinor_field2(dev_spinor* s1, dev_spinor* s2, float* erg);
__global__ void dev_squarenorm_spinor_field(dev_spinor* s1, float* erg);
__global__ void dev_skalarprod_spinor_field(dev_spinor* s1, dev_spinor* s2, float* erg);
__global__ void dev_zero_spinor_field(dev_spinor* s1);
__global__ void dev_copy_spinor_field(dev_spinor* s1, dev_spinor* s2);
__global__ void dev_zero_spinor_field_d(dev_spinor_d* s1);
__global__ void dev_copy_spinor_field_d(dev_spinor_d* s1, dev_spinor_d* s2);
__global__ void dev_skalarmult_add_assign_spinor_field(dev_spinor* s1, float lambda, dev_spinor* s2, dev_spinor* so);
__global__ void dev_skalarmult_spinor_field(dev_spinor* s1, float lambda, dev_spinor* so);
__global__ void dev_complexmult_spinor_field(dev_spinor* s1, dev_complex lambda, dev_spinor* so);
__global__ void he_cg_init (int* grid, float param_kappa, float param_mu, dev_complex k0, dev_complex k1, dev_complex k2, dev_complex k3);
__global__ void dev_gather_rand(dev_spinor * sin, dev_spinor * rand, int start, int size);
__global__ void dev_spread_rand(dev_spinor * sin, dev_spinor * rand, int start, int size);


__global__ void dev_gather_rand_relup(dev_spinor * sin, dev_spinor * rand, int start, int size);
__global__ void dev_gather_rand_reldn(dev_spinor * sin, dev_spinor * rand, int start, int size);
__global__ void dev_spread_rand_relup(dev_spinor * sin, dev_spinor * rand, int start, int size);
__global__ void dev_spread_rand_reldn(dev_spinor * sin, dev_spinor * rand, int start, int size);



//custom blas kernels for float 
void init_blas(int vol);
void finalize_blas();
void set_gpu_work_layout(int eoflag);
__device__ inline void dev_read_spinor(dev_spinor *i1, dev_spinor *i2);
__global__ void dev_axpy (float alpha, dev_spinor* x, dev_spinor* y);
__global__ void dev_xpay (float alpha, dev_spinor* x, dev_spinor* y);
__global__ void dev_blasscal (float alpha, dev_spinor* y);
__global__ void dev_blascopy (dev_spinor* x, dev_spinor* y);
__global__ void dev_dot( float* redfield, dev_spinor* x,dev_spinor* y);
extern "C" float float_dotprod(dev_spinor* x, dev_spinor* y, int N);


extern "C" int find_devices();
extern "C" int bind_texture_spin(dev_spinor* s, int i);
extern "C" int bind_texture_spin_dn(dev_spinor* s, int i);
extern "C" int unbind_texture_spin(int i);
extern "C" int unbind_texture_spin_dn(int i);
extern "C" int bind_texture_gf(dev_su3_2v * gf);
extern "C" int unbind_texture_gf();
extern "C" int bind_texture_nn(int* nn);
extern "C" int unbind_texture_nn();
extern "C" void test_operator(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, float* output,float* erg, int xsize, int ysize);
int dev_cg(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, int rescalekappa);
int dev_cg_eo(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, float epsfinal);
void initnn();
void initnn_eo();
void shownn_eo();
void show_su3(su3 gf1);
void show_dev_su3(dev_su3 gf1);
void lptovec(int k);
void shownn();
void su3to2vf4(su3** gf, dev_su3_2v* h2d_gf);
void su3to8(su3** gf, dev_su3_8* h2d_gf);
void reconstructgf_2v (dev_su3* gf);
__global__ void dev_check_gauge_reconstruction_8(dev_su3_2v* gf, int pos, dev_su3 * outgf1, dev_su3* outgf2);
void check_gauge_reconstruction_8(su3 ** gf1, dev_su3_2v * gf2, int ind1, int mu);
void reconstructgf_8 (dev_su3_8 * h2d_gf, dev_su3* gf);
void showcompare_gf(int t, int x, int y, int z, int mu);
void convert2double_spin(dev_spinor* spin, spinor* h2d);
void convert2REAL4_spin(spinor* spin, dev_spinor* h2d);
extern "C" void init_mixedsolve(su3** gf);
extern "C" void init_mixedsolve_eo(su3** gf, int use_eo);
extern "C" void finalize_mixedsolve(int use_eo);
extern "C" int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, double eps, const int rel_prec,const int N);
void init_mixedsolve_fields(int eo);
void finalize_mixedsolve_fields();
void alloc_single_gf();
void dealloc_single_gf();
void check_mixedsolve_params();
void update_constants(int *grid);

void dummy (dev_spinor* a, dev_spinor* b);
void benchmark_eo(spinor * const Q);
void benchmark_eo_mpi(spinor * const Q);
extern "C" int mixed_solve_eo (spinor * const P, spinor * const Q, const int max_iter, double eps, const int rel_prec, const int N);


extern "C" double double_dotprod(dev_spinor_d* x, dev_spinor_d* y, int N);
extern "C" void update_constants_d(int *grid);
extern "C" void update_gpu_gf_d(su3** gf);



 
				 
extern "C" void gpu_deriv_Sb(const int ieo, spinor * const l, spinor * const k,
                               hamiltonian_field_t * const hf, const double factor);				 
extern "C" void gpu_gauge_derivative(int withrectangles, hamiltonian_field_t * const hf, double c_gauge, double c_rect);				 



#define _MIXED_SOLVE_H_

#endif




