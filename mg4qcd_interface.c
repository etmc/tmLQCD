/***********************************************************************
 *
 * Copyright (C) 2016 Simone Bacchio, Jacob Finkenrath
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
 * Interface for MG4QCD
 *
 *******************************************************************************/
#include "mg4qcd_interface.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "boundary.h"
#include "linalg/convert_eo_to_lexic.h"
#include "solver/solver.h"
#include "solver/solver_field.h"
#include "gettime.h"
#include "read_input.h"
#include "mg4qcd.h"
#include "linalg_eo.h"

#define MGTEST

MG4QCD_Init mg_init;
MG4QCD_Parameters mg_params;
MG4QCD_Status mg_status;
int mg_do_setup=1; //if one do or redo the setup
int mg_update_gauge=1; //TODO set to zero if gaugefield is up to date, set to one if it has to be updated
int mg_update_setup=0; //Number of additional setup iteration 
int mg_initialized=0;
int mg_setup_iter=5;
int mg_Nvec=24;
int mg_lvl=3;
int mg_blk[4] = {0, 0, 0, 0};
double mg_cmu_factor = 5.0;
double mg_dtau = 0.0;


static int Cart_rank(MPI_Comm comm, const int coords[], int *rank) {
   int coords_l[4];
   
   coords_l[0]=coords[0];
   coords_l[1]=coords[3];
   coords_l[2]=coords[2];
   coords_l[3]=coords[1];
   
   return MPI_Cart_rank(comm, coords_l, rank);
}

static int Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) {
   int stat;
   stat=MPI_Cart_coords(comm, rank, maxdims, coords);
   int tmp=coords[1];
   coords[1]=coords[3];
   coords[3]=tmp;
   return stat;
}

static int conf_index_fct(int t, int z, int y, int x, int mu)
{
   int id;
   
   id=(g_ipt[t][x][y][z])*72; //9*2*4
   id+=18*((mu%2==0)?mu:((mu==1)?3:1));//9*2
   
   return id;
}

static int vector_index_fct(int t, int z, int y, int x )
{
   int id;
   
   id=24*(g_ipt[t][x][y][z]);
   
   return id;
}

void MG_init()
{

   mg_init.comm_cart=g_cart_grid;
   mg_init.Cart_rank=Cart_rank;
   mg_init.Cart_coords=Cart_coords;
   
   mg_init.global_lattice[0]=T*N_PROC_T;
   mg_init.global_lattice[1]=LZ*N_PROC_Z;
   mg_init.global_lattice[2]=LY*N_PROC_Y;
   mg_init.global_lattice[3]=LX*N_PROC_X;
   
   mg_init.procs[0]=N_PROC_T;
   mg_init.procs[1]=N_PROC_Z;
   mg_init.procs[2]=N_PROC_Y;
   mg_init.procs[3]=N_PROC_X;
   
   for(int i = 0; i<4; i++)
      if(mg_blk[i]==0)
	 mg_blk[i]=(((L/N_PROC_X)%2==0)?(((L/N_PROC_X)%4==0)?4:2):
				        (((L/N_PROC_X)%3==0)?3:1));
   
   mg_init.block_lattice[0]=mg_blk[0];
   mg_init.block_lattice[1]=mg_blk[1];
   mg_init.block_lattice[2]=mg_blk[2];
   mg_init.block_lattice[3]=mg_blk[3];
   
   if (X0==0)
      mg_init.bc=1;
   else
      mg_init.bc=3;
   
   mg_init.number_of_levels=mg_lvl;
#ifdef OMP
   mg_init.number_openmp_threads=omp_num_threads;
#else
   mg_init.number_openmp_threads=1;
#endif   
   mg_init.kappa=g_kappa;
   mg_init.mu=0.5*g_mu/g_kappa;
   
   if (g_c_sw==-1.00)
      mg_init.csw=0.0;
   else
      mg_init.csw=g_c_sw;
   
   MG4QCD_init( &mg_init, &mg_params, &mg_status);
   
   if (mg_status.success!=mg_lvl)
   {
      if (g_proc_id == 0)
      {
	 printf("MG WARNING: %d level initialized instead of %d\n",mg_status.success,mg_lvl);
         printf("MG WARNING: parameter: mg_lvl is changed to %d\n\n",mg_status.success);
      }
      mg_lvl=mg_status.success;
   }
   
   mg_params.conf_index_fct=conf_index_fct;
   mg_params.vector_index_fct=vector_index_fct;
   mg_params.print=1; // TODO
   /*** in MG4QCD
   *** Printing level:
   ***  -1: silent (errors or warnings)
   ***   0: minimal //default
   ***   1: with timings and iterations
   ***   2: verbose
   ***/
   
   mg_params.coarse_mu[mg_params.number_of_levels-1]=mg_cmu_factor*mg_params.mu; // input param mg_cmu_factor
   mg_params.mg_basis_vectors[0]=mg_Nvec;
   for (int j=1;j < mg_params.number_of_levels-1; j++)
      mg_params.mg_basis_vectors[j]=fmax(28,mg_params.mg_basis_vectors[j-1]);
   
   mg_params.setup_iterations[0]=mg_setup_iter;
  
   MG4QCD_update_parameters(&mg_params, &mg_status);
   
}

void MG_update_mu(double mu, double rho)
{
   // mu_oo = mu+rho -> Update this here
   if (mu!=mg_params.mu)
   {
      mg_params.mu = mu;
      for (int j=0;j < mg_params.number_of_levels-1; j++)
	 mg_params.coarse_mu[j] = mu;
      mg_params.coarse_mu[mg_params.number_of_levels-1] =  mg_cmu_factor*mu;
      MG4QCD_update_parameters(&mg_params, &mg_status);
   }
}

void MG_finalize()
{
   MG4QCD_finalize();
}

static inline void print_out_vector(  spinor * v)
{
   if (g_proc_id == 0)
      for (int j=0;j<12;j++)
	 printf("%d  %lf %lf\n",j, *(((double*) v) +2*j), *(((double*) v) +2*j+1) );
   
}



int MG_solver_degenerate(spinor * const phi_new, spinor * const phi_old,
		   const double precision, const int max_iter,
                   const int solver_flag, const int rel_prec,
                   const int even_odd_flag, su3 **gf, matrix_mult f)
{
   double mu;
   spinor ** solver_field = NULL;

   if (mg_initialized==0) {
      MG_init();
      mg_initialized = 1;
      if (g_proc_id == 0)
	 printf("MG4QCD initialized\n");
      MPI_Barrier(MPI_COMM_WORLD);
   }
   
   if (mg_update_gauge==1) {
      MG4QCD_set_configuration( (double*) &(gf[0][0]), &mg_status );
      mg_update_gauge = 0;
      if (mg_status.success == 1)
	 if (g_proc_id == 0) 
	    printf("MG4QCD cnfg set, plaquette %e\n", mg_status.info);
   }
   
   if (mg_do_setup==1) {
      if (g_proc_id == 0)
	 printf("MG4QCD running setup\n");
      MG4QCD_setup(&mg_status);
      mg_do_setup = 0;
      if (mg_status.success == 1)
	 if (g_proc_id == 0)	
	    printf("MG4QCD setup ran, time %.2f sec (%.2f %% on coarse grid)\n",
		     mg_status.time, 100.*(mg_status.coarse_time/mg_status.time));
   }
   
   if (mg_update_setup>0) {
      if (g_proc_id == 0)
	 printf("MG4QCD updating setup\n");
      MG4QCD_update_setup(mg_update_setup, &mg_status);
      mg_update_setup = 0;
      if (mg_status.success == 1)
	 if (g_proc_id == 0)	
	    printf("MG4QCD setup ran, time %.2f sec (%.2f %% on coarse grid)\n",
		     mg_status.time, 100.*(mg_status.coarse_time/mg_status.time));
   }
  
#ifndef MGTEST   
   init_solver_field(&solver_field, VOLUMEPLUSRAND, 3);
#else
   double differ[2];
   init_solver_field(&solver_field, VOLUMEPLUSRAND, 4);
#endif
   
  // for rescaling  convention in MG4QCD: (4+m)*\delta_{x,y} in tmLQCD: 1*\delta_{x,y} -> rescale by 1/4+m
   double mg_scale=0.25/(g_kappa*g_kappa);
   
   if (even_odd_flag==1)
   {
      convert_odd_to_lexic(solver_field[0], phi_old);
      mul_gamma5(solver_field[0],VOLUME);
      MG4QCD_solve( (double*) solver_field[1], (double*) solver_field[0], sqrt(precision), &mg_status );
      //MG4QCD_apply_operator( (double*) solver_field[1], (double*) solver_field[0], &mg_status );
      MG_update_mu(-mu, 0.0);
      // project Odd --> if even - odd
      set_even_to_zero(solver_field[1]);
      mul_gamma5(solver_field[1],VOLUME);
      MG4QCD_solve( (double*) solver_field[2], (double*) solver_field[1], sqrt(precision), &mg_status );
      mul_r(solver_field[2],mg_scale,solver_field[2],VOLUME);
      convert_lexic_to_odd(phi_new, solver_field[2]);
   }
   else
   {
      mul_gamma5(solver_field[0],VOLUME);
      MG4QCD_solve( (double*) solver_field[1], (double*) solver_field[0], sqrt(precision), &mg_status );
      //MG4QCD_apply_operator( (double*) solver_field[1], (double*) solver_field[0], &mg_status );
      MG_update_mu(-mu, 0.0);
      
      mul_gamma5(solver_field[1],VOLUME);
      MG4QCD_solve( (double*) solver_field[2], (double*) solver_field[1], sqrt(precision), &mg_status );
      //MG4QCD_apply_operator( (double*) solver_field[1], (double*) solver_field[0], &mg_status );
      mul_r(solver_field[2],mg_scale,solver_field[2],VOLUME);  
   }
   
   
   if (mg_status.success) {
      if (g_proc_id == 0)
	 printf("Solving time %.2f sec (%.1f %% on coarse grid)\n", mg_status.time,
	      100.*(mg_status.coarse_time/mg_status.time));
      if (g_proc_id == 0)
	 printf("Total iterations on fine grid %d\n", mg_status.iter_count);
      if (g_proc_id == 0)
	 printf("Total iterations on coarse grids %d\n", mg_status.coarse_iter_count);
    } 
   
#ifdef MGTEST
   double diff2[3]; 
   bicgstab_complex(solver_field[3], solver_field[0], max_iter, precision, rel_prec, VOLUME, f);
   
   diff2[0] = sqrt(square_norm(solver_field[0], VOLUME, 1));
   diff2[1] = sqrt(square_norm(solver_field[1], VOLUME, 1));
   diff2[3] = sqrt(square_norm(solver_field[3], VOLUME, 1));
   printf("%e %e %e\n",diff2[0],diff2[1],diff2[3]);
   diff(solver_field[3], solver_field[3], solver_field[1], VOLUME);
   
   differ[0] = sqrt(square_norm(solver_field[3], VOLUME, 1));
   differ[1] = sqrt(square_norm(solver_field[1], VOLUME, 1));
   if (g_proc_id == 0)
	 printf("Norm of the Difference of the Solution || D_{tmLQC}^{-1} s - D_{MG4QCD}^{-1}*s||/||s|| = %e/%e = %e \n", differ[0],differ[1],differ[0]/differ[1]);
   
   
   f(solver_field[3], solver_field[1]);
   diff(solver_field[1], solver_field[3], solver_field[0], VOLUME);
   differ[0] = sqrt(square_norm(solver_field[1], VOLUME, 1));
   differ[1] = sqrt(square_norm(solver_field[0], VOLUME, 1));
   if (g_proc_id == 0)
	 printf("Norm of the rel. Residual || s -D_{tmLQC} *D_{MG4QCD}^{-1}*s||/||s|| = %e/%e = %e \n", differ[0],differ[1],differ[0]/differ[1]);
   
#endif
   
#ifndef MGTEST   
   finalize_solver(solver_field, 2);
#else
   finalize_solver(solver_field, 3);
#endif   
   
   return mg_status.iter_count;
}




int MG_solver(spinor * const Even_new, spinor * const Odd_new,
                   spinor * const Even, spinor * const Odd,
                   const double precision, const int max_iter,
                   const int solver_flag, const int rel_prec,
                   const int even_odd_flag, su3 **gf, matrix_mult f)
{
   
   spinor ** solver_field = NULL;
#ifndef MGTEST   
   init_solver_field(&solver_field, VOLUMEPLUSRAND, 2);
#else
   double differ[2];
   init_solver_field(&solver_field, VOLUMEPLUSRAND, 3);
#endif
   
   
   if (mg_initialized==0) {
      MG_init();
      mg_initialized = 1;
      if (g_proc_id == 0)
	 printf("MG4QCD initialized\n");
      MPI_Barrier(MPI_COMM_WORLD);
   }
   
   if (mg_update_gauge==1) {
      MG4QCD_set_configuration( (double*) &(gf[0][0]), &mg_status );
      mg_update_gauge = 0;
      if (mg_status.success == 1)
	 if (g_proc_id == 0) 
	    printf("MG4QCD cnfg set, plaquette %e\n", mg_status.info);
   }
   
   if (mg_do_setup==1) {
      if (g_proc_id == 0)
	 printf("MG4QCD running setup\n");
      MG4QCD_setup(&mg_status);
      mg_do_setup = 0;
      if (mg_status.success == 1)
	 if (g_proc_id == 0)	
	    printf("MG4QCD setup ran, time %.2f sec (%.2f %% on coarse grid)\n",
		     mg_status.time, 100.*(mg_status.coarse_time/mg_status.time));
   }
   
   if (mg_update_setup>0) {
      if (g_proc_id == 0)
	 printf("MG4QCD updating setup\n");
      MG4QCD_update_setup(mg_update_setup, &mg_status);
      mg_update_setup = 0;
      if (mg_status.success == 1)
	 if (g_proc_id == 0)	
	    printf("MG4QCD setup ran, time %.2f sec (%.2f %% on coarse grid)\n",
		     mg_status.time, 100.*(mg_status.coarse_time/mg_status.time));
   }
  
  
  // for rescaling  convention in MG4QCD: (4+m)*\delta_{x,y} in tmLQCD: 1*\delta_{x,y} -> rescale by 1/4+m
   double mg_scale=0.5/g_kappa;
   
   
   convert_eo_to_lexic(solver_field[0],  Even, Odd);
   MG4QCD_solve( (double*) solver_field[1], (double*) solver_field[0], sqrt(precision), &mg_status );
   //MG4QCD_apply_operator( (double*) solver_field[1], (double*) solver_field[0], &mg_status );
   mul_r(solver_field[1],mg_scale,solver_field[1],VOLUME);
   convert_lexic_to_eo(Even_new, Odd_new, solver_field[1]);
   
   
   
   if (mg_status.success) {
      if (g_proc_id == 0)
	 printf("Solving time %.2f sec (%.1f %% on coarse grid)\n", mg_status.time,
	      100.*(mg_status.coarse_time/mg_status.time));
      if (g_proc_id == 0)
	 printf("Total iterations on fine grid %d\n", mg_status.iter_count);
      if (g_proc_id == 0)
	 printf("Total iterations on coarse grids %d\n", mg_status.coarse_iter_count);
    } 
   
#ifdef MGTEST
   double diff2[3]; 
   bicgstab_complex(solver_field[2], solver_field[0], max_iter, precision, rel_prec, VOLUME, f);
   
   diff2[0] = sqrt(square_norm(solver_field[0], VOLUME, 1));
   diff2[1] = sqrt(square_norm(solver_field[1], VOLUME, 1));
   diff2[2] = sqrt(square_norm(solver_field[2], VOLUME, 1));
   printf("%e %e %e\n",diff2[0],diff2[1],diff2[2]);
   diff(solver_field[2], solver_field[2], solver_field[1], VOLUME);
   
   differ[0] = sqrt(square_norm(solver_field[2], VOLUME, 1));
   differ[1] = sqrt(square_norm(solver_field[1], VOLUME, 1));
   if (g_proc_id == 0)
	 printf("Norm of the Difference of the Solution || D_{tmLQC}^{-1} s - D_{MG4QCD}^{-1}*s||/||s|| = %e/%e = %e \n", differ[0],differ[1],differ[0]/differ[1]);
   
   
   f(solver_field[2], solver_field[1]);
   diff(solver_field[1], solver_field[2], solver_field[0], VOLUME);
   differ[0] = sqrt(square_norm(solver_field[1], VOLUME, 1));
   differ[1] = sqrt(square_norm(solver_field[0], VOLUME, 1));
   if (g_proc_id == 0)
	 printf("Norm of the rel. Residual || s -D_{tmLQC} *D_{MG4QCD}^{-1}*s||/||s|| = %e/%e = %e \n", differ[0],differ[1],differ[0]/differ[1]);
   
#endif
   
#ifndef MGTEST   
   finalize_solver(solver_field, 2);
#else
   finalize_solver(solver_field, 3);
#endif   
   
   return mg_status.iter_count;
}