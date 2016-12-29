/*********************************************************************** 
 * Interface for solving the eigenvalue problem A*x=lambda*x 
 * for a complex matrix A using ARPACK and PARPACK. The matrix
 * is accessed through a matrix-vector multiplication function.
 *
 * Author: A.M. Abdel-Rehim, 2014
 *
 * For reference see the driver programs zndrv1 and in the EXAMPLES 
 * subdriectories of ARPACK and PARPACK.
 *
 * This file is part of tmLQCD software suite
 ***********************************************************************/
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#ifdef TM_USE_MPI
# include <mpi.h>
#endif
#include "global.h"
#include "linalg_eo.h"
#include "start.h"
#include "quicksort.h"
#include "linalg/blas.h"
#include "linalg/lapack.h"
#include "solver/eigenvalues_arpack.h"
void evals_arpack(
  int n, 
  int nev, 
  int ncv, 
  int which,
  int use_acc,
  int cheb_k,
  double amin,
  double amax,
  _Complex double *evals, 
  _Complex double *v,
  double tol, 
  int maxiter, 
  matrix_mult av, 
  int *info, 
  int *nconv,
  char *arpack_logfile)
/*
  compute nev eigenvectors using ARPACK and PARPACK
  n     : (IN) size of the local lattice
  nev   : (IN) number of eigenvectors requested.
  ncv   : (IN) size of the subspace used to compute eigenvectors (nev+1) =< ncv < 12*n
          where 12n is the size of the matrix under consideration
  which : (IN) which eigenvectors to compute. Choices are:
          0: smallest real part "SR"
          1: largest real part "LR"
          2: smallest absolute value "SM"
          3: largest absolute value "LM"
          4: smallest imaginary part "SI"
          5: largest imaginary part "LI"
  use_acc: (IN) specify the polynomial acceleration mode
                0 no acceleration
                1 use acceleration by computing the eigenvectors of a shifted-normalized chebyshev polynomial
  cheb_k   : (IN) degree of the chebyshev polynomial to be used for acceleration (irrelevant when use_acc=0 and init_resid_arpack=0)
  amin,amax: (IN) bounds of the interval [amin,amax] for the acceleration polynomial (irrelevant when use_acc=0 and init_resid_arpack=0)
  evals : (OUT) array of size nev+1 which has the computed nev Ritz values
  v     : computed eigenvectors. Size is n*ncv spinors.
  tol    : Requested tolerance for the accuracy of the computed eigenvectors.
           A value of 0 means machine precision.
  maxiter: maximum number of restarts (iterations) allowed to be used by ARPACK
  av     : operator for computing the action of the matrix on the vector
           av(vout,vin) where vout is output spinor and vin is input spinors.
  info   : output from arpack. 0 means that it converged to the desired tolerance. 
           otherwise, an error message is printed to stderr 
  nconv  : actual number of converged eigenvectors.
  arpack_logfile: name of the logfile to be used by arpack
*/ 
{

  //print the input:
  //================
   if(g_proc_id == g_stdio_proc)
   {       
      fprintf(stdout,"# [eigenvalues_arpack] Input to eigenvalues_arpack\n");
      fprintf(stdout,"# [eigenvalues_arpack] ===========================\n");
      fprintf(stdout,"# [eigenvalues_arpack] n= %d\n", n);
      fprintf(stdout,"# [eigenvalues_arpack] The number of Ritz values requested is %d\n", nev);
      fprintf(stdout,"# [eigenvalues_arpack] The number of Arnoldi vectors generated is %d\n", ncv);
      fprintf(stdout,"# [eigenvalues_arpack] What portion of the spectrum which: %d\n", which);
      fprintf(stdout,"# [eigenvalues_arpack] polynomial acceleartion option %d\n", use_acc );
      fprintf(stdout,"# [eigenvalues_arpack] chebyshev polynomial paramaters: degree %d amin %+e amx %+e\n",cheb_k,amin,amax); 
      fprintf(stdout,"# [eigenvalues_arpack] The convergence criterion is %+e\n", tol);
      fprintf(stdout,"# [eigenvalues_arpack] maximum number of iterations for arpack %d\n",maxiter);
   }

   //create the MPI communicator
#ifdef TM_USE_MPI
   MPI_Fint mpi_comm_f = MPI_Comm_c2f(g_cart_grid);
   //MPI_Comm comm; //communicator used when we call PARPACK
   //int comm_err = MPI_Comm_dup(MPI_COMM_WORLD,&comm); //duplicate the MPI_COMM_WORLD to create a communicator to be used with arpack
   //if(comm_err != MPI_SUCCESS) { //error when trying to duplicate the communicator
   //  if(g_proc_id == g_stdio_proc){
   //    fprintf(stderr,"[eigenvalues_arpack] MPI_Comm_dup return with an error. Exciting...\n");
   //    exit(-1);
   //  }
   //}
#endif

   int parallel;
#ifdef TM_USE_MPI
     parallel=1;
#else
     parallel=0;
#endif

   int ido=0;           //control of the action taken by reverse communications
                        //set initially to zero

   char *bmat=strdup("I");     /* Specifies that the right hand side matrix
                                  should be the identity matrix; this makes
                                  the problem a standard eigenvalue problem.
                               */

   //matrix dimensions 
   int ldv,N,LDV;
  
   if(n==VOLUME) //full 
     ldv= VOLUMEPLUSRAND;
   else          //even-odd
     ldv= VOLUMEPLUSRAND/2;

   //dimesnions as complex variables
   N  =12*n;       //dimension
   LDV=12*ldv;   //leading dimension (including communication buffers)
   

   char *which_evals;

   if(which==0)
     which_evals=strdup("SR");
   if(which==1)
     which_evals=strdup("LR");
   if(which==2)
     which_evals=strdup("SM");
   if(which==3)
     which_evals=strdup("LM");
   if(which==4)
     which_evals=strdup("SI");
   if(which==5)
     which_evals=strdup("LI");


    //check
    if (which_evals == NULL ||
            (  strcmp("SR", which_evals)
            && strcmp("LR", which_evals)
            && strcmp("SI", which_evals)
            && strcmp("LI", which_evals)
            && strcmp("SM", which_evals)
            && strcmp("LM", which_evals)))
    {
        if(g_proc_id == g_stdio_proc)
          {fprintf(stderr,"[eigenvalues_arpack] Error: invalid value for which_evals\n"); fflush(stderr); exit(1);}
    }

   //check input
   if(nev>=N){
       if(g_proc_id == g_stdio_proc)
          {fprintf(stderr,"[eigenvalues_arpack] number of eigenvalues requested should be less than the size of the matrix.\n"); fflush(stderr); exit(1);}
   }

   if(ncv < (nev+1)){
       if(g_proc_id == g_stdio_proc)
          {fprintf(stderr,"[eigenvalues_arpack] search subspace must be larger than the number of requested eiegnvalues.\n"); fflush(stderr); exit(1);}
   }

   _Complex double *resid  = (_Complex double *) malloc(N*sizeof(_Complex double));
   if(resid == NULL){
    if(g_proc_id == g_stdio_proc)
       { fprintf(stderr,"[eigenvalues_arpack] Error: not enough memory for resid in eigenvalues_arpack.\n"); fflush(stderr); exit(1);}
   }


   int *iparam = (int *) malloc(11*sizeof(int));
   if(iparam == NULL){
       if(g_proc_id == g_stdio_proc)
       { fprintf(stderr,"[eigenvalues_arpack] Error: not enough memory for iparam in eigenvalues_arpack.\n"); fflush(stderr); exit(1);}
   }


   iparam[0]=1;  //use exact shifts

   iparam[2]=maxiter;

   iparam[3]=1;

   iparam[6]=1;


   int *ipntr  = (int *) malloc(14*sizeof(int));

   _Complex double *workd  = (_Complex double *) malloc(3*N*sizeof(_Complex double)); 

   int lworkl=(3*ncv*ncv+5*ncv)*2; //just allocate more space

   _Complex double *workl=(_Complex double *) malloc(lworkl*sizeof(_Complex double));

   double *rwork  = (double *) malloc(ncv*sizeof(double));

   int rvec=1;       //always call the subroutine that computes orthonormal bais for the eigenvectors

   char *howmany=strdup("P");   //always compute orthonormal basis

   int *select = (int *) malloc(ncv*sizeof(int)); //since all Ritz vectors or Schur vectors are computed no need to initialize this array

   _Complex double sigma;
    
   _Complex double *workev = (_Complex double *) malloc(2*ncv*sizeof(_Complex double));

   double *sorted_evals = (double *) malloc(ncv*sizeof(double)); //will be used to sort the eigenvalues
   int *sorted_evals_index = (int *) malloc(ncv*sizeof(int)); 

   if((ipntr == NULL) || (workd==NULL) || (workl==NULL) || (rwork==NULL) || (select==NULL) || (workev==NULL) || (sorted_evals==NULL) || (sorted_evals_index==NULL)){
       if(g_proc_id == g_stdio_proc)
       { fprintf(stderr,"[eigenvalues_arpack] Error: not enough memory for ipntr,workd,workl,rwork,select,workev,sorted_evals,sorted_evals_index in eigenvalues_arpack.\n"); fflush(stderr); exit(1);}
   }

   double d1,d2,d3;

   
   (*info) = 0;                 //means use a random starting vector with Arnoldi

   void *_x,*_ax,*_r,*_tmps1,*_tmps2;   //spinors that might be needed
   spinor *x,*ax,*r,*tmps1,*tmps2;      //spinors that might be needed

   #if (defined SSE || defined SSE2 || defined SSE3)
   _x = malloc((ldv+ALIGN_BASE)*sizeof(spinor));
   if(_x==NULL){
       if(g_proc_id == g_stdio_proc)
       { fprintf(stderr,"[eigenvalues_arpack] Error: not enough memory for _x in eigenvalues_arpack.\n"); fflush(stderr); exit(1);}
   }
   else
      x  = (spinor *) ( ((unsigned long int)(_x)+ALIGN_BASE)&~ALIGN_BASE);


   _ax = malloc((ldv+ALIGN_BASE)*sizeof(spinor));
   if(_ax==NULL){
       if(g_proc_id == g_stdio_proc)
       { fprintf(stderr,"[eigenvalues_arpack] Error: not enough memory for _ax in eigenvalues_arpack.\n"); fflush(stderr); exit(1);}
   }
   else
     ax  = (spinor *) ( ((unsigned long int)(_ax)+ALIGN_BASE)&~ALIGN_BASE);


   _tmps1 = malloc((ldv+ALIGN_BASE)*sizeof(spinor));
   if(_tmps1==NULL){
       if(g_proc_id == g_stdio_proc)
       { fprintf(stderr,"[eigenvalues_arpack] Error: not enough memory for _tmps1 in eigenvalues_arpack.\n"); fflush(stderr); exit(1);}
   }
   else
     tmps1  = (spinor *) ( ((unsigned long int)(_tmps1)+ALIGN_BASE)&~ALIGN_BASE);


   _tmps2 = malloc((ldv+ALIGN_BASE)*sizeof(spinor));
   if(_tmps2==NULL){
       if(g_proc_id == g_stdio_proc)
       { fprintf(stderr,"[eigenvalues_arpack] Error: not enough memory for _tmps2 in eigenvalues_arpack.\n"); fflush(stderr); exit(1);}
   }
   else
     tmps2  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);


   #else

   x  = (spinor *) calloc(ldv,sizeof(spinor));
   ax = (spinor *) calloc(ldv,sizeof(spinor));
   tmps1  = (spinor *) calloc(ldv,sizeof(spinor));
   tmps2 = (spinor *) calloc(ldv,sizeof(spinor));
   if((x==NULL) || (ax==NULL) || (tmps1==NULL) || (tmps2==NULL)){
       if(g_proc_id == g_stdio_proc)
       { fprintf(stderr,"[eigenvalues_arpack] Error: not enough memory in x, ax, tmps1, tmps2 in eigenvalues_arpack.\n"); fflush(stderr); exit(1);}
   }

   #endif

   int i,j;

   /* Code added to print the log of ARPACK */
   //const char *arpack_logfile;
   int arpack_log_u = 9999;

#ifndef TM_USE_MPI
   //sprintf(arpack_logfile,"ARPACK_output.log");
   if ( NULL != arpack_logfile ) {
     /* correctness of this code depends on alignment in Fortran and C 
	being the same ; if you observe crashes, disable this part */
     _AFT(initlog)(&arpack_log_u, arpack_logfile, strlen(arpack_logfile));
     int msglvl0 = 0,
       msglvl1 = 1,
       msglvl2 = 2,
       msglvl3 = 3;
     _AFT(mcinitdebug)(
		  &arpack_log_u,      /*logfil*/
		  &msglvl3,           /*mcaupd*/
		  &msglvl3,           /*mcaup2*/
		  &msglvl0,           /*mcaitr*/
		  &msglvl3,           /*mceigh*/
		  &msglvl0,           /*mcapps*/
		  &msglvl0,           /*mcgets*/
		  &msglvl3            /*mceupd*/);
     
     fprintf(stdout,"# [arpack] *** ARPACK verbosity set to mcaup2=3 mcaupd=3 mceupd=3; \n"
	     "# [arpack] *** output is directed to '%s';\n"
	     "# [arpack] *** if you don't see output, your memory may be corrupted\n",
	     arpack_logfile);
  }
#else
   //if( g_proc_id == g_stdio_proc ){
   //  sprintf(arpack_logfile,"ARPACK_output.log");
   //}
   if ( NULL != arpack_logfile 
	&& (g_proc_id == g_stdio_proc) ) {
     /* correctness of this code depends on alignment in Fortran and C 
	being the same ; if you observe crashes, disable this part */
     _AFT(initlog)(&arpack_log_u, arpack_logfile, strlen(arpack_logfile));
     int msglvl0 = 0,
       msglvl1 = 1,
       msglvl2 = 2,
       msglvl3 = 3;
     _AFT(pmcinitdebug)(
		   &arpack_log_u,      /*logfil*/
		   &msglvl3,           /*mcaupd*/
		   &msglvl3,           /*mcaup2*/
		   &msglvl0,           /*mcaitr*/
		   &msglvl3,           /*mceigh*/
		   &msglvl0,           /*mcapps*/
		   &msglvl0,           /*mcgets*/
		   &msglvl3            /*mceupd*/);
    
     fprintf(stdout,"# [arpack] *** ARPACK verbosity set to mcaup2=3 mcaupd=3 mceupd=3; \n"
	    "# [arpack] *** output is directed to '%s';\n"
	    "# [arpack] *** if you don't see output, your memory may be corrupted\n",
	    arpack_logfile);
   }
#endif   




   /*
     M A I N   L O O P (Reverse communication)  
   */

   do
   {

#ifndef TM_USE_MPI 
      _AFT(znaupd)(&ido,"I", &N, which_evals, &nev, &tol,resid, &ncv,
                   v, &N, iparam, ipntr, workd, 
                  workl, &lworkl,rwork,info,1,2);
#else
      _AFT(pznaupd)(&mpi_comm_f, &ido,"I", &N, which_evals, &nev, &tol, resid, &ncv,
                   v, &N, iparam, ipntr, workd, 
                   workl, &lworkl,rwork,info,1,2);
#endif

      if (ido == 99 || (*info) == 1)
            break;

      if ((ido==-1)||(ido==1)){

         assign_complex_to_spinor(x,workd+ipntr[0]-1,N);
         if((use_acc==0))
           av(ax,x);
         else 
           cheb_poly_op(ax,x,av,n,amin,amax,cheb_k,tmps1,tmps2);

         assign_spinor_to_complex(workd+ipntr[1]-1, ax, n);
      }
   } while (ido != 99);
   
/*
 Check for convergence 
*/
     if ( (*info) < 0 ) 
     {
         if(g_proc_id == g_stdio_proc){
            fprintf(stderr,"[eigenvalues_arpack] Error with _naupd, info = %d\n", *info);
            fprintf(stderr,"[eigenvalues_arpack] Check the documentation of _naupd\n");}
     }
     else 
     {
        (*nconv) = iparam[4];
        if(g_proc_id == g_stdio_proc){
          fprintf(stderr,"[eigenvalues_arpack] number of converged eigenvectors = %d\n", *nconv);}

        //compute eigenvectors 
#ifndef TM_USE_MPI
        _AFT(zneupd) (&rvec,"P", select,evals,v,&N,&sigma, 
                     workev,"I",&N,which_evals,&nev,&tol,resid,&ncv, 
                     v,&N,iparam,ipntr,workd,workl,&lworkl, 
                     rwork,info,1,1,2);
#else
        _AFT(pzneupd) (&mpi_comm_f,&rvec,"P", select,evals, v,&N,&sigma, 
                       workev,"I",&N,which_evals,&nev,&tol, resid,&ncv, 
                       v,&N,iparam,ipntr,workd,workl,&lworkl, 
                       rwork,info,1,1,2);
#endif


        if( (*info)!=0) 
        {
           if(g_proc_id == g_stdio_proc){
             fprintf(stderr,"[eigenvalues_arpack] Error with _neupd, info = %d \n",(*info));
             fprintf(stderr,"[eigenvalues_arpack] Check the documentation of _neupd. \n");}
        }
        else //report eiegnvalues and their residuals
        {
             if(g_proc_id == g_stdio_proc){
               fprintf(stdout,"[eigenvalues_arpack] Ritz Values and their errors\n");
               fprintf(stdout,"[eigenvalues_arpack] ============================\n");
             }

             (*nconv) = iparam[4];
             for(j=0; j< (*nconv); j++)
             {
               /* print out the computed ritz values and their error estimates */
               if(g_proc_id == g_stdio_proc)
                  fprintf(stdout,"# [eigenvalues_arpack] RitzValue[%06d]  %+e  %+e  error= %+e \n",j,creal(evals[j]),cimag(evals[j]),cabs(*(workl+ipntr[10]-1+j)));
               sorted_evals_index[j] = j;
               sorted_evals[j] = cabs(evals[j]);
             }

             //SORT THE EIGENVALUES in ascending order based on their absolute value
             quicksort((*nconv),sorted_evals,sorted_evals_index);
             //Print sorted evals
             if(g_proc_id == g_stdio_proc)
                fprintf(stdout,"# [eigenvalues_arpack] Sorted eigenvalues based on their absolute values\n");

             for(j=0; j< (*nconv); j++)
             {
               /* print out the computed ritz values and their error estimates */
               if(g_proc_id == g_stdio_proc)
                  fprintf(stdout,"# [eigenvalues_arpack] RitzValue[%06d]  %+e  %+e  error= %+e \n",j,creal(evals[sorted_evals_index[j]]),cimag(evals[sorted_evals_index[j]]),cabs(*(workl+ipntr[10]-1+sorted_evals_index[j])));

             }

        }

        /*Print additional convergence information.*/
        if( (*info)==1)
        {
           if(g_proc_id == g_stdio_proc)
             fprintf(stderr,"[eigenvalues_arpack] Maximum number of iterations reached.\n");
        }
        else
        { 
          
           if(g_proc_id == g_stdio_proc)
           {
              if((*info)==3)
              {  
                 fprintf(stderr,"[eigenvalues_arpack] No shifts could be applied during implicit\n");
                 fprintf(stderr,"[eigenvalues_arpack] Arnoldi update, try increasing NCV.\n");
              }
         
              fprintf(stdout,"# [eigenvalues_arpack] _NDRV1\n");
              fprintf(stdout,"# [eigenvalues_arpack] =======\n");
              fprintf(stdout,"# [eigenvalues_arpack] Size of the matrix is %d\n", N);
              fprintf(stdout,"# [eigenvalues_arpack] The number of Ritz values requested is %d\n", nev);
              fprintf(stdout,"# [eigenvalues_arpack] The number of Arnoldi vectors generated is %d\n", ncv);
              fprintf(stdout,"# [eigenvalues_arpack] What portion of the spectrum: %s\n", which_evals);
              fprintf(stdout,"# [eigenvalues_arpack] The number of converged Ritz values is %d\n", (*nconv) ); 
              fprintf(stdout,"# [eigenvalues_arpack] The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
              fprintf(stdout,"# [eigenvalues_arpack] The number of OP*x is %d\n", iparam[8]);
              fprintf(stdout,"# [eigenvalues_arpack] The convergence criterion is %+e\n", tol);
           }
          
        }
     }  //if(info < 0) else part


#ifndef TM_USE_MPI
     if (NULL != arpack_logfile)
       _AFT(finilog)(&arpack_log_u);
#else
     if(g_proc_id == g_stdio_proc){
       if (NULL != arpack_logfile){
	 _AFT(finilog)(&arpack_log_u);
       }
     }
#endif     


     //free memory
     free(resid);
     free(iparam);
     free(ipntr);
     free(workd);
     //free(zv);
     free(select);
     free(workev);
     #if ( (defined SSE) || (defined SSE2) || (defined SSE3) )
     free(_x); free(_ax); free(_tmps1); free(_tmps2);
     #else
     free(x); free(ax); free(tmps1); free(tmps2);
     #endif

     return;
}





