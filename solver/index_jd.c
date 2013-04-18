/*******************************************************************************

 *
 * This routine computes the index with the Jacobi-Davidson method
 *
 * Author: Carsten Urbach, urbach@physik.fu-berlin.de
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include "global.h"
#include "start.h"
#include "sse.h"
#include "su3.h"
#include "linalg_eo.h"
#include "eigenvalues.h"
#include <io/eospinor.h>
#include "solver/solver.h"
#include "solver/jdher.h"
#include "solver/eigenvalues.h"
#include "operator/Dov_proj.h"
#include "gamma.h"
#include "index_jd.h"

#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>


double shift;

#define min(a,b)((a)<(b) ? (a) : (b))
#define max(a,b)((a)<(b) ? (b) : (a))

void index_jd(int * nr_of_eigenvalues_ov, 
	      const int max_iterations, const double precision_ov, char *conf_filename, 
	      const int nstore, const int method){
  
  _Complex double *eval;
  spinor  *eigenvectors_ov, *eigenvectors_ov_;
  spinor  *lowvectors, *lowvectors_;
  int i=0 , k=0, returncode=0, index = 0, determined = 0, signed_index = 0;
  char filename[120];
  FILE * ifs = NULL;
  matrix_mult Operator[2];
  double absdifference;
  const int N2 = VOLUMEPLUSRAND;

#ifdef MPI
  double atime, etime;
#endif
  double lowestmodes[20];
  int intsign, max_iter, first_blocksize = 1;
  int * idx = NULL;

  /**********************
   * For Jacobi-Davidson 
   **********************/
  int verbosity = 3, converged = 0, blocksize = 1, blockwise = 0;
  int solver_it_max = 50, j_max, j_min, v0dim = 0;
  double * eigenvalues_ov = NULL;
  double decay_min = 1.7, threshold_min = 1.e-3, prec;

  WRITER *writer=NULL;
  spinor *s;
  double sqnorm;
  paramsPropagatorFormat *propagatorFormat = NULL;
  
  double ap_eps_sq;
  int switch_on_adaptive_precision = 0;
  double ov_s = 0;

  /**********************                                                 
   * General variables                                                    
   **********************/

  eval= calloc((*nr_of_eigenvalues_ov),sizeof(_Complex double));
  shift = 0.0;

  //  ov_s = 0.5*(1./g_kappa - 8.) - 1.;
  ap_eps_sq = precision_ov*precision_ov; 

#if (defined SSE || defined SSE2 )
  eigenvectors_ov_= calloc(VOLUMEPLUSRAND*(*nr_of_eigenvalues_ov)+1, sizeof(spinor)); 
  eigenvectors_ov = (spinor *)(((unsigned long int)(eigenvectors_ov_)+ALIGN_BASE)&~ALIGN_BASE);
  lowvectors_ = calloc(2*first_blocksize*VOLUMEPLUSRAND+1, sizeof(spinor));
  lowvectors = (spinor *)(((unsigned long int)(lowvectors_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  //  eigenvectors_ov_ = calloc(VOLUMEPLUSRAND*(*nr_of_eigenvalues_ov), sizeof(spinor));
  eigenvectors_ov_ = calloc(VOLUMEPLUSRAND*(*nr_of_eigenvalues_ov), sizeof(spinor));
  lowvectors_ = calloc(2*first_blocksize*VOLUMEPLUSRAND, sizeof(spinor));
  eigenvectors_ov = eigenvectors_ov_;
  lowvectors = lowvectors_;
#endif

  //  idx = malloc((*nr_of_eigenvalues_ov)*sizeof(int));
  idx = malloc((*nr_of_eigenvalues_ov)*sizeof(int));
  Operator[0]=&Dov_proj_plus;
  Operator[1]=&Dov_proj_minus;
  
  if(g_proc_id == g_stdio_proc){
    printf("Computing first the two lowest modes in the positive and negative chirality sector, respectively\n");
    if(switch_on_adaptive_precision == 1) {
      printf("We have switched on adaptive precision with ap_eps_sq = %e!\n", ap_eps_sq);
    }
    printf("We have set the mass to zero within this computation!\n");
    fflush(stdout);
  }

  prec = precision_ov; 
  j_min = 8; j_max = 16;
  max_iter = 70;

#ifdef MPI
  atime = MPI_Wtime();
#endif

  v0dim = first_blocksize;
  blocksize = v0dim;
  for(intsign = 0; intsign < 2; intsign++){
    converged = 0;
    if(g_proc_id == g_stdio_proc){
      printf("%s chirality sector: \n", intsign ? "negative" : "positive");
      fflush(stdout);
    }
    if(max_iter == 70){
      /********************************************************************
       *
       * We need random start spinor fields, but they must be half zero,
       * that's why we apply the Projektor once
       *
       ********************************************************************/
      for(i = 0; i < first_blocksize; i++) {
	random_spinor_field(&lowvectors[(first_blocksize*intsign+i)*VOLUMEPLUSRAND],N2,0);
	Proj(&lowvectors[(first_blocksize*intsign+i)*VOLUMEPLUSRAND], 
	     &lowvectors[(first_blocksize*intsign+i)*VOLUMEPLUSRAND],N2, intsign);
      }
    }

    jdher(VOLUME*sizeof(spinor)/sizeof(_Complex double),
	  VOLUMEPLUSRAND*sizeof(spinor)/sizeof(_Complex double),
	  shift, prec, blocksize, j_max, j_min, 
	  max_iter, blocksize, blockwise, v0dim, (_Complex double*) &lowvectors[first_blocksize*intsign*VOLUMEPLUSRAND],
	  CG, solver_it_max,
	  threshold_min, decay_min, verbosity,
	  &converged, (_Complex double*) &lowvectors[first_blocksize*intsign*VOLUMEPLUSRAND], 
	  &lowestmodes[first_blocksize*intsign],
	  &returncode, JD_MINIMAL, 1,
	  Operator[intsign]);

    if(converged != blocksize && max_iter == 70){
      if(g_proc_id == g_stdio_proc){
	printf("Restarting %s chirality sector with more iterations!\n", intsign ? "negative" : "positive");
	fflush(stdout);
      }
      max_iter = 140;
      intsign-=1;
    }
    else {
      max_iter = 70;
      /* Save the allready computed eigenvectors_ov */
      for(i = 0; i< first_blocksize; i++) {
	sprintf(filename, "eigenvector_of_D%s.%.2d.%s.%.4d",((intsign==0)?"plus":"minus"),i , conf_filename, nstore);

	  construct_writer(&writer, filename, 0);
	  /* todo write propagator format */
	  propagatorFormat = construct_paramsPropagatorFormat(64, 1);
	  write_propagator_format(writer, propagatorFormat);
	  free(propagatorFormat);


	  s=(spinor*)&lowvectors[first_blocksize*intsign*VOLUMEPLUSRAND];
	  write_spinor(writer, &s,NULL, 1, 64);
	  destruct_writer(writer);
	  writer=NULL;
	  sqnorm=square_norm(s,VOLUME,1);
	  printf(" wrote eigenvector of overlap operator !!! | |^2 = %e \n",sqnorm);


      }
    }
  }

#ifdef MPI
  etime = MPI_Wtime();
  if(g_proc_id == g_stdio_proc){
    printf("It took %f sec to determine the sector with zero modes, if any!\n", etime-atime);
  }
#endif

  /*Compare the two lowest modes */
  absdifference = fabs(lowestmodes[0]-lowestmodes[first_blocksize]);
  if(absdifference < 0.1*max(lowestmodes[0],lowestmodes[first_blocksize])){
    /* They are equal within the errors */
    if(g_proc_id == g_stdio_proc){
      printf("Index is 0!\n");
      fflush(stdout);
      sprintf(filename, "eigenvalues_of_overlap_proj.%s.%.4d", conf_filename, nstore);
      ifs = fopen(filename, "w");  
      printf("\nThe following lowest modes have been computed:\n");
      fprintf(ifs, "Index is 0\n\n");
      fprintf(ifs, "Sector with positive chirality:\n");
      for(i = 0; i < first_blocksize; i++) {
	lowestmodes[i] = 2.*(1.+ov_s)*lowestmodes[i];
	fprintf(ifs, "%d %e positive\n", i, lowestmodes[i]);
	printf("%d %e positive\n", i, lowestmodes[i]);
      }
      fprintf(ifs, "Sector with negative chirality:\n");
      for(i = 0; i < first_blocksize; i++) {
	lowestmodes[i+first_blocksize] = 2.*(1.+ov_s)*lowestmodes[i+first_blocksize];
	fprintf(ifs, "%d %e negative\n", i, lowestmodes[i+first_blocksize]);
	printf("%d %e negative\n", i, lowestmodes[i+first_blocksize]);
      }
      fclose(ifs);
      for(k = 0; k < 2; k++) {
	sprintf(filename, "eigenvalues_of_D%s.%s.%.4d", 
		k ? "minus" : "plus", conf_filename, nstore);
	ifs = fopen(filename, "w");
	fwrite(&first_blocksize, sizeof(int), 1, ifs);
	index = 0;
	fwrite(&index, sizeof(int), 1, ifs);
	for(i = 0; i < first_blocksize; i++) {
	  fwrite(&lowestmodes[((intsign+1)%2)*first_blocksize+i], sizeof(double), 1, ifs);
	}
	fclose(ifs);
      }
    }
  }
  else{ 
    /* they are not equal */
    /* determine the sector with not trivial topology */
    if(lowestmodes[0] < lowestmodes[first_blocksize]){
      intsign = 0;
    }
    else{
      intsign = 1;
    }
    
    if(g_proc_id == g_stdio_proc){
      printf("Computing now up to %d modes in the sector with %s chirality\n", 
	     (*nr_of_eigenvalues_ov), intsign ? "negative" : "positive");
      fflush(stdout);
    }

    /* Here we set the (absolute) precision to be  */
    /* such that we can compare to the lowest mode */
    /* in the other sector                         */

    prec = (lowestmodes[first_blocksize*((intsign+1)%2)])*1.e-1;

    eigenvalues_ov = (double*)malloc((*nr_of_eigenvalues_ov)*sizeof(double));

    /* Copy the allready computed eigenvectors_ov */
    for(i = 0; i < first_blocksize; i++) { 
      assign(&eigenvectors_ov[i], &lowvectors[(first_blocksize*intsign+i)*VOLUMEPLUSRAND],N2);
      eigenvalues_ov[i] = lowestmodes[first_blocksize*intsign+i];
    }

#ifdef MPI
    atime = MPI_Wtime();
#endif

    blocksize = 3;
    j_min = 8; j_max = 16;
    converged = first_blocksize;
    for(i = first_blocksize; i < (*nr_of_eigenvalues_ov); i+=3) { 

      if((i + blocksize) > (*nr_of_eigenvalues_ov)) {
	blocksize = (*nr_of_eigenvalues_ov) - i;
      }

      /* Fill up the rest with random spinor fields  */
      /* and project it to the corresponding sector  */
      for(v0dim = i; v0dim < i+blocksize; v0dim++){
	random_spinor_field(&eigenvectors_ov[v0dim*VOLUMEPLUSRAND],N2,0);
	Proj(&eigenvectors_ov[v0dim*VOLUMEPLUSRAND], &eigenvectors_ov[v0dim*VOLUMEPLUSRAND],N2, intsign);
      }
      v0dim = blocksize;
      returncode = 0;

      /* compute minimal eigenvalues */
#ifdef MPI
      /*      pjdher(VOLUME*sizeof(spinor)/sizeof(_Complex double), VOLUMEPLUSRAND*sizeof(spinor)/sizeof(_Complex double),
	     shift, prec, omega, n_omega, ev_tr,
	     i+blocksize, j_max, j_min, 
	     max_iterations, blocksize, blockwise, v0dim, (_Complex double*)(&eigenvectors_ov[i*VOLUMEPLUSRAND]),
	     CG, solver_it_max,
	     threshold_min, decay_min, verbosity,
	     &converged, (_Complex double*) eigenvectors_ov, eigenvalues_ov,
	     &returncode, JD_MINIMAL, 1, use_AV,
	     Operator[intsign]);*/
#else
      jdher(VOLUME*sizeof(spinor)/sizeof(_Complex double),
	    VOLUMEPLUSRAND*sizeof(spinor)/sizeof(_Complex double),
	    shift, prec, blocksize, j_max, j_min,
	    max_iter, blocksize, blockwise, v0dim, (_Complex double*) &eigenvectors_ov[i*VOLUMEPLUSRAND],
	    CG, solver_it_max,
	    threshold_min, decay_min, verbosity,
	    &converged, (_Complex double*) eigenvectors_ov,
	    eigenvalues_ov,
	    &returncode, JD_MINIMAL, 1,
	    Operator[intsign]);
#endif
      /* Save eigenvectors_ov temporary    */
      /* in order to be able to restart */
      for (k=i; k < converged; k++){
	if(intsign == 0){
	  sprintf(filename, "eigenvector_of_Dplus.%.2d.%s.%.4d", k, conf_filename, nstore);
	}
	else{
	  sprintf(filename, "eigenvector_of_Dminus.%.2d.%s.%.4d", k, conf_filename, nstore);
	}
	/*	write_spinorfield(&eigenvectors_ov[k*VOLUMEPLUSRAND], filename);*/
      }

      /* order the eigenvalues_ov and vectors */
      for(k = 0; k < converged; k++) {
	idx[k] = k;
      }
      /*      quicksort(converged, eigenvalues_ov, idx);*/

      /* Check whether the index is detemined */
      index = 0;
      for(k = 0; k < converged; k++) { 
	absdifference = fabs(lowestmodes[first_blocksize*((intsign+1)%2)] - eigenvalues_ov[k]);
	if(absdifference < 0.1*lowestmodes[first_blocksize*((intsign+1)%2)]) {
	  /* We have found the first non zero */
	  if(k < converged-1) {
	    determined = 1;
	    break;
	  }
	  else {
	    blocksize = 1;
	    shift = eigenvalues_ov[converged-1];
	  }
	}
	else {
	  index++;
	}
      }
      /* If we have determined the index or */
      /* hit the maximal number of ev       */
      if(determined == 1 || converged == (*nr_of_eigenvalues_ov)) {
	break;
      }
      else if(g_proc_id == g_stdio_proc) {
	if(blocksize != 1) {
	  printf("Index %s (or equal) than %s%d, continuing!\n\n", 
		 intsign ? "lower" : "bigger", 
		 intsign ? "-" : "+", index);
	  fflush( stdout );
	}
	else {
	  printf("Index is %s%d, one non zero is missing, continuing!\n\n", 
		 intsign ? "-" : "+", index);
	  fflush( stdout );
	}
      }
    }

#ifdef MPI
    etime = MPI_Wtime();
#endif

    /* Save the eigenvectors_ov */
    for(i = 0; i < converged; i++){
      eval[i] = 2.*(1.+ov_s)*eigenvalues_ov[i];
      if(intsign == 0){
	sprintf(filename, "eigenvector_of_Dplus.%.2d.%s.%.4d", i, conf_filename, nstore);
      }
      else{
	sprintf(filename, "eigenvector_of_Dminus.%.2d.%s.%.4d", i, conf_filename, nstore);
      }
      /*      write_spinorfield(&eigenvectors_ov[idx[i]*VOLUMEPLUSRAND], filename);*/
    }

    /* Some Output */
    if(g_proc_id == g_stdio_proc) {
      printf("Index is %s%d!\n", intsign ? "-" : "+", index);
#ifdef MPI
      printf("Zero modes determined in %f sec!\n", etime-atime);
#endif
    }
    if(g_proc_id == 0) {
      sprintf(filename, "eigenvalues_of_overlap_proj.%s.%.4d", conf_filename, nstore);
      ifs = fopen(filename, "w");
      printf("\nThe following lowest modes have been computed:\n");
      fprintf(ifs, "Index is %s%d!\n\n", intsign ? "-" : "+", index);
      for(k = 0; k < 2; k++) {
	if(k == intsign) {
	  for (i=0; i < converged; i++) {
	    fprintf(ifs, "%d %e %s\n", i, creal(eval[i]), intsign ? "negative" : "positive");
	    printf("%d %e %s\n", i, creal(eval[i]), intsign ? "negative" : "positive");
	  }
	}
	else {
	  for(i = 0; i < first_blocksize; i++) {
	    lowestmodes[((intsign+1)%2)*first_blocksize+i] = 2.*(1.+ov_s)*lowestmodes[((intsign+1)%2)*first_blocksize+i];
	    fprintf(ifs, "%d %e %s\n", i, lowestmodes[((intsign+1)%2)*first_blocksize+i], intsign ? "positive" : "negative");
	    printf("%d %e %s\n", i, lowestmodes[((intsign+1)%2)*first_blocksize+i], intsign ? "positive" : "negative");
	  }
	}
      }
      fclose(ifs);
      if(intsign != 0) signed_index = -index;
      else signed_index = index;
      for(k = 0; k < 2; k++) {
	sprintf(filename, "eigenvalues_of_D%s.%s.%.4d", 
		k ? "minus" : "plus", conf_filename, nstore);
	ifs = fopen(filename, "w");
	if(k == intsign) {
	  fwrite(&converged, sizeof(int), 1, ifs);
	  fwrite(&signed_index, sizeof(int), 1, ifs);
	  for (i=index; i < converged; ++i)
	  {
	    double eval_re = creal(eval[i]);
	    fwrite(&eval_re, sizeof(double), 1, ifs);
	  }
	}
	else {
	  fwrite(&first_blocksize, sizeof(int), 1, ifs);
	  fwrite(&signed_index, sizeof(int), 1, ifs);
	  for(i = 0; i < first_blocksize; i++) {
	    fwrite(&lowestmodes[((intsign+1)%2)*first_blocksize+i], sizeof(double), 1, ifs);
	  }
	}
	fclose(ifs);
      }
    }
  }

  switch_on_adaptive_precision = 0; 
  /* Free memory */
  free(eigenvectors_ov_);
  free(lowvectors_);
  free(eval);
  free(eigenvalues_ov);
  free(idx);
}
