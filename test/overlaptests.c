#include"lime.h"
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "su3spinor.h"
#include "linalg_eo.h"
#include "start.h"
#ifdef MPI
# include "xchange/xchange.h"
#endif
#include "read_input.h"
#include "boundary.h"
#include "linalg/convert_eo_to_lexic.h"
#include "operator/Dov_psi.h"

#include "overlaptests.h"
#include "gamma.h"

void ov_check_alloc(void * pS) {

  if (pS==NULL) {
    fprintf(stderr, "Error: could not allocate memory for spinor field");
    exit(EXIT_FAILURE);
  }

}

spinor** ov_alloc_spinors(int n) {

  spinor *_s, **s;
  int i;

  s  = (spinor**)calloc(n, sizeof(spinor*));
  ov_check_alloc(s);
#if ( defined SSE || defined SSE2 || defined SSE3)
  _s = malloc((n*VOLUMEPLUSRAND)*sizeof(spinor)+ALIGN_BASE+sizeof(spinor*));
  ov_check_alloc(_s);
  s[0] = (spinor *)(((unsigned long int)(_s) + sizeof(spinor*) + ALIGN_BASE)&~ALIGN_BASE);
  *(((spinor**)s[0])-1) = _s;
#else
  s[0] = malloc(n*VOLUMEPLUSRAND*sizeof(spinor));
  ov_check_alloc(s[0]);
#endif

  for(i = 1; i < n; i++)
    s[i] = s[0]+VOLUMEPLUSRAND;

  return s;
}

void ov_free_spinors(spinor** s) {
	
#if (defined SSE3 || defined SSE2 || defined SSE)
  free(*(((spinor**)s[0])-1));
#else
  free(s[0]);
#endif
  free(s);

}

spinor* ov_alloc_spinor(void) {

  spinor *s, *_s;

#if ( defined SSE || defined SSE2 || defined SSE3)
  _s = malloc(sizeof(spinor*)+VOLUMEPLUSRAND*sizeof(spinor)+ALIGN_BASE);
  ov_check_alloc(_s);
  s = (spinor *)(((unsigned long int)(_s) + sizeof(spinor*) + ALIGN_BASE)&~ALIGN_BASE);
  *(((spinor**)s)-1) = _s;
#else
  s = malloc(VOLUMEPLUSRAND*sizeof(spinor));
  ov_check_alloc(s);
#endif

  return s;
}

void ov_free_spinor(spinor *s) {

#if ( defined SSE || defined SSE2 || defined SSE3)
  free(*(((spinor**)s)-1));
#else
  free(s);
#endif

}


/* col sum norm of operator in colour and spinor space */
double ov_operator_colsumnorm(spinor *s[4][3], int k)
{
  double norm = 0.0, nrm;

  for (int i=0; i<4; ++i)
    for (int j=0; j<3; ++j)
    {
      _spinor_norm_l1(nrm, s[i][j][k]);
      if (nrm > norm)
	norm = nrm;
    }
  return norm;
}

void ov_check_locality() {

  double norm, *maxnorm, *minnorm, *avgnorm;
  int i, j, k, x, x_taxi, y, y_taxi, z, z_taxi, t, t_taxi, maxtaxi, *samples, taxi;
  spinor *s[4][3];

  /* evaluate Dov(psi) */
  for (i=0; i<4; i++)
    for (j=0; j<3; j++) {

      /* get memory for the spinor */
      s[i][j] = ov_alloc_spinor();

      /* create delta source at origin */
      source_spinor_field(g_spinor_field[1], g_spinor_field[0], i, j);
      convert_eo_to_lexic(g_spinor_field[2], g_spinor_field[1], g_spinor_field[0]);

      /* apply Dov */
      Dov_psi(s[i][j], g_spinor_field[2]);
    }

  /* init locality table */
  maxtaxi = (LX/2)+(LY/2)+(LZ/2)+T/2;
  maxnorm = (double*)calloc(maxtaxi+1, sizeof(double));
  minnorm = (double*)calloc(maxtaxi+1, sizeof(double));
  avgnorm = (double*)calloc(maxtaxi+1, sizeof(double));
  samples = (int*)calloc(maxtaxi+1, sizeof(int));
  for(i = 0; i <= maxtaxi; i++) {
    maxnorm[i] = 0.;
    minnorm[i] = 1.0e100;
    avgnorm[i] = 0.;
    samples[i] = 0;
  }

  /* fill locality table */
  printf("// beginning locality test\n");
  for(x=0; x<LX; x++) {
    x_taxi = (x > LX/2) ? LX-x : x;
    for(y = 0; y < LY; y++){
      y_taxi =  (y > LY/2) ? LY-y : y;
      for(z = 0; z < LZ; z++){
	z_taxi = (z > LZ/2) ? LZ-z : z;
	for(t = 0; t < T; t++){
	  t_taxi = (t > T/2) ? T - t : t;
	  taxi = x_taxi + y_taxi + z_taxi + t_taxi;
	  k = g_ipt[t][x][y][z];

	  norm = ov_operator_colsumnorm(s, k);

	  // statistics
	  if (norm > maxnorm[taxi])
	    maxnorm[taxi] = norm;
	  if (norm < minnorm[taxi])
	    minnorm[taxi] = norm;
	  avgnorm[taxi] += norm;
	  samples[taxi]++;
	}
      }
    }
  }

  /* print locality table */
  printf("// locality check of overlap operator\n");
  printf("// taxi | max norm     | avg norm     | min norm     | #samples\n");
  for(i = 0; i <= maxtaxi; i++)
    printf("%7d   %10.6le   %10.6le   %10.6le   %8d\n", i, maxnorm[i], (double)(avgnorm[i]/samples[i]), minnorm[i], samples[i]);
  printf("\n");

  /* free memory */
  free(maxnorm);
  free(minnorm);
  free(avgnorm);
  free(samples);
  for (i=0; i<4; i++)
    for (j=0; j<3; j++)
      ov_free_spinor(s[i][j]);

}

void ov_matrix4x4_diff(matrix4x4 result, matrix4x4 left, matrix4x4 right)
{
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      result[i][j] = left[i][j] - right[i][j];
}


double ov_matrix4x4_rowsumnorm(matrix4x4 A) {
	
  double norm, nrm;
  int i, j;

  norm = 0.0;
  for (i=0; i<4; i++) {

    nrm = 0.0;
    for (j=0; j<4; j++)
      nrm += cabs(A[i][j]);

    if (nrm > norm)
      norm = nrm;
  }

  return norm;
}

void ov_matrix12x12_diff(matrix12x12 result, matrix12x12 left, matrix12x12 right)
{
  for (int i = 0; i < 12; ++i)
    for (int j = 0; j < 12; ++j)
      result[i][j] = left[i][j] - right[i][j];
}

double ov_matrix12x12_rowsumnorm(matrix12x12 A) {
	
  double norm, nrm;
  int i, j;

  norm = 0.0;
  for (i=0; i<12; i++) {

    nrm = 0.0;
    for (j=0; j<12; j++)
      nrm += cabs(A[i][j]);

    if (nrm > norm)
      norm = nrm;
  }

  return norm;
}
/* compares the operator with the one given in pFileName */
void ov_compare_4x4(const char * pFileName) {

  double norm, rel, *max_rel, *max_abs, Max_rel = 0.0, Max_abs = 0.0;
  int i, j, k, x, x_taxi, y, y_taxi, z, z_taxi, t, t_taxi, maxtaxi, taxi;
  spinor *s[4];
  matrix4x4 mat, mat2, diff;
  FILE * pCompare;

  /* evaluate Dov(psi) */
  for (i=0; i<4; i++) {

    /* get memory for the spinor */
    s[i] = ov_alloc_spinor();

    /* create delta source at origin */
    source_spinor_field(g_spinor_field[1], g_spinor_field[0], i, 0);
    convert_eo_to_lexic(g_spinor_field[2], g_spinor_field[1], g_spinor_field[0]);

    /* apply Dov */
    Dov_psi(s[i], g_spinor_field[2]);
  }

  /* init locality table */
  maxtaxi = (LX/2)+(LY/2)+(LZ/2)+T/2;
  max_abs = (double*)calloc(maxtaxi+1, sizeof(double));
  max_rel = (double*)calloc(maxtaxi+1, sizeof(double));
  for(i = 0; i <= maxtaxi; i++) {
    max_abs[i] = 0.0;
    max_rel[i] = 0.0;
  }

  /* open file containing operator for comparison */
  pCompare = fopen(pFileName, "r");
  if (pCompare == NULL) {
    fprintf(stderr, "Error: could not open '%s' for comparison of operator\n", pFileName);
    exit(1);
  }

  /* fill locality table */
  if (g_debug_level > 0) {
    printf("// beginning comparison\n");
    fflush(stdout);
  }
  for(t = 0; t < T; t++){
    t_taxi = (t > T/2) ? T - t : t;
    for(x = 0; x < LX; x++){
      x_taxi =  (x > LX/2) ? LX-x : x;
      for(y = 0; y < LY; y++){
	y_taxi = (y > LY/2) ? LY-y : y;
	for(z=0; z<LZ; z++) {
	  z_taxi = (z > LZ/2) ? LZ-z : z;
	  taxi = x_taxi + y_taxi + z_taxi + t_taxi;
	  k = g_ipt[t][x][y][z];

	  for (i=0; i<4; i++) {
	    mat[0][i] = s[i][k].s0.c0;
	    mat[1][i] = s[i][k].s1.c0;
	    mat[2][i] = s[i][k].s2.c0;
	    mat[3][i] = s[i][k].s3.c0;
	  }

	  for (i=0;i<4; i++)
	    for (j=0; j<4; j++)
	      fscanf(pCompare, "%le %le", (double*)&mat2[i][j], (double*)&mat2[i][j] + 1);

	  ov_matrix4x4_diff(diff, mat, mat2);

	  /* statistics */
	  norm = ov_matrix4x4_rowsumnorm(diff);
	  if (norm > max_abs[taxi]) {
	    max_abs[taxi] = norm;
	    if (norm > Max_abs)
	      Max_abs = norm;
	  }
	  rel = (ov_matrix4x4_rowsumnorm(mat) + ov_matrix4x4_rowsumnorm(mat2))/2;
	  if (rel>0.0) {
	    rel = norm/rel;
	    if (rel > max_rel[taxi]) {
	      max_rel[taxi] = rel;
	      if (rel > Max_rel)
		Max_rel = rel;
	    }
	  }
	}
      }
    }
  }

  /* print locality table */
  printf("// comparison of overlap operator to %s\n", pFileName);
  printf(" - maximum absolute deviation: %.4le\n", Max_abs);
  printf(" - maximum relative deviation: %.4le\n", Max_rel);
  printf("// taxi | max abs     | max rel\n");
  for(i = 0; i <= maxtaxi; i++)
    printf("%7d   %10.6le   %10.6le\n", i, max_abs[i], max_rel[i]);
  printf("\n");

  /* close file */
  fclose(pCompare);

  /* free memory */
  free(max_abs);
  free(max_rel);
  for (i=0; i<4; i++)
    ov_free_spinor(s[i]);

}

/* compares the operator with the one given in pFileName */
void ov_compare_12x12(const char * pFileName) {

  double norm, rel, *max_rel, *max_abs, Max_rel = 0.0, Max_abs = 0.0;
  int i, j, k, x, x_taxi, y, y_taxi, z, z_taxi, t, t_taxi, maxtaxi, taxi;
  spinor *s[4][3];
  matrix12x12 mat, mat2, diff;
  FILE * pCompare;

  /* evaluate Dov(psi) */
  for (i=0; i<4; i++)
    for (j=0; j<3; j++) {

      /* get memory for the spinor */
      s[i][j] = ov_alloc_spinor();

      /* create delta source at origin */
      source_spinor_field(g_spinor_field[1], g_spinor_field[0], i, j);
      convert_eo_to_lexic(g_spinor_field[2], g_spinor_field[1], g_spinor_field[0]);

      /* apply Dov */
      Dov_psi(s[i][j], g_spinor_field[2]);
    }

  /* init locality table */
  maxtaxi = (LX/2)+(LY/2)+(LZ/2)+T/2;
  max_abs = (double*)calloc(maxtaxi+1, sizeof(double));
  max_rel = (double*)calloc(maxtaxi+1, sizeof(double));
  for(i = 0; i <= maxtaxi; i++) {
    max_abs[i] = 0.0;
    max_rel[i] = 0.0;
  }

  /* open file containing operator for comparison */
  pCompare = fopen(pFileName, "r");
  if (pCompare == NULL) {
    fprintf(stderr, "Error: could not open '%s' for comparison of operator\n", pFileName);
    exit(1);
  }

  /* fill locality table */
  if (g_debug_level > 0) {
    printf("// beginning comparison\n");
    fflush(stdout);
  }
  for(t = 0; t < T; t++){
    t_taxi = (t > T/2) ? T - t : t;
    for(x = 0; x < LX; x++){
      x_taxi =  (x > LX/2) ? LX-x : x;
      for(y = 0; y < LY; y++){
	y_taxi = (y > LY/2) ? LY-y : y;
	for(z=0; z<LZ; z++) {
	  z_taxi = (z > LZ/2) ? LZ-z : z;
	  taxi = x_taxi + y_taxi + z_taxi + t_taxi;
	  k = g_ipt[t][x][y][z];

	  for (j=0; j<3; j++)
	    for (i=0; i<4; i++) {
	      mat[0][i+4*j] = s[i][j][k].s0.c0;
	      mat[1][i+4*j] = s[i][j][k].s1.c0;
	      mat[2][i+4*j] = s[i][j][k].s2.c0;
	      mat[3][i+4*j] = s[i][j][k].s3.c0;
	      mat[4][i+4*j] = s[i][j][k].s0.c1;
	      mat[5][i+4*j] = s[i][j][k].s0.c1;
	      mat[6][i+4*j] = s[i][j][k].s1.c1;
	      mat[7][i+4*j] = s[i][j][k].s2.c1;
	      mat[8][i+4*j] = s[i][j][k].s3.c2;
	      mat[9][i+4*j] = s[i][j][k].s1.c2;
	      mat[10][i+4*j] = s[i][j][k].s2.c2;
	      mat[11][i+4*j] = s[i][j][k].s3.c2;
	    }

	  for (i=0;i<12; i++)
	    for (j=0; j<12; j++)
	      fscanf(pCompare, "%le %le", (double*)&mat2[i][j], (double*)&mat2[i][j] + 1);

	  ov_matrix12x12_diff(diff, mat, mat2);

	  /* statistics */
	  norm = ov_matrix12x12_rowsumnorm(diff);
	  if (norm > max_abs[taxi]) {
	    max_abs[taxi] = norm;
	    if (norm > Max_abs)
	      Max_abs = norm;
	  }
	  rel = (ov_matrix12x12_rowsumnorm(mat) + ov_matrix12x12_rowsumnorm(mat2))/2;
	  if (rel>0.0) {
	    rel = norm/rel;
	    if (rel > max_rel[taxi]) {
	      max_rel[taxi] = rel;
	      if (rel > Max_rel)
		Max_rel = rel;
	    }
	  }
	}
      }
    }
  }

  /* print locality table */
  printf("// comparison of overlap operator to %s\n", pFileName);
  printf(" - maximum absolute deviation: %.4le\n", Max_abs);
  printf(" - maximum relative deviation: %.4le\n", Max_rel);
  printf("// taxi | max abs     | max rel\n");
  for(i = 0; i <= maxtaxi; i++)
    printf("%7d   %10.6le   %10.6le\n", i, max_abs[i], max_rel[i]);
  printf("\n");

  /* close file */
  fclose(pCompare);

  /* free memory */
  free(max_abs);
  free(max_rel);
  for (i=0; i<4; i++)
    for (j=0; j<3; j++)
      ov_free_spinor(s[i][j]);

}

/* saves the operator to the given filename */
void ov_save_12x12(const char * pFileName) {

  int i, j, k, x, y, z, t;
  spinor *s[4][3];
  matrix12x12 mat;
  FILE * pCompare;

  /* evaluate Dov(psi) */
  for (i=0; i<4; i++)
    for (j=0; j<3; j++) {

      /* get memory for the spinor */
      s[i][j] = ov_alloc_spinor();

      /* create delta source at origin */
      source_spinor_field(g_spinor_field[1], g_spinor_field[0], i, j);
      convert_eo_to_lexic(g_spinor_field[2], g_spinor_field[1], g_spinor_field[0]);

      /* apply Dov */
      Dov_psi(s[i][j], g_spinor_field[2]);
    }

  /* open file for storing the operator */
  pCompare = fopen(pFileName, "w");
  if (pCompare == NULL) {
    fprintf(stderr, "Error: could not open '%s' for writing the operator\n", pFileName);
    exit(1);
  }

  for(t = 0; t < T; t++){
    for(x = 0; x < LX; x++){
      for(y = 0; y < LY; y++){
	for(z=0; z<LZ; z++) {
	  k = g_ipt[t][x][y][z];

	  for (j=0; j<3; j++)
	    for (i=0; i<4; i++) {
	      mat[0][i+4*j] = s[i][j][k].s0.c0;
	      mat[1][i+4*j] = s[i][j][k].s1.c0;
	      mat[2][i+4*j] = s[i][j][k].s2.c0;
	      mat[3][i+4*j] = s[i][j][k].s3.c0;
	      mat[4][i+4*j] = s[i][j][k].s0.c1;
	      mat[5][i+4*j] = s[i][j][k].s0.c1;
	      mat[6][i+4*j] = s[i][j][k].s1.c1;
	      mat[7][i+4*j] = s[i][j][k].s2.c1;
	      mat[8][i+4*j] = s[i][j][k].s3.c2;
	      mat[9][i+4*j] = s[i][j][k].s1.c2;
	      mat[10][i+4*j] = s[i][j][k].s2.c2;
	      mat[11][i+4*j] = s[i][j][k].s3.c2;
	    }

	  for (i=0;i<12; i++)
	    for (j=0; j<12; j++)
	      fprintf(pCompare, "%.20le %.20le ", creal(mat[i][j]), cimag(mat[i][j]));
	}
      }
    }
  }

  /* close file */
  fclose(pCompare);

  /* free memory */
  for (i=0; i<4; i++)
    for (j=0; j<3; j++)
      ov_free_spinor(s[i][j]);

}
void ov_check_q_over_sqrt_q_sqr(spinor * const P, spinor * const S) {
  /*
    spinor *s;
    static int n_cheby = 0;
    static int rec_coefs = 1;

    if(n_cheby != ov_n_cheby || rec_coefs) {
    if(ov_cheby_coef != NULL)
    free(ov_cheby_coef);
    ov_cheby_coef = (double*)malloc(ov_n_cheby*sizeof(double));
    chebyshev_coefs(ev_minev, 1., ov_cheby_coef, ov_n_cheby, -0.5);
    rec_coefs = 0;
    n_cheby = ov_n_cheby;
    }

    s = ov_alloc_spinor();

    Q_over_sqrt_Q_sqr(s, ov_cheby_coef, ov_n_cheby, S, ev_qnorm, ev_minev);
    Q_over_sqrt_Q_sqr(P, ov_cheby_coef, ov_n_cheby, s, ev_qnorm, ev_minev);

    diff(P, P, S, VOLUME);

    free(s);
  */
}

void ov_print_spinor(spinor * pS) {

  printf("Color 0:                            | Color 1:                            | Color 2:                            \n");
  printf("------------------------------------+-------------------------------------+-------------------------------------\n");
  printf("%16.19le %+16.19le I | %16.9le %+16.9le I | %16.9le %+16.9le I\n", (double)creal(pS->s0.c0), (double)cimag(pS->s0.c0), (double)creal(pS->s0.c1), (double)cimag(pS->s0.c1), (double)creal(pS->s0.c2), (double)cimag(pS->s0.c2));
  printf("%16.19le %+16.19le I | %16.9le %+16.9le I | %16.9le %+16.9le I\n", (double)creal(pS->s1.c0), (double)cimag(pS->s1.c0), (double)creal(pS->s1.c1), (double)cimag(pS->s1.c1), (double)creal(pS->s1.c2), (double)cimag(pS->s1.c2));
  printf("%16.19le %+16.19le I | %16.9le %+16.9le I | %16.9le %+16.9le I\n", (double)creal(pS->s2.c0), (double)cimag(pS->s2.c0), (double)creal(pS->s2.c1), (double)cimag(pS->s2.c1), (double)creal(pS->s2.c2), (double)cimag(pS->s2.c2));
  printf("%16.19le %+16.19le I | %16.9le %+16.9le I | %16.9le %+16.9le I\n", (double)creal(pS->s3.c0), (double)cimag(pS->s3.c0), (double)creal(pS->s3.c1), (double)cimag(pS->s3.c1), (double)creal(pS->s3.c2), (double)cimag(pS->s3.c2));

}


void ov_check_operator(int t, int x, int y, int z) {

  /* Create delta source at origin */
  source_spinor_field(g_spinor_field[1], g_spinor_field[0], 0, 0);
  convert_eo_to_lexic(g_spinor_field[2], g_spinor_field[1], g_spinor_field[0]);

  /* Evaluate Dov(psi) */
  Dov_psi(g_spinor_field[3], g_spinor_field[2]);
  ov_print_spinor(&g_spinor_field[3][g_ipt[t][x][y][z]]);

}

/* Check GW relation with operator norm over the full lattice */
void ov_check_ginsparg_wilson_relation_strong(void) {

  double norm_diff, norm_left, norm_right, norm, max_rel = 0.0, min_left = 1.0e100, min_right = 1.0e100,  max_diff = 0.0, min_norm = 1.0e100;
  int x, y, z, t, i, j, k;
  spinor *S_left[4][3], *S_right[4][3], *S_diff[4][3];

  if (g_debug_level>0) {
    printf("// creating spinor fields and calculating {gamma_5,D} psi and a D gamma_5 D psi\n");
    fflush(stdout);
  }
  for (i=0; i<4; i++)
    for (j=0; j<3; j++) {
			
      if (g_debug_level>1) {
	printf("// spinor field: delta_dirac at %d, delta_color at %d\n", i, j);
	fflush(stdout);
      }
			
      /* get memory for the spinor */
      S_left[i][j]  = ov_alloc_spinor();
      S_right[i][j] = ov_alloc_spinor();
      S_diff[i][j]  = ov_alloc_spinor();

      /* Create delta source at origin */
      source_spinor_field(g_spinor_field[1], g_spinor_field[0], i, j);
      convert_eo_to_lexic(g_spinor_field[2], g_spinor_field[1], g_spinor_field[0]);

      /* S_right = D gamma_5 D psi */
      Dov_psi(g_spinor_field[3], g_spinor_field[2]);
      gamma5(S_left[i][j], g_spinor_field[3], VOLUME);
      Dov_psi(S_right[i][j], S_left[i][j]);

      /* S_left = {gamma_5, D} psi */
      gamma5(g_spinor_field[3], g_spinor_field[2], VOLUME);
      Dov_psi(g_spinor_field[4], g_spinor_field[3]);
      add(S_left[i][j], S_left[i][j], g_spinor_field[4], VOLUME);

      /* S_diff = (S_left-S_right) psi, should be zero (GW relation) */
      diff(S_diff[i][j], S_left[i][j], S_right[i][j], VOLUME);
    }

  /* scan the whole lattice and check GW relation */
  printf("// test of the Ginsparg-Wilson relation:\n");
  if (g_debug_level>0)
    fflush(stdout);
  for(x=0; x<LX; x++)
    for(y = 0; y < LY; y++)
      for(z = 0; z < LZ; z++)
	for(t = 0; t < T; t++) {
	  k = g_ipt[t][x][y][z];
	  norm_diff  = ov_operator_colsumnorm(S_diff, k);
	  norm_left  = ov_operator_colsumnorm(S_left, k);
	  norm_right = ov_operator_colsumnorm(S_right, k);
	  norm = (norm_left+norm_right)/2.;
	  if (norm < min_norm)
	    min_norm = norm;
	  if (norm > 0.0) {
	    norm = norm_diff/norm;
	    if (norm > max_rel)
	      max_rel = norm;
	    if ((norm > 1.8) && (g_debug_level)>=5) {
	      printf("(%d,%d,%d,%d): taxi = %d, rel = %.20le, lr = [%.4le, %.4le], diff = %.4le\n", t, x, y, z, ((x>LX/2) ? LX-x : x)+((y>LY/2) ? LY-y : y)+((z>LZ/2) ? LZ-z : z)+((t>T/2) ? T-t : t), norm, norm_left, norm_right, norm_diff);
	      printf("// left[0][0]:\n");
	      ov_print_spinor(&S_left[0][0][k]);
	      printf("// right[0][0]:\n");
	      ov_print_spinor(&S_right[0][0][k]);
	      printf("// diff[0][0]:\n");
	      ov_print_spinor(&S_diff[0][0][k]);
	    }
	  }
	  if (norm_left < min_left)
	    min_left = norm_left;
	  if (norm_right < min_right)
	    min_right = norm_right;
	  if (norm_diff > max_diff)
	    max_diff = norm_diff;
	}

  /* print results */
  printf(" - maximum absolute deviation: %.4le\n", max_diff);
  printf(" - maximum relative deviation: %.4le\n", max_rel);
  printf(" - minimum mean norm: %.4le\n", min_norm);
  printf(" - minimum norm {gamma_5, D}: %.4le\n", min_left);
  printf(" - minimum norm D gamma_5 D: %.4le\n", min_right);

  /* free memory */
  for (i=0; i<4; i++)
    for (j=0; j<3; j++) {
      ov_free_spinor(S_left[i][j]);
      ov_free_spinor(S_right[i][j]);
      ov_free_spinor(S_diff[i][j]);
    }
}

/* Checks GW relation only by applying Dov to delta(0,0) */
void ov_check_ginsparg_wilson_relation(void) {

  double norm_diff, norm_left, norm_right, norm, max_rel = 0.0, min_left = 1.0e100, min_right = 1.0e100, max_diff = 0.0, min_norm = 1.0e100;
  int x, y, z, t, i;
  spinor *S_left, *S_right, *S_diff;

  /* get memory for the spinor fields */
  S_left  = ov_alloc_spinor();
  S_right = ov_alloc_spinor();
  S_diff  = ov_alloc_spinor();

  /* Create delta source at origin */
  source_spinor_field(g_spinor_field[1], g_spinor_field[0], 0, 0);
  convert_eo_to_lexic(g_spinor_field[2], g_spinor_field[1], g_spinor_field[0]);

  /* S_right = D gamma_5 D */
  Dov_psi(g_spinor_field[3], g_spinor_field[2]);
  gamma5(S_left, g_spinor_field[3], VOLUME);
  Dov_psi(S_right, S_left);

  /* S_left = {gamma_5, D} */
  gamma5(g_spinor_field[3], g_spinor_field[2], VOLUME);
  Dov_psi(g_spinor_field[4], g_spinor_field[3]);
  add(S_left, S_left, g_spinor_field[4], VOLUME);

  /* S_diff = S_left-S_right */
  diff(S_diff, S_left, S_right, VOLUME);

  /* scan the whole lattice */
  printf("// test of the Ginsparg-Wilson relation\n");
  for(x=0; x<LX; x++)
    for(y = 0; y < LY; y++)
      for(z = 0; z < LZ; z++)
	for(t = 0; t < T; t++) {
	  i = g_ipt[t][x][y][z];
 	  _spinor_norm_sq(norm_diff, S_diff[i]); 
 	  _spinor_norm_sq(norm_left, S_left[i]); 
 	  _spinor_norm_sq(norm_right, S_right[i]); 
	  norm_diff = sqrt(norm_diff);
	  norm_left = sqrt(norm_left);
	  norm_right = sqrt(norm_right);
	  norm = norm_left+norm_right;
	  if (norm < min_norm)
	    min_norm = norm;
	  if (norm > 0.0) {
	    norm = 2.*norm_diff/norm;
	    if (norm > max_rel)
	      max_rel = norm;
	  }
	  if (norm_left < min_left)
	    min_left = norm_left;
	  if (norm_right < min_right)
	    min_right = norm_right;
	  if (norm_diff > max_diff)
	    max_diff = norm_diff;
	}

  /* print results */
  printf(" - maximum absoulte deviation: %.4le\n", max_diff);
  printf(" - maximum relative deviation: %.4le\n", max_rel);
  printf(" - minimum mean norm: %4.le\n", min_norm);
  printf(" - minimum norm {gamma_5, D}: %.4le\n", min_left);
  printf(" - minimum norm D gamma_5 D: %.4le\n", min_right);

  /* free memory */
  ov_free_spinor(S_left);
  ov_free_spinor(S_right);
  ov_free_spinor(S_diff);
}

