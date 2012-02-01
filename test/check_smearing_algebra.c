#include <su3.h>
#include <su3adj.h>
#include <expo.h>
#include <smearing/utils.h>

#include <stdio.h>

int main()
{
  su3 Q, U;
  su3adj T;
  _Complex double f0, f1, f2;
  
  /* Positive determinant */
  Q.c00 = -0.2994;
  Q.c01 =  0.5952 + 1.3123 * I;
  Q.c02 = -0.7943 + 0.0913 * I;
  Q.c11 = -1.1430;
  Q.c12 = -2.0025 + 0.2978 * I;
  Q.c22 = +1.4424;
  Q.c10 = conj(Q.c01);
  Q.c20 = conj(Q.c02);
  Q.c21 = conj(Q.c12);
  
  /* Matlab's solution for U = exp(i * Q) */
  U.c00 = +0.3391 - 0.1635 * I;
  U.c01 = -0.2357 + 0.5203 * I;
  U.c02 = +0.5609 + 0.4663 * I;
  U.c10 = -0.0740 - 0.4204 * I;
  U.c11 = -0.7706 - 0.1863 * I;
  U.c12 = +0.1191 - 0.4185 * I;
  U.c20 = +0.5351 - 0.6243 * I;
  U.c21 = +0.1825 + 0.1089 * I;
  U.c22 = -0.5279 - 0.0022 * I;

  fprintf(stdout, "Input matrix -- positive determinant:\n");
  print_su3(&Q);
  fprintf(stdout, "Known result:\n");
  print_su3(&U);
  cayley_hamilton_exponent(&U, &f0, &f1, &f2, &Q);
  fprintf(stdout, "Calculated result:\n");
  print_su3(&U);
  fprintf(stdout, "Produced coefficients:\n");
  fprintf(stdout, "f0: %6.4f + %6.4f * I,   f1: %6.4f + %6.4f * I,   f2: %6.4f + %6.4f * I\n\n", creal(f0), cimag(f0), creal(f1), cimag(f1), creal(f2), cimag(f2));
  
  /* Negative determinant */
  Q.c00 = -0.0186;
  Q.c01 =  0.4756 + 0.4783 * I;
  Q.c02 =  0.3640 + 0.5255 * I;
  Q.c11 =  1.1207;
  Q.c12 = -1.0275 - 0.2772 * I;
  Q.c22 = -1.1021;
  Q.c10 = conj(Q.c01);
  Q.c20 = conj(Q.c02);
  Q.c21 = conj(Q.c12);
  
  /* Matlab's solution for U = exp(i * Q) */
  U.c00 = +0.6722 + 0.1020 * I;
  U.c01 = -0.2735 + 0.2111 * I;
  U.c02 = -0.0124 + 0.6467 * I;
  U.c10 = +0.2438 + 0.3043 * I;
  U.c11 = -0.1096 - 0.7077 * I;
  U.c12 = +0.0254 - 0.5783 * I;
  U.c20 = +0.5715 - 0.2432 * I;
  U.c21 = -0.2758 - 0.5400 * I;
  U.c22 = -0.0149 - 0.4963 * I;

  fprintf(stdout, "Input matrix -- negative determinant:\n");
  print_su3(&Q);
  fprintf(stdout, "Known result:\n");
  print_su3(&U);
  cayley_hamilton_exponent(&U, &f0, &f1, &f2, &Q);
  fprintf(stdout, "Calculated result:\n");
  print_su3(&U);
  fprintf(stdout, "Produced coefficients:\n");
  fprintf(stdout, "f0: %6.4f + %6.4f * I,   f1: %6.4f + %6.4f * I,   f2: %6.4f + %6.4f * I\n\n", creal(f0), cimag(f0), creal(f1), cimag(f1), creal(f2), cimag(f2));
  _trace_lambda(T, Q);
  Q = exposu3(T);
  fprintf(stdout, "With old routine:\n");
  print_su3(&Q);

  /* Almost singular */
  Q.c00 =  0.0000;
  Q.c01 = -0.0009 - 0.0006 * I;
  Q.c02 =  0.0014 - 0.0003 * I;
  Q.c11 =  0.0005;
  Q.c12 =  0.0004 - 0.0010 * I;
  Q.c22 = -0.0005;
  Q.c10 = conj(Q.c01);
  Q.c20 = conj(Q.c02);
  Q.c21 = conj(Q.c12);
  
  /* Matlab's solution for U = exp(i * Q) */
  U.c00 =  1.0000;
  U.c01 =  0.0006 - 0.0009 * I;
  U.c02 =  0.0003 + 0.0014 * I;
  U.c10 = -0.0006 - 0.0009 * I;
  U.c11 =  1.0000;
  U.c12 =  0.0010 + 0.0004 * I;
  U.c20 = -0.0003 + 0.0014 * I;
  U.c21 = -0.0010 + 0.0004 * I;
  U.c22 =  1.0000;

  fprintf(stdout, "Input matrix -- nearly singular:\n");
  print_su3(&Q);
  fprintf(stdout, "Known result:\n");
  print_su3(&U);
  cayley_hamilton_exponent(&U, &f0, &f1, &f2, &Q);
  fprintf(stdout, "Calculated result:\n");
  print_su3(&U);
  fprintf(stdout, "Produced coefficients:\n");
  fprintf(stdout, "f0: %6.4f + %6.4f * I,   f1: %6.4f + %6.4f * I,   f2: %6.4f + %6.4f * I\n\n", creal(f0), cimag(f0), creal(f1), cimag(f1), creal(f2), cimag(f2));

  return 0;
}
