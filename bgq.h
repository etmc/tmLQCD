#ifndef _BGQ_H
#define _BGQ_H

#define regtype vector4double

inline void vec_load2(vector4double * r, su3_vector * const phi) {
  r[0] = vec_ld2(0, (double*) &phi->c0); 
  r[1] = vec_ld2(0, (double*) &phi->c1);
  r[2] = vec_ld2(0, (double*) &phi->c2);
  return;
}

inline void vec_store2(su3_vector * const phi, vector4double * r) {
  vec_st2(r[0], 0, (double*) &phi->c0);
  vec_st2(r[1], 0, (double*) &phi->c1);
  vec_st2(r[2], 0, (double*) &phi->c2);
  return;
}

inline void vec_add2(vector4double * restrict r, vector4double * restrict s) {
#pragma disjoint(*s, *r)
#pragma unroll(6)
  for(int i = 0; i < 6; i++) {
    r[i] = vec_add(r[i], s[i]);
  }
  return;
}

inline void vec_cmplx_mul_double2(complex double * restrict c, vector4double * restrict r, vector4double * restrict tmp) {
  __alignx(32, c);
  tmp[0] = vec_ld2(0, (double*) c);
#pragma unroll(6)
  for(int i = 0; i < 6; i++) {
    r[i] = vec_xmul(r[6+i], tmp[0]);
  }
#pragma unroll(6)
  for(int i = 0; i < 6; i++) {
    r[i] = vec_xxnpmadd(tmp[0], r[6+i], r[i]);
  }
  return;
}

// multiplies one su3 matrix with two su3_vectors
// the first of which stored in r[0-2]
// and the second one in r[3-5]
//
// the resulting two vectors are stored in
// r[0-5]
//
// this routine uses only half of the 4 doubles in vector4double

inline void vec_su3_multiply_double2(su3 * const restrict u, vector4double * restrict U, vector4double * restrict r) {
#pragma disjoint(*U, *r)
  __alignx(32, u);
  __alignx(32, U);
  __alignx(32, r);

  U[0] = vec_ld2(0, (double*) &u->c00);
  U[3] = vec_ld2(0, (double*) &u->c01);
  U[6] = vec_ld2(0, (double*) &u->c02);
  U[1] = vec_ld2(0, (double*) &u->c10);
  U[4] = vec_ld2(0, (double*) &u->c11);
  U[7] = vec_ld2(0, (double*) &u->c12);
  U[2] = vec_ld2(0, (double*) &u->c20);
#pragma unroll(3)
  for(int i = 0; i < 3; i++) {
    r[6+i] = vec_xmul(r[0], U[i]);
    r[9+i] = vec_xmul(r[3], U[i]);
  }
#pragma unroll(3)
  for(int i = 0; i < 3; i++) {
    r[6+i] = vec_xxnpmadd(U[i], r[0], r[6+i]);
    r[9+i] = vec_xxnpmadd(U[i], r[3], r[9+i]);
  }
  U[5] = vec_ld2(0, (double*) &u->c21);
#pragma unroll(3)
  for(int i = 0; i < 3; i++) {
    r[6+i] = vec_xmadd(r[1], U[3+i], r[6+i]);
    r[9+i] = vec_xmadd(r[4], U[3+i], r[9+i]);
  }
#pragma unroll(3)
  for(int i = 0; i < 3; i++) {
    r[6+i] = vec_xxnpmadd(U[3+i], r[1], r[6+i]);
    r[9+i] = vec_xxnpmadd(U[3+i], r[4], r[9+i]);
  }
  U[8] = vec_ld2(0, (double*) &u->c22);
#pragma unroll(3)
  for(int i = 0; i < 3; i++) {
    r[6+i] = vec_xmadd(r[2], U[6+i], r[6+i]);
    r[9+i] = vec_xmadd(r[5], U[6+i], r[9+i]);
  }
#pragma unroll(3)
  for(int i = 0; i < 3; i++) {
    r[6+i] = vec_xxnpmadd(U[6+i], r[2], r[6+i]);
    r[9+i] = vec_xxnpmadd(U[6+i], r[5], r[9+i]);
  }
  return;
}


// might not be optimal for pipeline as result is 
// re-used in the next line.
inline void vec_su3_multiply_double2b(su3 * const u, vector4double * U, vector4double * r) {
  __alignx(32, u);
  __alignx(32, U);
  __alignx(32, r);
  U[0] = vec_ld2(0, (double*) &u->c00);
  U[1] = vec_ld2(0, (double*) &u->c01);
  U[2] = vec_ld2(0, (double*) &u->c02);

  r[6] = vec_xmul(r[0], U[0]);
  r[6] = vec_xxnpmadd(U[0], r[0], r[6]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[6] = vec_xmadd(r[i], U[i], r[6]);
    r[6] = vec_xxnpmadd(U[i], r[i], r[6]);
  }
  r[9] = vec_xmul(r[3], U[0]);
  r[9] = vec_xxnpmadd(U[0], r[3], r[9]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[9] = vec_xmadd(r[3+i], U[i], r[9]);
    r[9] = vec_xxnpmadd(U[i], r[3+i], r[9]);
  }

  U[0] = vec_ld2(0, (double*) &u->c10);
  U[1] = vec_ld2(0, (double*) &u->c11);
  U[2] = vec_ld2(0, (double*) &u->c12);

  r[7] = vec_xmul(r[0], U[0]);
  r[7] = vec_xxnpmadd(U[0], r[0], r[7]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[7] = vec_xmadd(r[i], U[i], r[7]);
    r[7] = vec_xxnpmadd(U[i], r[i], r[7]);
  }
  r[10] = vec_xmul(r[3], U[0]);
  r[10] = vec_xxnpmadd(U[0], r[3], r[10]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[10] = vec_xmadd(r[3+i], U[i], r[10]);
    r[10] = vec_xxnpmadd(U[i], r[3+i], r[10]);
  }

  U[0] = vec_ld2(0, (double*) &u->c20);
  U[1] = vec_ld2(0, (double*) &u->c21);
  U[2] = vec_ld2(0, (double*) &u->c22);

  r[8] = vec_xmul(r[0], U[0]);
  r[8] = vec_xxnpmadd(U[0], r[0], r[8]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[8] = vec_xmadd(r[i], U[i], r[8]);
    r[8] = vec_xxnpmadd(U[i], r[i], r[8]);
  }
  r[11] = vec_xmul(r[3], U[0]);
  r[11] = vec_xxnpmadd(U[0], r[3], r[11]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[11] = vec_xmadd(r[3+i], U[i], r[11]);
    r[11] = vec_xxnpmadd(U[i], r[3+i], r[11]);
  }
  return;
}


#endif