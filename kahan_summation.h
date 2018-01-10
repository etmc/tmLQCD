/***********************************************************************
 *
 * Copyright (C) 2018 Bartosz Kostrzewa 
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
 ************************************************************************/

#ifndef KAHAN_SUMMATION_H
#define KAHAN_SUMMATION_H

typedef struct kahan_re_t {
  double kc;
  double ks;
  double ts;
  double tr;
  double tt;
} kahan_re_t;

typedef struct kahan_cplx_t {
  complex double kc;
  complex double ks;
  complex double ts;
  complex double tr;
  complex double tt;
} kahan_cplx_t;

static inline void kahan_sum_re_step(const double in, kahan_re_t * const acc){
  acc->tr = in + acc->kc;
  acc->ts = acc->tr + acc->ks;
  acc->tt = acc->ts - acc->ks;
  acc->ks = acc->ts;
  acc->kc = acc->tr - acc->tt;
}

static inline double kahan_sum_re_final(const kahan_re_t * acc){
  return( acc->kc + acc->ks );
}

static inline void kahan_sum_cplx_step(const complex double in, kahan_cplx_t * acc){
  acc->tr = in + acc->kc;
  acc->ts = acc->tr + acc->ks;
  acc->tt = acc->ts - acc->ks;
  acc->ks = acc->ts;
  acc->kc = acc->tr - acc->tt;
}

static inline complex double kahan_sum_cplx_final(const kahan_cplx_t * const acc){
  return( acc->kc + acc->ks );
}

static inline kahan_re_t new_kahan_re(){
  kahan_re_t ret;
  ret.kc = 0.0;
  ret.ks = 0.0;
  ret.ts = 0.0;
  ret.tr = 0.0;
  ret.tt = 0.0;
  return(ret);
}

static inline kahan_cplx_t new_kahan_cplx(){
  kahan_cplx_t ret;
  ret.kc = 0.0;
  ret.ks = 0.0;
  ret.ts = 0.0;
  ret.tr = 0.0;
  ret.tt = 0.0;
  return(ret);
}

static inline void reset_kahan_re(kahan_re_t * const in){
  in->kc = 0.0;
  in->ks = 0.0;
  in->ts = 0.0;
  in->tr = 0.0;
  in->tt = 0.0;
}

static inline void reset_kahan_cplx(kahan_cplx_t * const in){
  in->kc = 0.0;
  in->ks = 0.0;
  in->ts = 0.0;
  in->tr = 0.0;
  in->tt = 0.0;
}

#endif
