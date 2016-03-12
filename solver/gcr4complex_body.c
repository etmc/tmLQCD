/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *               2010 claude Tadonki
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
 ***********************************************************************/



static void _PSWITCH(init_lgcr)(const int _M, const int _V);
static void _PSWITCH(free_lgcr)();
static _C_TYPE ** _PSWITCH(a)     = NULL; 
static _C_TYPE *  _PSWITCH(_a)    = NULL;
static _F_TYPE *  _PSWITCH(b)     = NULL;
static _C_TYPE *  _PSWITCH(c)     = NULL;
static _C_TYPE ** _PSWITCH(chi)   = NULL;
static _C_TYPE *  _PSWITCH(_chi)  = NULL;
static _C_TYPE ** _PSWITCH(xi)    = NULL;
static _C_TYPE *  _PSWITCH(_xi)   = NULL;
static _C_TYPE *  _PSWITCH(tmp)   = NULL;
static _C_TYPE *  _PSWITCH(rho)   = NULL;
static int _PSWITCH(lgcr_init)   = 0;


int _PSWITCH(gcr4complex)(_C_TYPE * const P, _C_TYPE * const Q, 
                          const int m, const int max_restarts,
                          const double eps_sq, const int rel_prec,
                          const int N, const int parallel, 
                          const int lda, const int precon, _PSWITCH(c_matrix_mult) f) {

  int k, l, restart, i, p=0;
  double norm_sq, err;
  _C_TYPE ctmp;

  _PSWITCH(init_lgcr)(m, lda);

  norm_sq = _PSWITCH(lsquare_norm)(Q, N, parallel);
  if(norm_sq < 1.e-20) {
    norm_sq = 1.;
  }

  for(restart = 0; restart < max_restarts; restart++) {
    f(_PSWITCH(tmp), P);
    _PSWITCH(ldiff)(_PSWITCH(rho), Q, _PSWITCH(tmp), N);
    err = _PSWITCH(lsquare_norm)(_PSWITCH(rho), N, parallel);
    if(g_proc_id == g_stdio_proc && g_debug_level > 2){
      printf("lGCR: %d\t%g true residue\n", p, err); 
      fflush(stdout);
    }
    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq * norm_sq) && (rel_prec == 1))) {
      if(g_proc_id == 0 && g_debug_level > 2) printf("lGCR: %d %e %e %e %e\n", p, err, norm_sq, err/norm_sq, eps_sq);
      return (p);
    }
    for(k = 0; k < m ; k++) {
      if(precon == 0) {
        memcpy(_PSWITCH(xi)[k], _PSWITCH(rho), N*sizeof(_C_TYPE));
      }
      else {
	_PSWITCH(little_mg_precon)(_PSWITCH(xi)[k], _PSWITCH(rho));
      }
      f(_PSWITCH(tmp), _PSWITCH(xi)[k]); 
      /* tmp will become chi[k] */
      for(l = 0; l < k; l++) {
        _PSWITCH(a)[l][k] = _PSWITCH(lscalar_prod)(_PSWITCH(chi)[l], _PSWITCH(tmp), N, parallel);
        _PSWITCH(lassign_diff_mul)(_PSWITCH(tmp), _PSWITCH(chi)[l], _PSWITCH(a)[l][k], N);
      }
      _PSWITCH(b)[k] = sqrt(_PSWITCH(lsquare_norm)(_PSWITCH(tmp), N, parallel));
      _PSWITCH(lmul_r)(_PSWITCH(chi)[k], 1./_PSWITCH(b)[k], _PSWITCH(tmp), N);
      _PSWITCH(c)[k] = _PSWITCH(lscalar_prod)(_PSWITCH(chi)[k], _PSWITCH(rho), N, parallel);
      _PSWITCH(lassign_diff_mul)(_PSWITCH(rho), _PSWITCH(chi)[k], _PSWITCH(c)[k], N);
      err = _PSWITCH(lsquare_norm)(_PSWITCH(rho), N, parallel);
      if(g_proc_id == g_stdio_proc && g_debug_level > 2){
        printf("lGCR: %d\t%g iterated residue\n", p, err); 
        fflush(stdout);
      }
      p++;
      /* Precision reached? */
      if((k == m-1) || ((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*norm_sq) && (rel_prec == 1))) {
        break;
      }
    }
    /* prepare for restart */
    _PSWITCH(c)[k] /= _PSWITCH(b)[k];
    _PSWITCH(lassign_add_mul)(P, _PSWITCH(xi)[k], _PSWITCH(c)[k], N);
    for(l = k-1; l >= 0; --l) {
      for(i = l+1; i <= k; ++i) {
        ctmp  = _PSWITCH(a)[l][i] * _PSWITCH(c)[i];
        _PSWITCH(c)[l] -= ctmp;
      }
      _PSWITCH(c)[l] /= _PSWITCH(b)[l];
      _PSWITCH(lassign_add_mul)(P, _PSWITCH(xi)[l], _PSWITCH(c)[l], N);
    }
  }
  return(max_restarts*m);
}

static void _PSWITCH(init_lgcr)(const int _M, const int _V){
  static int Vo = -1;
  static int M = -1;

  int i;
  if((M != _M) || (_PSWITCH(lgcr_init) == 0) || (Vo != _V)){
    if(_PSWITCH(lgcr_init) == 1) _PSWITCH(free_lgcr)();
    Vo = _V;
    M = _M;
    _PSWITCH(a) = calloc(M+1, sizeof(_C_TYPE *));
    _PSWITCH(chi) = calloc(M, sizeof(_C_TYPE *));
    _PSWITCH(xi) = calloc(M, sizeof(_C_TYPE *));
    _PSWITCH(tmp) = calloc(Vo, sizeof(_C_TYPE));
    _PSWITCH(rho) = calloc(Vo, sizeof(_C_TYPE));
    _PSWITCH(_a) = calloc((M+1)*M, sizeof(_C_TYPE));
    _PSWITCH(a)[0] = _PSWITCH(_a);
    _PSWITCH(_chi) = calloc(M*Vo, sizeof(_C_TYPE));
    _PSWITCH(chi)[0] = _PSWITCH(_chi);
    _PSWITCH(_xi) = calloc(M*Vo, sizeof(_C_TYPE));
    _PSWITCH(xi)[0] = _PSWITCH(_xi);

    _PSWITCH(b) = calloc(M, sizeof(_F_TYPE));
    _PSWITCH(c) = calloc(M, sizeof(_C_TYPE));
    for(i = 1; i < M; i++) { 
      _PSWITCH(chi)[i] = _PSWITCH(chi)[i-1] + Vo;
      _PSWITCH(xi)[i] = _PSWITCH(xi)[i-1] + Vo;
      _PSWITCH(a)[i] = _PSWITCH(a)[i-1] + M;
    }
    _PSWITCH(a)[M] = _PSWITCH(a)[M-1] + M;
    _PSWITCH(lgcr_init) = 1;
  }
}

static void _PSWITCH(free_lgcr)() 
{
  _PSWITCH(lgcr_init) = 0;
  free(_PSWITCH(a));
  free(_PSWITCH(chi));
  free(_PSWITCH(_a));
  free(_PSWITCH(_chi));
  free(_PSWITCH(b));
  free(_PSWITCH(c));
  free(_PSWITCH(_xi));
  free(_PSWITCH(xi));
  free(_PSWITCH(rho));
  free(_PSWITCH(tmp));
  return;
}


void _PSWITCH(ldiff)(_C_TYPE * const Q, _C_TYPE * const R, _C_TYPE * const S, const int N) 
{
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N; ++i)
    Q[i] = R[i] - S[i];
  return;
}

void _PSWITCH(lassign)(_C_TYPE * const R, _C_TYPE * const S, const int N)
{
  memcpy(R, S, N*sizeof(_C_TYPE));
  return;
}

void _PSWITCH(ldiff_assign)(_C_TYPE * const Q, _C_TYPE * const S, const int N) 
{
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N; ++i)
    Q[i] -= S[i];
  return;
}

void _PSWITCH(ladd)(_C_TYPE * const Q, _C_TYPE * const R, _C_TYPE * const S, const int N) 
{
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N; ++i)
    Q[i] = R[i] + S[i];
  return;
}

void _PSWITCH(ladd_assign)(_C_TYPE * const Q, _C_TYPE * const S, const int N) 
{
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N; ++i)
    Q[i] += S[i];
  return;
}

_F_TYPE _PSWITCH(lsquare_norm)(_C_TYPE * const Q, const int N, const int parallel) 
{
  double nrm = 0.0;
  /* #ifdef OMP */
  /* #pragma omp parallel for */
  /* #endif */
  for(int i = 0; i < N; ++i)
    nrm += conj(Q[i]) * Q[i];

#ifdef MPI
  if(parallel)
    {
      double nrm2 = nrm;
      MPI_Allreduce(&nrm2, &nrm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif

  return((_F_TYPE)nrm);
}

_C_TYPE _PSWITCH(lscalar_prod)(_C_TYPE * const R, _C_TYPE * const S, const int N, const int parallel) 
{
  _Complex double res = 0.0;

  /* #ifdef OMP */
  /* #pragma omp parallel for */
  /* #endif */
  for(int i = 0; i < N; ++i)
    res += conj(R[i]) * S[i];

#ifdef MPI
  if(parallel)
    {
      _Complex double res2 = res;
      MPI_Allreduce(&res2, &res, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    }
#endif

  return((_C_TYPE)res);
}

_F_TYPE _PSWITCH(lscalar_prod_r)(_C_TYPE * const R, _C_TYPE * const S, const int N, const int parallel) 
{
  double res = 0.0;

  /* #ifdef OMP */
  /* #pragma omp parallel for */
  /* #endif */
  for(int i = 0; i < N; ++i) {
    res += creal(conj(R[i]) * S[i]);
  }

#ifdef MPI
  if(parallel) {
    double res2 = res;
    MPI_Allreduce(&res2, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif

  return((_F_TYPE)res);
}


void _PSWITCH(lmul_r)(_C_TYPE * const R, const _F_TYPE c, _C_TYPE * const S, const int N) 
{
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N; ++i)
    R[i] = c * S[i];
}

void _PSWITCH(lmul)(_C_TYPE * const R, const _C_TYPE c, _C_TYPE * const S, const int N) 
{
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N; ++i)
    R[i] = c * S[i];
}

void _PSWITCH(lassign_add_mul)(_C_TYPE * const R, _C_TYPE * const S, const _C_TYPE c, const int N)
{
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N; ++i)
    R[i] += c * S[i];
}

void _PSWITCH(lassign_add_mul_r)(_C_TYPE * const R, _C_TYPE * const S, const _F_TYPE c, const int N)
{
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N; ++i)
    R[i] += c * S[i];
}

void _PSWITCH(lassign_mul_add_r)(_C_TYPE * const R, const _F_TYPE c, _C_TYPE * const S, const int N)
{
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N; ++i)
    R[i] = c * R[i] + S[i];
}

void _PSWITCH(lassign_diff_mul)(_C_TYPE * const R, _C_TYPE * const S, const _C_TYPE c, const int N) 
{
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N; i++)
    R[i] -= c * S[i];
}
