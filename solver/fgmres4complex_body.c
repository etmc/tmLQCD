static void _PSWITCH(init_lgmres)(const int _M, const int _V);

static _Complex _F_TYPE ** _PSWITCH(H);
static _Complex _F_TYPE ** _PSWITCH(V);
static _Complex _F_TYPE ** _PSWITCH(Z);
static _Complex _F_TYPE * _PSWITCH(alpha);
static _Complex _F_TYPE * _PSWITCH(c);
static _F_TYPE * _PSWITCH(s);
extern void little_mg_precon_32(_Complex float *, _Complex float *);

int _PSWITCH(fgmres4complex)(_Complex _F_TYPE * const P, _Complex _F_TYPE * const Q,
			     const int m, const int max_restarts,
			     const double eps_sq, const int rel_prec,
			     const int N, const int parallel,
			     const int lda, const int precon, _PSWITCH(c_matrix_mult) f) {

  int restart, i, j, k;
  double beta, eps, norm;
  _Complex _F_TYPE tmp1, tmp2;
  _Complex _F_TYPE * r0;
  _Complex _F_TYPE ** solver_field = NULL;
  const int nr_sf = 3;

  _PSWITCH(init_lsolver_field)(&solver_field, /*why not N?*/ lda, nr_sf);/* #ifdef HAVE_LAPACK */

  eps=sqrt(eps_sq);
  _PSWITCH(init_lgmres)(m, lda);
  r0 = solver_field[0];
  
  norm = sqrt(_PSWITCH(lsquare_norm)(Q, N, parallel));

  _PSWITCH(lassign)(solver_field[2], P, N);
  for(restart = 0; restart < max_restarts; restart++){
    /* r_0=Q-AP  (b=Q, x+0=P) */
    f(r0, solver_field[2]);
    _PSWITCH(ldiff)(r0, Q, r0, N);

    /* v_0=r_0/||r_0|| */
    _PSWITCH(alpha)[0] = sqrt(_PSWITCH(lsquare_norm)(r0, N, parallel));

    if(g_proc_id == g_stdio_proc && g_debug_level > 2){
      printf("lFGMRES %d\t%g true residue\n", restart*m, creal(_PSWITCH(alpha)[0])*creal(_PSWITCH(alpha)[0]));
      fflush(stdout);
    }

    if(creal(_PSWITCH(alpha)[0])==0.){ 
      _PSWITCH(lassign)(P, solver_field[2], N);
      _PSWITCH(finalize_lsolver)(solver_field, nr_sf);
      return(restart*m);
    }

    _PSWITCH(lmul_r)(_PSWITCH(V)[0], 1./creal(_PSWITCH(alpha)[0]), r0, N);

    for(j = 0; j < m; j++){
      /* solver_field[0]=A*M^-1*v_j */

      if(precon == 0) {
	_PSWITCH(lassign)(_PSWITCH(Z)[j], _PSWITCH(V)[j], N);
      }
      else {
	_PSWITCH(little_mg_precon)(_PSWITCH(Z)[j], _PSWITCH(V)[j]);
      }

      f(r0, _PSWITCH(Z)[j]); 
      /* Set h_ij and omega_j */
      /* solver_field[1] <- omega_j */
      _PSWITCH(lassign)(solver_field[1], solver_field[0], N);
      for(i = 0; i <= j; i++){
	_PSWITCH(H)[i][j] = _PSWITCH(lscalar_prod)(_PSWITCH(V)[i], solver_field[1], N, parallel);
	_PSWITCH(lassign_diff_mul)(solver_field[1], _PSWITCH(V)[i], _PSWITCH(H)[i][j], N);
      }

      _PSWITCH(H)[j+1][j] = sqrt(_PSWITCH(lsquare_norm)(solver_field[1], N, parallel));
      for(i = 0; i < j; i++){
	tmp1 = _PSWITCH(H)[i][j];
	tmp2 = _PSWITCH(H)[i+1][j];
	(_PSWITCH(H)[i][j]) = (tmp2) * (_PSWITCH(s)[i]);
	(_PSWITCH(H)[i][j]) += conj(_PSWITCH(c)[i]) * (tmp1);
	(_PSWITCH(H)[i+1][j]) = (tmp1) * (_PSWITCH(s)[i]);
	(_PSWITCH(H)[i+1][j]) -= (_PSWITCH(c)[i]) * (tmp2);
      }

      /* Set beta, s, _PSWITCH(c), _PSWITCH(alpha)[j],[j+1] */
      beta = sqrt(creal(_PSWITCH(H)[j][j] * conj(_PSWITCH(H)[j][j])) + creal(_PSWITCH(H)[j+1][j] * conj(_PSWITCH(H)[j+1][j])));
      _PSWITCH(s)[j] = creal(_PSWITCH(H)[j+1][j]) / beta;
      (_PSWITCH(c)[j]) = (_PSWITCH(H)[j][j]) / beta;
      (_PSWITCH(H)[j][j]) = beta;
      (_PSWITCH(alpha)[j+1]) = (_PSWITCH(alpha)[j]) * (_PSWITCH(s)[j]);
      tmp1 = _PSWITCH(alpha)[j];
      (_PSWITCH(alpha)[j]) = conj(_PSWITCH(c)[j]) * (tmp1);

      /* precision reached? */
      if(g_proc_id == g_stdio_proc && g_debug_level > 2){
	printf("lFGMRES\t%d\t%g iterated residue\n", restart*m+j, creal(_PSWITCH(alpha)[j+1])*creal(_PSWITCH(alpha)[j+1]));
	fflush(stdout);
      }
      if(((creal(_PSWITCH(alpha)[j+1]) <= eps) && (rel_prec == 0)) || ((creal(_PSWITCH(alpha)[j+1]) <= eps*norm) && (rel_prec == 1))){
	(_PSWITCH(alpha)[j]) = (_PSWITCH(alpha)[j]) * (1./creal(_PSWITCH(H)[j][j]));
	_PSWITCH(lassign_add_mul)(solver_field[2], _PSWITCH(Z)[j], _PSWITCH(alpha)[j], N);
	for(i = j-1; i >= 0; i--){
	  for(k = i+1; k <= j; k++){
 	    (tmp1) = (_PSWITCH(H)[i][k]) * (_PSWITCH(alpha)[k]); 
	    (_PSWITCH(alpha)[i]) -= tmp1;
	  }
	  (_PSWITCH(alpha)[i]) = (_PSWITCH(alpha)[i]) * (1./creal(_PSWITCH(H)[i][i]));
	  _PSWITCH(lassign_add_mul)(solver_field[2], _PSWITCH(Z)[i], _PSWITCH(alpha)[i], N);
	}
	for(i = 0; i < m; i++){
	  _PSWITCH(alpha)[i] = creal(_PSWITCH(alpha)[i]);
	}
	_PSWITCH(lassign)(P, solver_field[2], N);
	_PSWITCH(finalize_lsolver)(solver_field, nr_sf);
	return(restart*m+j);
      }
      /* if not */
      else{
	if(j != m-1){
	  _PSWITCH(lmul_r)(_PSWITCH(V)[(j+1)], 1./creal(_PSWITCH(H)[j+1][j]), solver_field[1], N);
	}
      }

    }
    j=m-1;
    /* prepare for restart */
    (_PSWITCH(alpha)[j]) = (_PSWITCH(alpha)[j]) * (1./creal(_PSWITCH(H)[j][j]));
    _PSWITCH(lassign_add_mul)(solver_field[2], _PSWITCH(Z)[j], _PSWITCH(alpha)[j], N);
    for(i = j-1; i >= 0; i--){
      for(k = i+1; k <= j; k++){
	(tmp1) = (_PSWITCH(H)[i][k]) * (_PSWITCH(alpha)[k]);
	(_PSWITCH(alpha)[i]) -= tmp1;
      }
      (_PSWITCH(alpha)[i]) = (_PSWITCH(alpha)[i]) * (1./creal(_PSWITCH(H)[i][i]));
      _PSWITCH(lassign_add_mul)(solver_field[2], _PSWITCH(Z)[i], _PSWITCH(alpha)[i], N);
    }
    for(i = 0; i < m; i++){
      _PSWITCH(alpha)[i] = creal(_PSWITCH(alpha)[i]);
    }
  }

  /* If maximal number of restarts is reached */
  _PSWITCH(lassign)(P, solver_field[2], N);
  _PSWITCH(finalize_lsolver)(solver_field, nr_sf);
  return(-1);
}

static void _PSWITCH(init_lgmres)(const int _M, const int _V){
  static int Vo = -1;
  static int M = -1;
  static int init = 0;
  static _Complex _F_TYPE * _v;
  static _Complex _F_TYPE * _z;
  static _Complex _F_TYPE * _h;
  

  int i;
  if((M != _M)||(init == 0)||(Vo != _V)){
    if(init == 1){
      free(_PSWITCH(H));
      free(_PSWITCH(V));
      free(_h);
      free(_v);
      free(_PSWITCH(alpha));
      free(_PSWITCH(c));
      free(_PSWITCH(s));
    }
    Vo = _V;
    M = _M;
    _PSWITCH(H) = calloc(M+1, sizeof(_Complex _F_TYPE *));
    _PSWITCH(V) = calloc(M, sizeof(_Complex _F_TYPE *));
    _PSWITCH(Z) = calloc(M, sizeof(_Complex _F_TYPE *));
#if (defined SSE || defined SSE2)
    _h = calloc((M+2)*M, sizeof(_Complex _F_TYPE));
    _PSWITCH(H)[0] = (_Complex _F_TYPE *)(((unsigned long int)(_h)+ALIGN_BASE)&~ALIGN_BASE); 
    _v = calloc(M*Vo+1, sizeof(_Complex _F_TYPE));
    _PSWITCH(V)[0] = (_Complex _F_TYPE *)(((unsigned long int)(_v)+ALIGN_BASE)&~ALIGN_BASE);
    _z = calloc(M*Vo+1, sizeof(_Complex _F_TYPE));
    _PSWITCH(Z)[0] = (_Complex _F_TYPE *)(((unsigned long int)(_z)+ALIGN_BASE)&~ALIGN_BASE);
#else
    _h = calloc((M+1)*M, sizeof(_Complex _F_TYPE));
    _PSWITCH(H)[0] = _h;
    _v = calloc(M*Vo, sizeof(_Complex _F_TYPE));
    _PSWITCH(V)[0] = _v;
    _z = calloc(M*Vo, sizeof(_Complex _F_TYPE));
    _PSWITCH(Z)[0] = _z;
#endif
    _PSWITCH(s) = calloc(M, sizeof(_F_TYPE));
    _PSWITCH(c) = calloc(M, sizeof(_Complex _F_TYPE));
    _PSWITCH(alpha) = calloc(M+1, sizeof(_Complex _F_TYPE));
    for(i = 1; i < M; i++){
      _PSWITCH(V)[i] = _PSWITCH(V)[i-1] + Vo;
      _PSWITCH(H)[i] = _PSWITCH(H)[i-1] + M;
      _PSWITCH(Z)[i] = _PSWITCH(Z)[i-1] + Vo;
    }
    _PSWITCH(H)[M] = _PSWITCH(H)[M-1] + M;
    init = 1;
  }
  return;
}
