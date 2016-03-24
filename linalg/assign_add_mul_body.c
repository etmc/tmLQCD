void _PSWITCH(assign_add_mul)(_PTSWITCH(spinor) * const R, _PTSWITCH(spinor) * const S, 
			      const _C_TYPE c, const int N)
{
#ifdef OMP
#pragma omp parallel
  {
#endif
    _PTSWITCH(spinor) *r,*s;

#ifdef OMP
#pragma omp for
#endif
  for (int ix=0; ix<N; ix++)
  {
    r=(_PTSWITCH(spinor) *) R + ix;
    s=(_PTSWITCH(spinor) *) S + ix;

    r->s0.c0 += c * s->s0.c0;
    r->s0.c1 += c * s->s0.c1;
    r->s0.c2 += c * s->s0.c2;

    r->s1.c0 += c * s->s1.c0;
    r->s1.c1 += c * s->s1.c1;
    r->s1.c2 += c * s->s1.c2;

    r->s2.c0 += c * s->s2.c0;
    r->s2.c1 += c * s->s2.c1;
    r->s2.c2 += c * s->s2.c2;

    r->s3.c0 += c * s->s3.c0;
    r->s3.c1 += c * s->s3.c1;
    r->s3.c2 += c * s->s3.c2;
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

void _PSWITCH(assign_add_mul_ts)(_PTSWITCH(spinor) * const R, _PTSWITCH(spinor) * const S, 
			      const _C_TYPE c, const int N)
{
  _PTSWITCH(spinor) *r,*s;

  for (int ix=0; ix<N; ix++)
  {
    r=(_PTSWITCH(spinor) *) R + ix;
    s=(_PTSWITCH(spinor) *) S + ix;

    r->s0.c0 += c * s->s0.c0;
    r->s0.c1 += c * s->s0.c1;
    r->s0.c2 += c * s->s0.c2;

    r->s1.c0 += c * s->s1.c0;
    r->s1.c1 += c * s->s1.c1;
    r->s1.c2 += c * s->s1.c2;

    r->s2.c0 += c * s->s2.c0;
    r->s2.c1 += c * s->s2.c1;
    r->s2.c2 += c * s->s2.c2;

    r->s3.c0 += c * s->s3.c0;
    r->s3.c1 += c * s->s3.c1;
    r->s3.c2 += c * s->s3.c2;
  }

  return;
}
