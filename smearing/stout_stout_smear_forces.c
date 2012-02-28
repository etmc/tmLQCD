#include "stout.ih"

void stout_smear_forces(stout_control *control, adjoint_field_t in)
{
  if (!control->calculate_force_terms)
    fatal_error("Stout control structure not setup for calculating force terms.", "stout_smear_forces");

  if (!control->smearing_performed)
    fatal_error("Stout smearing not yet performed.", "stout_smear_forces");
  
  for (unsigned int ctr = 0; ctr < 4; ++ctr)
    control->scratch[ctr] = get_gauge_field();
  
  /* We'll need the forces in their tangent space representation, so let's first build this up. */
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
      _make_su3(control->scratch[0].field[x][mu], in.field[x][mu]);
    }
  
  _Complex double trace;
  /* The modifications are done backwards, all the time peeling off one layer of stouting... */
  for (int iter = control->iterations - 1; iter > 0; --iter)
  {
    /* NOTE Write routines operating on whole gauge fields? This is clunky and tedious. */
    /* NOTE Or the opposite for cache coherence... */
    
    /* Calculate Lambda as defined by Peardon and Morningstar */
    for (unsigned int x = 0; x < VOLUME; ++x)
      for (unsigned int mu = 0; mu < 4; ++mu)
      {
        /* scratch[3] == Tr[Sigma' * B2 * U] * Q */
        _su3_times_su3(control->scratch[2].field[x][mu], control->scratch[0].field[x][mu], control->B2[iter].field[x][mu]);
        _trace_su3_times_su3(trace, control->scratch[2].field[x][mu], control->U[iter].field[x][mu]); /* NOTE Iteration of U before last! */
        _complex_times_su3(control->scratch[3].field[x][mu], trace, control->Q[iter].field[x][mu]);
        
        /* scratch[1] == ( Tr[Sigma' * B1 * U] + Tr[Sigma' * B2 * U] * Q ) * Q */
        _su3_times_su3(control->scratch[2].field[x][mu], control->scratch[0].field[x][mu], control->B1[iter].field[x][mu]);
        _trace_su3_times_su3(trace, control->scratch[2].field[x][mu], control->U[iter].field[x][mu]); /* NOTE Iteration of U before last! */
        _su3_add_equals_complex_identity(control->scratch[3].field[x][mu], trace)
        _su3_times_su3(control->scratch[1].field[x][mu], control->scratch[3].field[x][mu], control->Q[iter].field[x][mu]);

        /* scratch[1] == ( Tr[Sigma' * B1 * U] + Tr[Sigma' * B2 * U] * Q ) * Q  + f_1 * U * Sigma' */
        _su3_times_su3(control->scratch[2].field[x][mu], control->U[iter].field[x][mu], control->scratch[0].field[x][mu]);
        _complex_times_su3(control->scratch[3].field[x][mu], control->f1[iter][x][mu], control->scratch[2].field[x][mu]);
        _su3_plus_su3(control->scratch[1].field[x][mu], control->scratch[3].field[x][mu], control->scratch[1].field[x][mu]);
        
        /* scratch[1] == ( Tr[Sigma' * B1 * U] + Tr[Sigma' * B2 * U] * Q ) * Q  + f_1 * U * Sigma' + f2 * Q * U * Sigma' */
        _su3_times_su3(control->scratch[3].field[x][mu], control->Q[iter].field[x][mu], control->scratch[2].field[x][mu]);
        _complex_times_su3(control->scratch[3].field[x][mu], control->f2[iter][x][mu], control->scratch[3].field[x][mu]);
        _su3_plus_su3(control->scratch[1].field[x][mu], control->scratch[3].field[x][mu], control->scratch[1].field[x][mu]);
        
        /* scratch[1] == ( Tr[Sigma' * B1 * U] + Tr[Sigma' * B2 * U] * Q ) * Q  + f_1 * U * Sigma' + f2 * Q * U * Sigma' + f2 * U * Sigma' * Q* == Gamma */
        _su3_times_su3(control->scratch[3].field[x][mu], control->scratch[2].field[x][mu], control->Q[iter].field[x][mu]);
        _complex_times_su3(control->scratch[3].field[x][mu], control->f2[iter][x][mu], control->scratch[3].field[x][mu]);
        _su3_plus_su3(control->scratch[1].field[x][mu], control->scratch[3].field[x][mu], control->scratch[1].field[x][mu]);
        
        /* scratch[1] = 0.5 * ((Gamma + Gamma^dag) - Tr(Gamma + Gamma^dag) / N_c) == Lambda */
        project_herm(&control->scratch[1].field[x][mu]);

        /* scratch[3] = Sigma' * U' * U^dag = Sigma' * exp[i Q] * U * U^dag = Sigma' * exp[i Q] */
        /* NOTE Check indices on U from here on downwards -- iterations look a little fishy. */
        _su3_times_su3d(control->scratch[3].field[x][mu], control->U[iter + 1].field[x][mu], control->U[iter].field[x][mu]);
        _su3_times_su3(control->scratch[2].field[x][mu], control->scratch[0].field[x][mu], control->scratch[3].field[x][mu]);
      }
        
    /********************************************************************************************************/
    /* At this point, scratch[0] contains Sigma', scratch[1] contains Lambda, scratch[2] contains exp[i Q]. */
    /********************************************************************************************************/
    
    _Complex double minus_i_rho = -I * control->rho; /* NOTE ... or fold this into the macros below? */
    generic_staples(control->scratch[4], control->U[iter]);
    
    for (unsigned int x = 0; x < VOLUME; ++x)
      for (unsigned int mu = 0; mu < 4; ++mu)
      {
        /* scratch[3] = i * C^dag * Lambda */
        _su3d_times_su3_acc(control->scratch[3].field[x][mu], control->scratch[4].field[x][mu], control->scratch[1].field[x][mu]);
        _complex_times_su3(control->scratch[3].field[x][mu], I, control->scratch[3].field[x][mu]);
      }  

    su3_tuple aggr;
    su3 t1;
    su3 t2;
      
#define _STAPLE_DERIV_TERM_1(mu, nu) \
  _su3_times_su3d(t1, control->U[iter].field[x][nu], control->U[iter].field[g_iup[x][nu]][mu]); \
  _su3_times_su3d(t2, t1, control->U[iter].field[x][nu]); \
  _su3_times_su3_acc(aggr[mu], t2, control->scratch[1].field[x][nu]);

/* NOTE The first term below wants to multiply two daggered matrices, which we don't have code for. Hence the rewrite */
#define _STAPLE_DERIV_TERM_2(mu, nu) \
  _su3_times_su3(t1, control->U[iter].field[g_idn[x][nu]][mu], control->U[iter].field[g_idn[g_iup[x][mu]][nu]][nu]); \
  _su3d_times_su3(t2, t1, control->scratch[1].field[g_idn[x][nu]][mu]); \
  _su3_times_su3_acc(aggr[mu], t2, control->U[iter].field[g_idn[x][nu]][nu]);

#define _STAPLE_DERIV_TERM_3(mu, nu) \
  _su3d_times_su3(t1, control->U[iter].field[g_idn[g_iup[x][mu]][nu]][nu], control->scratch[1].field[g_idn[g_iup[x][mu]][nu]][nu]); \
  _su3_times_su3d(t2, t1, control->U[iter].field[g_idn[x][nu]][mu]); \
  _su3_times_su3_acc(aggr[mu], t2, control->U[iter].field[g_idn[x][nu]][nu]);

/* NOTE Derivative terms 2 and 4 could be easily combined by plugging in a commutator on Lambda. */

#define _STAPLE_DERIV_TERM_4(mu, nu) \
  _su3_times_su3(t1, control->U[iter].field[g_idn[x][nu]][mu], control->U[iter].field[g_idn[g_iup[x][mu]][nu]][nu]); \
  _su3d_times_su3(t2, t1, control->scratch[1].field[g_idn[x][nu]][nu]); \
  _su3_times_su3(t1, t2, control->U[iter].field[g_idn[x][nu]][nu]); \
  _su3_refac_acc(aggr[mu], -1.0, t1);

#define _STAPLE_DERIV_TERM_5(mu, nu) \
  _su3_times_su3(t1, control->scratch[1].field[g_iup[x][mu]][nu], control->U[iter].field[g_iup[x][mu]][nu]); \
  _su3_times_su3d(t2, t1, control->U[iter].field[g_iup[x][nu]][mu]); \
  _su3_times_su3d(aggr[mu], t2, control->U[iter].field[x][nu]); \
  _su3_refac_acc(aggr[mu], -1.0, t1);

#define _STAPLE_DERIV_TERM_6(mu, nu) \
  _su3_times_su3d(t1, control->U[iter].field[g_iup[x][mu]][nu], control->U[iter].field[g_iup[x][nu]][mu]); \
  _su3_times_su3(t2, t1, control->scratch[1].field[g_iup[x][nu]][mu]); \
  _su3_times_su3d_acc(aggr[mu], t2, control->U[iter].field[x][nu]);
  
    for (unsigned int x = 0; x < VOLUME; ++x)
    {
      for (unsigned int mu = 0; mu < 4; ++mu)
      {
        _su3_zero(aggr[mu]);
      }
      
      /* NOTE The stuff below can obviously be done better, with less repetition. */
      
      _STAPLE_DERIV_TERM_1(0, 1);
      _STAPLE_DERIV_TERM_1(0, 2);
      _STAPLE_DERIV_TERM_1(0, 3);
      
      _STAPLE_DERIV_TERM_2(0, 1);
      _STAPLE_DERIV_TERM_2(0, 2);
      _STAPLE_DERIV_TERM_2(0, 3);

      _STAPLE_DERIV_TERM_3(0, 1);
      _STAPLE_DERIV_TERM_3(0, 2);
      _STAPLE_DERIV_TERM_3(0, 3);
      
      _STAPLE_DERIV_TERM_4(0, 1);
      _STAPLE_DERIV_TERM_4(0, 2);
      _STAPLE_DERIV_TERM_4(0, 3);
      
      _STAPLE_DERIV_TERM_5(0, 1);
      _STAPLE_DERIV_TERM_5(0, 2);
      _STAPLE_DERIV_TERM_5(0, 3);
      
      _STAPLE_DERIV_TERM_6(0, 1);
      _STAPLE_DERIV_TERM_6(0, 2);
      _STAPLE_DERIV_TERM_6(0, 3);
      
      
      _STAPLE_DERIV_TERM_1(1, 0);
      _STAPLE_DERIV_TERM_1(1, 2);
      _STAPLE_DERIV_TERM_1(1, 3);
      
      _STAPLE_DERIV_TERM_2(1, 0);
      _STAPLE_DERIV_TERM_2(1, 2);
      _STAPLE_DERIV_TERM_2(1, 3);

      _STAPLE_DERIV_TERM_3(1, 0);
      _STAPLE_DERIV_TERM_3(1, 2);
      _STAPLE_DERIV_TERM_3(1, 3);
      
      _STAPLE_DERIV_TERM_4(1, 0);
      _STAPLE_DERIV_TERM_4(1, 2);
      _STAPLE_DERIV_TERM_4(1, 3);
      
      _STAPLE_DERIV_TERM_5(1, 0);
      _STAPLE_DERIV_TERM_5(1, 2);
      _STAPLE_DERIV_TERM_5(1, 3);
      
      _STAPLE_DERIV_TERM_6(1, 0);
      _STAPLE_DERIV_TERM_6(1, 2);
      _STAPLE_DERIV_TERM_6(1, 3);

      
      _STAPLE_DERIV_TERM_1(2, 0);
      _STAPLE_DERIV_TERM_1(2, 1);
      _STAPLE_DERIV_TERM_1(2, 3);
      
      _STAPLE_DERIV_TERM_2(2, 0);
      _STAPLE_DERIV_TERM_2(2, 1);
      _STAPLE_DERIV_TERM_2(2, 3);

      _STAPLE_DERIV_TERM_3(2, 0);
      _STAPLE_DERIV_TERM_3(2, 1);
      _STAPLE_DERIV_TERM_3(2, 3);
      
      _STAPLE_DERIV_TERM_4(2, 0);
      _STAPLE_DERIV_TERM_4(2, 1);
      _STAPLE_DERIV_TERM_4(2, 3);
      
      _STAPLE_DERIV_TERM_5(2, 0);
      _STAPLE_DERIV_TERM_5(2, 1);
      _STAPLE_DERIV_TERM_5(2, 3);
      
      _STAPLE_DERIV_TERM_6(2, 0);
      _STAPLE_DERIV_TERM_6(2, 1);
      _STAPLE_DERIV_TERM_6(2, 3);
      
      
      _STAPLE_DERIV_TERM_1(3, 0);
      _STAPLE_DERIV_TERM_1(3, 1);
      _STAPLE_DERIV_TERM_1(3, 2);
      
      _STAPLE_DERIV_TERM_2(3, 0);
      _STAPLE_DERIV_TERM_2(3, 1);
      _STAPLE_DERIV_TERM_2(3, 2);

      _STAPLE_DERIV_TERM_3(3, 0);
      _STAPLE_DERIV_TERM_3(3, 1);
      _STAPLE_DERIV_TERM_3(3, 2);
      
      _STAPLE_DERIV_TERM_4(3, 0);
      _STAPLE_DERIV_TERM_4(3, 1);
      _STAPLE_DERIV_TERM_4(3, 2);
      
      _STAPLE_DERIV_TERM_5(3, 0);
      _STAPLE_DERIV_TERM_5(3, 1);
      _STAPLE_DERIV_TERM_5(3, 2);
      
      _STAPLE_DERIV_TERM_6(3, 0);
      _STAPLE_DERIV_TERM_6(3, 1);
      _STAPLE_DERIV_TERM_6(3, 2);
      
      for (unsigned int mu = 0; mu < 4; ++mu)
      {
	_su3_times_su3(t1, control->scratch[0].field[x][mu], control->scratch[2].field[x][mu]);
	_su3_plus_su3(control->scratch[0].field[x][mu], t1, control->scratch[3].field[x][mu]);
	/* NOTE _su3_refac_acc is abused slightly by providing a complex value here... */
	_su3_refac_acc(control->scratch[0].field[x][mu], minus_i_rho, aggr[mu]);
      }
    }
    
#undef _STAPLE_DERIV_TERM_1
#undef _STAPLE_DERIV_TERM_2
#undef _STAPLE_DERIV_TERM_3
#undef _STAPLE_DERIV_TERM_4
#undef _STAPLE_DERIV_TERM_5
#undef _STAPLE_DERIV_TERM_6

  }
      
  /* The force terms are still in the tangent space representation, so project them back to the adjoint one */
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
      _trace_lambda(control->force_result.field[x][mu], control->scratch[0].field[x][mu]);
    }
  
  for (unsigned int ctr = 0; ctr < 4; ++ctr)
    return_gauge_field(&control->scratch[ctr]);
}
