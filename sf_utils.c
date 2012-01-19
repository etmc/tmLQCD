/*******************************************************************************
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"

double calc_sq_plaq( void )
{
  int x;
  int x_p_mu;
  int x_p_nu;
  int mu;
  int nu;
  su3* u_mu_x;
  su3* u_nu_x_p_mu;
  su3* u_mu_x_p_nu;
  su3* u_nu_x;
  su3 tmp1;
  su3 tmp2;
  double tr;
  double sum = 0;

  for( x=0; x<VOLUME; x++ )
  for( mu=0; mu<4;  mu++ )
  for( nu=0; nu<mu; nu++ )
  {
    x_p_mu = g_iup[ x ][ mu ];
    x_p_nu = g_iup[ x ][ nu ];

    u_mu_x      = &g_gauge_field[ x      ][ mu ];
    u_nu_x_p_mu = &g_gauge_field[ x_p_mu ][ nu ];
    u_mu_x_p_nu = &g_gauge_field[ x_p_nu ][ mu ];
    u_nu_x      = &g_gauge_field[ x      ][ nu ];

    _su3_times_su3( tmp1, *u_nu_x, *u_mu_x_p_nu );
    _su3_times_su3( tmp2, *u_mu_x, *u_nu_x_p_mu );
    _trace_su3_times_su3d( tr, tmp1, tmp2 );

    sum += tr;
  }

  sum /= ((double)3);

  return sum;
}

double calc_bulk_sq_plaq( void )
{
  int x;
  int x_p_mu;
  int x_p_nu;
  int mu;
  int nu;
  su3* u_mu_x;
  su3* u_nu_x_p_mu;
  su3* u_mu_x_p_nu;
  su3* u_nu_x;
  su3 tmp1;
  su3 tmp2;
  double tr;
  double sum = 0;

  /* None of the forward plaquettes for t=T-1 contribute. */
  /* None of the forward plaquettes for t=0 contribute. */
  /* Also the space-time square for t=T-2 does not contribute. */

  for( x=0; x<VOLUME; x++ )
  for( mu=0; mu<4;  mu++ )
  for( nu=0; nu<mu; nu++ )
  if( ( g_t[ x ] != (T-1) ) &
      ( ( g_t[ x ] != (T-2) ) || ( ( mu != 3 ) && ( nu != 3 ) ) ) &
      ( g_t[ x ] != 0 ) )
  {
    x_p_mu = g_iup[ x ][ mu ];
    x_p_nu = g_iup[ x ][ nu ];

    u_mu_x      = &g_gauge_field[ x      ][ mu ];
    u_nu_x_p_mu = &g_gauge_field[ x_p_mu ][ nu ];
    u_mu_x_p_nu = &g_gauge_field[ x_p_nu ][ mu ];
    u_nu_x      = &g_gauge_field[ x      ][ nu ];

    _su3_times_su3( tmp1, *u_nu_x, *u_mu_x_p_nu );
    _su3_times_su3( tmp2, *u_mu_x, *u_nu_x_p_mu );
    _trace_su3_times_su3d( tr, tmp1, tmp2 );

    sum += tr;
  }

  sum /= ((double)3);

  return sum;
}

double calc_boundary_space_space_sq_plaq( void )
{
  int x;
  int x_p_mu;
  int x_p_nu;
  int mu;
  int nu;
  su3* u_mu_x;
  su3* u_nu_x_p_mu;
  su3* u_mu_x_p_nu;
  su3* u_nu_x;
  su3 tmp1;
  su3 tmp2;
  double tr;
  double sum = 0;

  /* We need the space-space plaquettes for t=T-1 and t=0. */

  for( x=0; x<VOLUME; x++ )
  for( mu=0; mu<4;  mu++ )
  for( nu=0; nu<mu; nu++ )
  if( ( ( g_t[ x ] == (T-1) ) & ( mu != 3 ) & ( nu != 3 ) ) ||
      ( ( g_t[ x ] ==     0 ) & ( mu != 3 ) & ( nu != 3 ) ) )
  {
    x_p_mu = g_iup[ x ][ mu ];
    x_p_nu = g_iup[ x ][ nu ];

    u_mu_x      = &g_gauge_field[ x      ][ mu ];
    u_nu_x_p_mu = &g_gauge_field[ x_p_mu ][ nu ];
    u_mu_x_p_nu = &g_gauge_field[ x_p_nu ][ mu ];
    u_nu_x      = &g_gauge_field[ x      ][ nu ];

    _su3_times_su3( tmp1, *u_nu_x, *u_mu_x_p_nu );
    _su3_times_su3( tmp2, *u_mu_x, *u_nu_x_p_mu );
    _trace_su3_times_su3d( tr, tmp1, tmp2 );

    sum += tr;
  }

  sum /= ((double)3);

  return sum;
}

double calc_boundary_space_time_sq_plaq( void )
{
  int x;
  int x_p_mu;
  int x_p_nu;
  int mu;
  int nu;
  su3* u_mu_x;
  su3* u_nu_x_p_mu;
  su3* u_mu_x_p_nu;
  su3* u_nu_x;
  su3 tmp1;
  su3 tmp2;
  double tr;
  double sum = 0;

  /* The space-time square for t=T-2 contributes. */
  /* The space-time square for t=0 contributes. */

  for( x=0; x<VOLUME; x++ )
  for( mu=0; mu<4;  mu++ )
  for( nu=0; nu<mu; nu++ )
  if( ( ( g_t[ x ] == (T-2) ) & ( ( mu == 3 ) || ( nu == 3 ) ) ) ||
      ( ( g_t[ x ] ==     0 ) & ( ( mu == 3 ) || ( nu == 3 ) ) ) )
  {
    x_p_mu = g_iup[ x ][ mu ];
    x_p_nu = g_iup[ x ][ nu ];

    u_mu_x      = &g_gauge_field[ x      ][ mu ];
    u_nu_x_p_mu = &g_gauge_field[ x_p_mu ][ nu ];
    u_mu_x_p_nu = &g_gauge_field[ x_p_nu ][ mu ];
    u_nu_x      = &g_gauge_field[ x      ][ nu ];

    _su3_times_su3( tmp1, *u_nu_x, *u_mu_x_p_nu );
    _su3_times_su3( tmp2, *u_mu_x, *u_nu_x_p_mu );
    _trace_su3_times_su3d( tr, tmp1, tmp2 );

    sum += tr;
  }

  sum /= ((double)3);

  return sum;
}

double calc_wrapped_sq_plaq( void )
{
  int x;
  int x_p_mu;
  int x_p_nu;
  int mu;
  int nu;
  su3* u_mu_x;
  su3* u_nu_x_p_mu;
  su3* u_mu_x_p_nu;
  su3* u_nu_x;
  su3 tmp1;
  su3 tmp2;
  double tr;
  double sum = 0;

  for( x=0; x<VOLUME; x++ )
  for( mu=0; mu<4;  mu++ )
  for( nu=0; nu<mu; nu++ )
  if( ( g_t[ x ] == (T-1) ) & ( ( mu == 3 ) || ( nu == 3 ) ) )
  {
    x_p_mu = g_iup[ x ][ mu ];
    x_p_nu = g_iup[ x ][ nu ];

    u_mu_x      = &g_gauge_field[ x      ][ mu ];
    u_nu_x_p_mu = &g_gauge_field[ x_p_mu ][ nu ];
    u_mu_x_p_nu = &g_gauge_field[ x_p_nu ][ mu ];
    u_nu_x      = &g_gauge_field[ x      ][ nu ];

    _su3_times_su3( tmp1, *u_nu_x, *u_mu_x_p_nu );
    _su3_times_su3( tmp2, *u_mu_x, *u_nu_x_p_mu );
    _trace_su3_times_su3d( tr, tmp1, tmp2 );

    sum += tr;
  }

  sum /= ((double)3);

  return sum;
}

double calc_rect_plaq( void )
{
  int x;
  int x_p_mu;
  int x_p_nu;
  int x_p_nu_m_mu;
  int x_m_mu;
  int mu;
  int nu;
  su3* u_mu_x;
  su3* u_nu_x_p_mu;
  su3* u_mu_x_p_nu;
  su3* u_mu_x_p_nu_m_mu;
  su3* u_nu_x_m_mu;
  su3* u_mu_x_m_mu;
  su3 tmp1;
  su3 tmp2;
  su3 tmp3;
  double tr;
  double sum = 0;

  for( x=0; x<VOLUME; x++ )
  for( mu=0; mu<4;  mu++ )
  for( nu=0; nu<4; nu++ )
  if( mu != nu )
  {
    x_p_mu      = g_iup[ x      ][ mu ];
    x_p_nu      = g_iup[ x      ][ nu ];
    x_p_nu_m_mu = g_idn[ x_p_nu ][ mu ];
    x_m_mu      = g_idn[ x      ][ mu ];

    u_mu_x           = &g_gauge_field[ x           ][ mu ];
    u_nu_x_p_mu      = &g_gauge_field[ x_p_mu      ][ nu ];
    u_mu_x_p_nu      = &g_gauge_field[ x_p_nu      ][ mu ];
    u_mu_x_p_nu_m_mu = &g_gauge_field[ x_p_nu_m_mu ][ mu ];
    u_nu_x_m_mu      = &g_gauge_field[ x_m_mu      ][ nu ];
    u_mu_x_m_mu      = &g_gauge_field[ x_m_mu      ][ mu ];

    _su3_times_su3( tmp1, *u_mu_x_p_nu_m_mu, *u_mu_x_p_nu );
    _su3_times_su3( tmp2, *u_nu_x_m_mu, tmp1 );

    _su3_times_su3( tmp1, *u_mu_x, *u_nu_x_p_mu );
    _su3_times_su3( tmp3, *u_mu_x_m_mu, tmp1 );

    _trace_su3_times_su3d( tr, tmp2, tmp3 );

    sum += tr;
  }

  sum /= ((double)3);

  return sum;
}
