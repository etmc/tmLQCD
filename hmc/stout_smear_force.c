/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2007,2008 Jan Volkholz
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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <complex.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "stout_smear.h"
#include "get_staples.h"
#include "read_input.h"

#include "lime.h"

/*void stout_smear_force(int x, int mu, su3 force){*/
void stout_smear_force()
{
  su3adj before_force;

  extern su3 ** g_gauge_field_smeared;
  extern su3 ** g_stout_force_field;
  extern su3 ** g_previous_stout_force_field;
  extern su3 ** g_C_smearing;
  extern su3 ** g_Q_smearing;
  extern su3 ** g_Q_squared_smearing;
  extern su3 ** g_B1_smearing;
  extern su3 ** g_B2_smearing;
  extern su3 ** g_Gamma_smearing;
  extern su3 ** g_Lambda_smearing;

  extern double * g_c0_smearing;
  extern double * g_c1_smearing;

  extern complex * g_f0_smearing;
  extern complex * g_f1_smearing;
  extern complex * g_f2_smearing;

  extern complex * g_b10_smearing;
  extern complex * g_b11_smearing;
  extern complex * g_b12_smearing;
  extern complex * g_b20_smearing;
  extern complex * g_b21_smearing;
  extern complex * g_b22_smearing;

  extern complex * g_r10_smearing;
  extern complex * g_r11_smearing;
  extern complex * g_r12_smearing;
  extern complex * g_r20_smearing;
  extern complex * g_r21_smearing;
  extern complex * g_r22_smearing;

  const int dim = 4;

  int c0_is_negative;

  int x, mu, nu, stout_iter_counter;
  double c0_magnitude, theta, eps, sqrt_two, u, u_squared, w, w_squared, xi_0, xi_1;
  double cos_u, sin_u, cos_w, sin_w, sin_2u, cos_2u, u_cos_u, u_sin_u, u_cos_2u, u_sin_2u, denominator;

  double tmp_d_0, tmp_d_1, tmp_d_2, c0_max;
  complex tmp_c_0, tmp_c_1;
  su3 tmp_su3_0, tmp_su3_1, tmp_su3_2, tmp_su3_3;

  printf("Running stout_smear_force()\n");
  before_force.d1 = df0[0][0].d1;
  before_force.d2 = df0[0][0].d2;
  before_force.d3 = df0[0][0].d3;
  before_force.d4 = df0[0][0].d4;
  before_force.d5 = df0[0][0].d5;
  before_force.d6 = df0[0][0].d6;
  before_force.d7 = df0[0][0].d7;
  before_force.d8 = df0[0][0].d8;

  /*###########################################################################*/

  /*read_lime_gauge_field("config_saved.dat");
    printf("DUMP OF g_gauge_field\n");
    print_config_to_screen(g_gauge_field);

    load_config_from_file(g_stout_force_field, "force_saved.dat");
    printf("DUMP OF g_stout_force_field\n");
    print_config_to_screen(g_stout_force_field);*/

  /*###########################################################################*/

  /*
   *  first we copy the gauge field into a new variable
   */
  for(x = 0; x < VOLUME; x++)
    for(mu = 0; mu < dim; mu++)
    {
      _su3_assign(g_gauge_field_smeared[x][mu], g_gauge_field[x][mu]);
    }

  /*for(stout_iter_counter = stout_no_iter-1; stout_iter_counter > -1; stout_iter_counter--)*/
  for(stout_iter_counter = 0; stout_iter_counter < stout_no_iter; stout_iter_counter++)
  {
    for(mu = 0; mu < dim; mu++)
    {
      for(x = 0; x < VOLUME; x++)
      {
        printf("x=%d  mu=%d\n", x, mu);
        print_su3(&(g_gauge_field_smeared[x][mu]));
        /*
         *  first we save the original force field, so we have \Sigma'
         *  available when iterating thru  eqtn(75) hep-lat/0311018
         */
        _su3_assign(g_previous_stout_force_field[x][mu], g_stout_force_field[x][mu]);

        /*
         *  get C from eqtn(1) hep-lat/0311018
         */
        g_C_smearing[x][mu] = get_staples(x, mu, g_gauge_field_smeared);
        scale_su3(&(g_C_smearing[x][mu]), stout_rho);

        /*printf("X=%d  MU=%d\n", x, mu);
          print_su3(&(g_C_smearing[x][mu]));*/

        /*
         *  get Q from eqtn(1) hep-lat/0311018
         */
        _su3_times_su3d(g_Q_smearing[x][mu], g_C_smearing[x][mu], g_gauge_field_smeared[x][mu]);
        project_anti_herm(&(g_Q_smearing[x][mu])) ;

        /*
         *  we save the square of the above in order to
         *  some CPU effort
         */
        _su3_times_su3(g_Q_squared_smearing[x][mu], g_Q_smearing[x][mu], g_Q_smearing[x][mu]);

        /*printf("X=%d  MU=%d\n", x, mu);
          print_su3(&(g_Q_smearing[x][mu]));*/

        /*
         * c0 from eqtn (14) hep-lat/0311018
         *
         * c0 = (Tr Q^3)/3
         */
        _trace_su3_times_su3(tmp_c_0, g_Q_smearing[x][mu], g_Q_squared_smearing[x][mu]);
        g_c0_smearing[x] = -0.33333333333333333333333333333333333333333 * tmp_c_0.im;

        /*printf("\n");
          printf("site = %d   mu = %d\n", x, mu);
          printf("trQQQ = %12.14lf + i * %12.14lf    c0 = %12.14lf\n", tmp_c_0.re, tmp_c_0.im, g_c0_smearing[x]);*/

        /*
         * c1 from eqtn (15) hep-lat/0311018
         *
         * c1 = (Tr Q^2)/2
         */
        g_c1_smearing[x] = -0.5 * (g_Q_squared_smearing[x][mu].c00.re + g_Q_squared_smearing[x][mu].c11.re + g_Q_squared_smearing[x][mu].c22.re);

        /*printf("trQQ  = %12.14lf + i * %12.14lf     c1 = %12.14lf\n", g_Q_squared_smearing[x][mu].c00.re + g_Q_squared_smearing[x][mu].c11.re + g_Q_squared_smearing[x][mu].c22.re, g_Q_squared_smearing[x][mu].c00.im + g_Q_squared_smearing[x][mu].c11.im + g_Q_squared_smearing[x][mu].c22.im,  g_c1_smearing[x]);*/

        /*
         *  as stated in the chroma sources ("stout_utils.cc")
         *  we need to watch out for the case c0 \approx 0
         *  this case is ommitted in hep-lat/0311018
         */
        if (g_c1_smearing[x] < 4.0e-3)
        {
          g_f0_smearing[x].re = 1-g_c0_smearing[x] * g_c0_smearing[x] / 720.0;
          g_f0_smearing[x].im =  -(g_c0_smearing[x] / 6.0) * (1.0 - (g_c1_smearing[x] / 20)
              * ( 1.0 - (g_c1_smearing[x] / 42.0)));

          g_f1_smearing[x].re =  g_c0_smearing[x] / 24.0 * (1.0 - g_c1_smearing[x] / 15.0 * (1.0 - 3.0 * g_c1_smearing[x] / 112.0)) ;
          g_f1_smearing[x].im =  1.0 - g_c1_smearing[x] / 6.0 * (1.0 - g_c1_smearing[x] / 20.0 * (1.0 - g_c1_smearing[x] / 42.0))
            - g_c0_smearing[x] * g_c0_smearing[x] / 5040.0;

          g_f2_smearing[x].re = 0.5 * (-1.0 + g_c1_smearing[x] / 12.0 * (1.0 - g_c1_smearing[x] / 30.0
                * (1.0 - g_c1_smearing[x] / 56.0)) + g_c0_smearing[x] * g_c0_smearing[x] / 20160.0);
          g_f2_smearing[x].im = 0.5 * (g_c0_smearing[x] / 60.0 * (1.0 - g_c1_smearing[x] / 21.0 * (1.0 - g_c1_smearing[x] / 48.0)));

          g_b20_smearing[x].re = -g_c0_smearing[x] / 360.0;
          g_b20_smearing[x].im =  -(1.0 / 6.0) * (1.0 - (g_c1_smearing[x] / 20.0) * (1.0 - g_c1_smearing[x] / 42.0));

          g_b10_smearing[x].re = 0.0;
          g_b10_smearing[x].im = (g_c0_smearing[x] / 120.0)*(1.0 - g_c1_smearing[x] / 21.0);

          g_b21_smearing[x].re = (1.0 / 24.0) * (1.0 - g_c1_smearing[x] / 15.0 * (1.0 - 3.0 * g_c1_smearing[x] / 112.0));
          g_b21_smearing[x].im = -g_c0_smearing[x] / 2520.0;

          g_b11_smearing[x].re = -g_c0_smearing[x] / 360.0 * (1.0 - 3.0 * g_c1_smearing[x] / 56.0);
          g_b11_smearing[x].im = -1.0 / 6.0 * (1.0 - g_c1_smearing[x] / 10.0 * (1.0 - g_c1_smearing[x] / 28.0));

          g_b22_smearing[x].re = 0.5 * g_c0_smearing[x] / 10080.0;
          g_b22_smearing[x].im = 0.5 * (1.0 / 60.0 * (1.0 - g_c1_smearing[x] / 21.0 * (1.0 - g_c1_smearing[x] / 48.0)) );

          g_b12_smearing[x].re = 0.5 * (1.0 / 12.0 * (1.0 - (2.0 * g_c1_smearing[x] / 30.0) * (1.0 - 3.0*g_c1_smearing[x]/112.0)) );
          g_b12_smearing[x].im = 0.5*( -g_c0_smearing[x]/1260.0 * (1.0-g_c1_smearing[x]/24.0) );
        }

        /*
         *  theta from eqtn (25) hep-lat/0311018
         */

        /*
         *  this is the case as discussed in the paper
         */
        else
        {
          c0_magnitude = fabs(g_c0_smearing[x]);

          if(g_c0_smearing[x] < 0)
            c0_is_negative = 1;
          else
            c0_is_negative = 0;
          /*
           * c0_max from eqtn (17) hep-lat/0311018
           *
           * c0_max = 2*(c1/3)^1.5
           *
           * 2/(3)^1.5 = 0.384900179459750...
           */
          c0_max = 0.38490017945975058411969329005490109652452737968 
            * sqrt(g_c1_smearing[x]) * g_c1_smearing[x];

          /*printf("c0_magnitude = %12.14lf c0_max = %12.14lf\n", c0_magnitude, c0_max);*/
          /*
           *  CORNER CASE 1:
           *  now we need to consoider a cornercase not
           *  discussed in the paper but in chroma's stout_utils.cc
           */
          eps =  (c0_max - c0_magnitude)/c0_max;

          /*
           *  if by rounding error c0_magnitude > c0_max
           *  they are actaully equal, meaning theta=0
           */
          if(eps < 0)
          {
            theta=0;
          }
          else

            /*
             *  CORNER CASE 2:
             *  now c0 \approx  c0_max, meaning
             *  c0_magnitude / c0_max approx 1
             *  therefore this is handled best by a
             *  series expansion
             *
             */
            if ( eps < 1.0e-3 )
            {
              sqrt_two =  1.414213562373095048801688724209698078569671875376948;
              theta = sqrt_two * sqrt(eps)*( 1 + ( (1.0/12.0) + ( (3.0/160.0) + ( (5.0/896.0) + ( (35.0/18432.0) + (63.0/90112.0) * eps) * eps) * eps) * eps) * eps);
            }
            else
              /*
               *  this is the regular no strings attached case
               */
            {
              theta = acos(c0_magnitude / c0_max);
            }

          /*printf("theta = %12.14lf\n", theta);*/

          /*
           * u from eqtn (23) hep-lat/0311018
           *
           * u = sqrt(c1 / 3) * cos(theta / 3)
           */
          u = sqrt(0.333333333333333333333333333333333
              * g_c1_smearing[x])*cos(0.333333333333333333333333333333333*theta);
          u_squared = u * u;

          /*
           * w from eqtn (24) hep-lat/0311018
           *
           * w = sqrt(c1) * sin(theta / 3)
           *
           */
          w = sqrt(g_c1_smearing[x])*sin(0.333333333333333333333333333333333*theta);
          w_squared = w * w;

          /*printf("u = %12.14lf w = %12.14lf\n", u, w);*/

          /*
           *  CORNER CASE 3:
           *  this case has been discussed in the paper
           *  and is spelled out in eqtn (33) hep-lat/0311018
           *
           *  1/6  = 0.1666666666...
           *  1/42 = 0.0238095238...
           */
          if(fabs(w) < 0.05)
          {
            xi_0 = 1.0 - 0.1666666666666666666666666666666666666666666666666667*w_squared*
              (1.0 - 0.05 * w_squared *
               (1.0 - 0.0238095238095238095238095238095238095238095238095238 * w_squared));
            xi_1 = -1.0 * ( 1.0/3.0) - (1.0/30.0)*w_squared*( 1.0 - (1.0/28.0 * w_squared*( 1.0 - (1.0/54.0)*w_squared)));

          }
          else
          {
            xi_0 = sin(w)/w;
            xi_1 = (cos(w) - xi_0)/w_squared;
          }

          /*printf("xi_0 = %12.14lf xi_1 = %12.14lf\n", xi_0, xi_1);*/

          /*
           * store a few values to save repeated calculations
           */
          cos_u = cos(u);
          sin_u = sin(u);
          cos_w = cos(w);
          sin_w = sin(w);
          cos_2u = cos(2.0 * u);
          sin_2u = sin(2.0 * u);
          u_cos_u = u * cos_u;
          u_sin_u = u * sin_u;
          u_cos_2u = u * cos_2u;
          u_sin_2u = u * sin_2u;
          denominator = 1.0/(9.0*u_squared-w_squared);

          /*
           * f_0, f_1 and f_2  from eqtn (29) hep-lat/0311018
           */
          tmp_d_0 = u_squared - w_squared;
          tmp_d_1 = 8.0 * u_squared * cos_w;
          tmp_d_2 = (3.0 * u_squared + w_squared) * xi_0;

          g_f0_smearing[x].re = (tmp_d_0 * cos_2u + cos_u*tmp_d_1 + 2*u_sin_u*tmp_d_2)*denominator;
          g_f0_smearing[x].im = (tmp_d_0 * sin_2u - sin_u*tmp_d_1 + 2*u_cos_u*tmp_d_2)*denominator;

          tmp_d_0 = (3.0 * u_squared - w_squared) * xi_0;

          g_f1_smearing[x].re = (2.0*(u_cos_2u - u_cos_u*cos_w) + tmp_d_0*sin_u) * denominator;
          g_f1_smearing[x].im = (2.0*(u_sin_2u + u_sin_u*cos_w) + tmp_d_0*cos_u) * denominator;

          tmp_d_0 = 3.0 * xi_0;

          g_f2_smearing[x].re = (cos_2u - cos_u * cos_w - u_sin_u * tmp_d_0) * denominator;
          g_f2_smearing[x].im = (sin_2u + sin_u * cos_w - u_cos_u * tmp_d_0) * denominator;

          /*printf("f0 = %12.14lf += i * %12.14lf\n", g_f0_smearing[x].re, g_f0_smearing[x].im);
            printf("f1 = %12.14lf + i * %12.14lf\n", g_f1_smearing[x].re, g_f1_smearing[x].im);
            printf("f2 = %12.14lf + i * %12.14lf\n", g_f2_smearing[x].re, g_f2_smearing[x].im);*/

          /*
           *  r10, r11, r12, r20, r21 and r22 
           *  from eqtns (60) to (65) in hep-lat/0311018
           */
          tmp_d_0 = u_squared - w_squared;
          tmp_d_1 = 8.0 * cos_w + (3.0*u_squared + w_squared)*xi_0;
          tmp_d_2 = 4.0 * u_squared * cos_w -(9.0 * u_squared + w_squared) * xi_0;

          g_r10_smearing[x].re = 2.0 * (u_cos_2u - sin_2u * tmp_d_0 + u_cos_u * tmp_d_1 - sin_u * tmp_d_2);
          g_r10_smearing[x].im = 2.0 * (u_sin_2u + cos_2u * tmp_d_0 - u_sin_u * tmp_d_1 - cos_u * tmp_d_2);

          tmp_d_0 = cos_w + 3.0 * xi_0;
          tmp_d_1 = 2.0 * cos_w + (w_squared - 3.0*u_squared)*xi_0;

          g_r11_smearing[x].re = 2.0 * ((cos_2u - 2.0 * u_sin_2u) + u_sin_u * tmp_d_0) - cos_u * tmp_d_1;
          g_r11_smearing[x].im = 2.0 * ((sin_2u + 2.0 * u_cos_2u) + u_cos_u * tmp_d_0) + sin_u * tmp_d_1;

          tmp_d_0 = cos_w - 3.0 * xi_0;

          g_r12_smearing[x].re = -2.0 * sin_2u - 3.0 * u_cos_u * xi_0 + sin_u * tmp_d_0;
          g_r12_smearing[x].im =  2.0 * cos_2u + 3.0 * u_sin_u * xi_0 + cos_u * tmp_d_0;

          tmp_d_0 = cos_w + xi_0 + 3.0 * u_squared * xi_1;

          g_r20_smearing[x].re = -2.0 * (cos_2u + u * (4.0 * u_cos_u * xi_0 - sin_u * tmp_d_0));
          g_r20_smearing[x].im = -2.0 * (sin_2u - u * (4.0 * u_sin_u * xi_0 + cos_u * tmp_d_0));

          tmp_d_0 = cos_w + xi_0 - 3.0 * u_squared * xi_1;

          g_r21_smearing[x].re =  2.0 * u_cos_u * xi_0 - sin_u * tmp_d_0;
          g_r21_smearing[x].im = -2.0 * u_sin_u * xi_0 - cos_u * tmp_d_0;

          tmp_d_0 = 3.0 * xi_1;

          g_r22_smearing[x].re =  cos_u * xi_0 - u_sin_u * tmp_d_0;
          g_r22_smearing[x].im = -sin_u * xi_0 - u_cos_u * tmp_d_0;

          /*printf("r10 = %12.14lf += i * %12.14lf\n", g_r10_smearing[x].re, g_r10_smearing[x].im);
            printf("r11 = %12.14lf += i * %12.14lf\n", g_r11_smearing[x].re, g_r11_smearing[x].im);
            printf("r12 = %12.14lf += i * %12.14lf\n", g_r12_smearing[x].re, g_r12_smearing[x].im);
            printf("r20 = %12.14lf += i * %12.14lf\n", g_r20_smearing[x].re, g_r20_smearing[x].im);
            printf("r21 = %12.14lf += i * %12.14lf\n", g_r21_smearing[x].re, g_r21_smearing[x].im);
            printf("r22 = %12.14lf += i * %12.14lf\n", g_r22_smearing[x].re, g_r22_smearing[x].im);*/

          denominator = 0.5 * denominator * denominator;
          g_b10_smearing[x].re = 

            tmp_d_0 = 2.0 * u;
          tmp_d_1 = 3.0 * u_squared - w_squared;
          tmp_d_2 = 2.0 * (15.0 * u_squared + w_squared);

          g_b10_smearing[x].re 
            =( tmp_d_0 * g_r10_smearing[x].re + tmp_d_1 * g_r20_smearing[x].re - tmp_d_2 * g_f0_smearing[x].re) * denominator;
          g_b10_smearing[x].im 
            =( tmp_d_0 * g_r10_smearing[x].im + tmp_d_1 * g_r20_smearing[x].im - tmp_d_2 * g_f0_smearing[x].im ) * denominator;

          g_b11_smearing[x].re 
            =( tmp_d_0 * g_r11_smearing[x].re + tmp_d_1 * g_r21_smearing[x].re - tmp_d_2 * g_f1_smearing[x].re ) * denominator;
          g_b11_smearing[x].im 
            =( tmp_d_0 * g_r11_smearing[x].im + tmp_d_1 * g_r21_smearing[x].im - tmp_d_2 * g_f1_smearing[x].im ) * denominator;

          g_b12_smearing[x].re 
            =( tmp_d_0 * g_r12_smearing[x].re + tmp_d_1 * g_r22_smearing[x].re - tmp_d_2 * g_f2_smearing[x].re ) * denominator;
          g_b12_smearing[x].im 
            =( tmp_d_0 * g_r12_smearing[x].im + tmp_d_1 * g_r22_smearing[x].im - tmp_d_2 * g_f2_smearing[x].im ) * denominator;

          tmp_d_0 = 3.0 * u;
          tmp_d_1 = 24.0 * u;

          g_b20_smearing[x].re=( g_r10_smearing[x].re 
              - tmp_d_0 * g_r20_smearing[x].re - tmp_d_1 * g_f0_smearing[x].re) * denominator;
          g_b20_smearing[x].im=( g_r10_smearing[x].im 
              - tmp_d_0 * g_r20_smearing[x].im - tmp_d_1 * g_f0_smearing[x].im ) * denominator;

          g_b21_smearing[x].re=( g_r11_smearing[x].re 
              - tmp_d_0 * g_r21_smearing[x].re - tmp_d_1 * g_f1_smearing[x].re ) * denominator;
          g_b21_smearing[x].im=( g_r11_smearing[x].im 
              - tmp_d_0 * g_r21_smearing[x].im - tmp_d_1 * g_f1_smearing[x].im ) * denominator;

          g_b22_smearing[x].re=( g_r12_smearing[x].re 
              - tmp_d_0 * g_r22_smearing[x].re - tmp_d_1 * g_f2_smearing[x].re ) * denominator;
          g_b22_smearing[x].im=( g_r12_smearing[x].im 
              - tmp_d_0 * g_r22_smearing[x].im - tmp_d_1 * g_f2_smearing[x].im ) * denominator;

          /*printf("b_denominator = %12.14lf\n", 1.0/denominator);
            printf("b10 = %12.14lf + i * %12.14lf\n", g_b10_smearing[x].re, g_b10_smearing[x].im);
            printf("b20 = %12.14lf + i * %12.14lf\n", g_b20_smearing[x].re, g_b20_smearing[x].im);
            printf("b11 = %12.14lf + i * %12.14lf\n", g_b11_smearing[x].re, g_b11_smearing[x].im);
            printf("b21 = %12.14lf + i * %12.14lf\n", g_b21_smearing[x].re, g_b21_smearing[x].im);
            printf("b12 = %12.14lf + i * %12.14lf\n", g_b12_smearing[x].re, g_b12_smearing[x].im);
            printf("b22 = %12.14lf + i * %12.14lf\n", g_b22_smearing[x].re, g_b22_smearing[x].im);*/

          if(c0_is_negative == 1)
          {
            /*printf("c0 is negative\n");*/
            g_b10_smearing[x].im *= -1.0;
            g_b11_smearing[x].re *= -1.0;
            g_b12_smearing[x].im *= -1.0;
            g_b20_smearing[x].re *= -1.0;
            g_b21_smearing[x].im *= -1.0;
            g_b22_smearing[x].re *= -1.0;

            g_f0_smearing[x].im *= -1.0;
            g_f1_smearing[x].re *= -1.0;
            g_f2_smearing[x].im *= -1.0;

          }
        }

        /*
         *  here we calculate B1, B2  eqtn (69) hep-lat/0311018
         *  
         *  B_i = bi0 + bi1 * Q + bi2 Q^2
         *
         *  we need to correct for the i factor we have in our matrices
         *  in order to get the same as chroma
         */
        _su3_one(tmp_su3_0);
        _complex_times_su3(g_B1_smearing[x][mu], g_b10_smearing[x], tmp_su3_0);

        tmp_c_0.re = 0.0;
        tmp_c_0.im = -1.0;
        _complex_times_su3(tmp_su3_0, g_b11_smearing[x], g_Q_smearing[x][mu]);
        _complex_times_su3(tmp_su3_1, tmp_c_0, tmp_su3_0);

        _su3_plus_su3(g_B1_smearing[x][mu], g_B1_smearing[x][mu], tmp_su3_1);

        tmp_c_0.re = -1.0;
        tmp_c_0.im = 0.0;
        _complex_times_su3(tmp_su3_0, g_b12_smearing[x], g_Q_squared_smearing[x][mu]);
        _complex_times_su3(tmp_su3_1, tmp_c_0, tmp_su3_0);
        _su3_plus_su3(g_B1_smearing[x][mu], g_B1_smearing[x][mu], tmp_su3_1);

        _su3_one(tmp_su3_0);
        _complex_times_su3(g_B2_smearing[x][mu], g_b20_smearing[x], tmp_su3_0);

        tmp_c_0.re = 0.0;
        tmp_c_0.im = -1.0;
        _complex_times_su3(tmp_su3_0, g_b21_smearing[x], g_Q_smearing[x][mu]);
        _complex_times_su3(tmp_su3_1, tmp_c_0, tmp_su3_0);

        _su3_plus_su3(g_B2_smearing[x][mu], g_B2_smearing[x][mu], tmp_su3_1);

        tmp_c_0.re = -1.0;
        tmp_c_0.im = 0.0;
        _complex_times_su3(tmp_su3_0, g_b22_smearing[x], g_Q_squared_smearing[x][mu]);
        _complex_times_su3(tmp_su3_1, tmp_c_0, tmp_su3_0);
        _su3_plus_su3(g_B2_smearing[x][mu], g_B2_smearing[x][mu], tmp_su3_1);

        /*if(x==0 && mu == 0)
          {
          printf("B1 =\n");
          print_su3_octave(&B1);
          printf("B2 =\n");
          print_su3_octave(&B2);
          }*/

        /*
         *  now we can calculate \Gamma eqtn (74) hep-lat/0311018
         */

        /*
         *  Tr(\Sigma' * B1 * U ) * Q
         */
        _su3_times_su3(tmp_su3_0, g_B1_smearing[x][mu], g_gauge_field_smeared[x][mu]);
        _trace_su3_times_su3(tmp_c_0, g_previous_stout_force_field[x][mu], tmp_su3_0);
        /* times -i in order to account for our conventions*/
        tmp_c_1.re = tmp_c_0.im;
        tmp_c_1.im = -tmp_c_0.re;
        _complex_times_su3(g_Gamma_smearing[x][mu], tmp_c_1, g_Q_smearing[x][mu]);

        /*if(x==0 )
          {
          printf("Gamma dummy_0 mu = %d\n", mu);
          print_su3(&(g_Gamma_smearing[0][0]));
          }*/

        /*
         *  Tr(\Sigma' * B2 * U ) * Q^2
         */
        _su3_times_su3(tmp_su3_0, g_B2_smearing[x][mu], g_gauge_field_smeared[x][mu]);
        _trace_su3_times_su3(tmp_c_0, g_previous_stout_force_field[x][mu], tmp_su3_0);
        /* times -1 in order to account for our conventions*/
        tmp_c_1.re = -tmp_c_0.re;
        tmp_c_1.im = -tmp_c_0.im;
        _complex_times_su3(tmp_su3_0, tmp_c_1, g_Q_squared_smearing[x][mu]);

        _su3_plus_su3(g_Gamma_smearing[x][mu], g_Gamma_smearing[x][mu], tmp_su3_0);

        /*if(x==0)
          {
          printf("Gamma dummy_1 mu = %d\n", mu);
          print_su3(&(tmp_su3_0));
          print_su3(&(g_Gamma_smearing[0][0]));
          }*/


        /*
         *  f_1 * U  *\Sigma'
         */
        _su3_times_su3(tmp_su3_0, g_gauge_field_smeared[x][mu], g_previous_stout_force_field[x][mu]);
        _complex_times_su3(tmp_su3_1, g_f1_smearing[x], tmp_su3_0);

        _su3_plus_su3(g_Gamma_smearing[x][mu], g_Gamma_smearing[x][mu], tmp_su3_1);

        /*if(x==0)
          {
          printf("Gamma dummy_2 mu = %d\n", mu);
          print_su3(&(tmp_su3_1));
          print_su3(&(g_Gamma_smearing[0][0]));
          }*/


        /*
         *  f_2 * Q * U  *\Sigma' + f_2 * U * \Sigma' * Q
         */
        _su3_times_su3(tmp_su3_0, g_gauge_field_smeared[x][mu], g_previous_stout_force_field[x][mu]);
        _su3_times_su3(tmp_su3_1, g_Q_smearing[x][mu], tmp_su3_0);
        _su3_times_su3(tmp_su3_2, tmp_su3_0,  g_Q_smearing[x][mu]);
        _su3_plus_su3(tmp_su3_0, tmp_su3_1, tmp_su3_2);
        _complex_times_su3(tmp_su3_1, g_f2_smearing[x], tmp_su3_0);
        /* times -i in order to account for our conventions*/
        tmp_c_1.re = 0.0;
        tmp_c_1.im = -1.0;
        _complex_times_su3(tmp_su3_0, tmp_c_1, tmp_su3_1);
        _su3_plus_su3(g_Gamma_smearing[x][mu], g_Gamma_smearing[x][mu], tmp_su3_0);

        /*if(x==0)
          {
          printf("Gamma dummy_3 mu = %d\n", mu);
          print_su3(&(tmp_su3_0));
          print_su3(&(g_Gamma_smearing[0][0]));
          }*/

        /*if(x==0 && mu == 0)
          {
          printf("Gamma Teil E\n");
          print_su3_octave(&(g_previous_stout_force_field[0][0]));
          printf("\n");
          print_su3_octave(&(g_stout_Gamma_field[0][0]));
          }*/

        /*
         *  now we can calculate \Lambda eqtn (73) hep-lat/0311018
         */

        /*
         *  first we store Gamma + Gamma^\dagger in tmp_su3_0
         */
        tmp_su3_0.c00.re = g_Gamma_smearing[x][mu].c00.re + g_Gamma_smearing[x][mu].c00.re;
        tmp_su3_0.c00.im = 0.0;
        tmp_su3_0.c11.re = g_Gamma_smearing[x][mu].c11.re + g_Gamma_smearing[x][mu].c11.re;
        tmp_su3_0.c11.im = 0.0;
        tmp_su3_0.c22.re = g_Gamma_smearing[x][mu].c22.re + g_Gamma_smearing[x][mu].c22.re;
        tmp_su3_0.c22.im = 0.0;

        tmp_su3_0.c01.re = g_Gamma_smearing[x][mu].c01.re + g_Gamma_smearing[x][mu].c10.re;
        tmp_su3_0.c01.im = g_Gamma_smearing[x][mu].c01.im - g_Gamma_smearing[x][mu].c10.im;
        tmp_su3_0.c10.re = tmp_su3_0.c01.re;
        tmp_su3_0.c10.im = -tmp_su3_0.c01.im;

        tmp_su3_0.c02.re = g_Gamma_smearing[x][mu].c02.re + g_Gamma_smearing[x][mu].c20.re;
        tmp_su3_0.c02.im = g_Gamma_smearing[x][mu].c02.im - g_Gamma_smearing[x][mu].c20.im;
        tmp_su3_0.c20.re = tmp_su3_0.c02.re;
        tmp_su3_0.c20.im = -tmp_su3_0.c02.im;

        tmp_su3_0.c21.re = g_Gamma_smearing[x][mu].c21.re + g_Gamma_smearing[x][mu].c21.re;
        tmp_su3_0.c21.im = g_Gamma_smearing[x][mu].c21.im - g_Gamma_smearing[x][mu].c21.im;
        tmp_su3_0.c12.re = tmp_su3_0.c21.re;
        tmp_su3_0.c12.im = -tmp_su3_0.c21.im;

        /* tmp_c_0 holds one third of the trace*/
        tmp_c_0.re = 0.333333333333333333333333 * (tmp_su3_0.c00.re + tmp_su3_0.c11.re + tmp_su3_0.c22.re); 
        tmp_c_0.im = 0.333333333333333333333333 * (tmp_su3_0.c00.im + tmp_su3_0.c11.im + tmp_su3_0.c22.im); 
        tmp_su3_0.c00.re -= tmp_c_0.re;
        tmp_su3_0.c00.im -= tmp_c_0.im;
        tmp_su3_0.c11.re -= tmp_c_0.re;
        tmp_su3_0.c11.im -= tmp_c_0.im;
        tmp_su3_0.c22.re -= tmp_c_0.re;
        tmp_su3_0.c22.im -= tmp_c_0.im;

        tmp_c_0.re = 0.5;
        tmp_c_0.im = 0.0;

        _complex_times_su3(g_Lambda_smearing[x][mu], tmp_c_0, tmp_su3_0);

        /*if(x==0 && mu == 0)
          {
          printf("Lambda Teil A\n");
          print_su3_octave(&(g_stout_Lambda_field[0][0]));
          printf("\n");
          print_su3_octave(&(g_stout_Gamma_field[0][0]));
          }*/
      }
      /*printf("DUMP OF c0 at mu = %d\n", mu);
        print_scalar_real_field_to_screen(g_c0_smearing);
        printf("DUMP OF c1 at mu = %d\n", mu);
        print_scalar_real_field_to_screen(g_c1_smearing);
        printf("DUMP OF f0 at mu = %d\n", mu);
        print_scalar_complex_field_to_screen(g_f0_smearing);
        printf("DUMP OF f1 at mu = %d\n", mu);
        print_scalar_complex_field_to_screen(g_f1_smearing);
        printf("DUMP OF f2 at mu = %d\n", mu);
        print_scalar_complex_field_to_screen(g_f2_smearing);*/
    }

    /*printf("DUMP OF g_C_smearing\n");
      print_config_to_screen(g_C_smearing);
      printf("DUMP OF g_Q_smearing\n");
      print_config_to_screen(g_Q_smearing);
      printf("DUMP OF g_Q_squared_smearing\n");
      print_config_to_screen(g_Q_squared_smearing);
      printf("DUMP OF g_B1_smearing\n");
      print_config_to_screen(g_B1_smearing);
      printf("DUMP OF g_B2_smearing\n");
      print_config_to_screen(g_B2_smearing);
      printf("DUMP OF g_Gamma_smearing\n");
      print_config_to_screen(g_Gamma_smearing);
      printf("DUMP OF g_Lambda_smearing\n");
      print_config_to_screen(g_Lambda_smearing);*/

    /*
     *  now we can calculate the forcefield \Signma_\mu(x)
     *  eqtn (75) in hep-lat/0311018
     */
    for(x = 0; x < VOLUME; x++)
      for(mu = 0; mu < dim; mu++)
      {


        /*
         *  first we do f0 * 1 + f1 * Q + f2 * q^2
         */
        _su3_zero(tmp_su3_0);

        tmp_su3_0.c00.re += g_f0_smearing[x].re;
        tmp_su3_0.c00.im += g_f0_smearing[x].im;
        tmp_su3_0.c11.re += g_f0_smearing[x].re;
        tmp_su3_0.c11.im += g_f0_smearing[x].im;
        tmp_su3_0.c22.re += g_f0_smearing[x].re;
        tmp_su3_0.c22.im += g_f0_smearing[x].im;

        _complex_times_su3(tmp_su3_1, g_f1_smearing[x], g_Q_smearing[x][mu]);
        /* here we compensate for our i factors*/
        tmp_c_0.re =  0.0;
        tmp_c_0.im = -1.0;
        _complex_times_su3(tmp_su3_2, tmp_c_0, tmp_su3_1);
        _su3_plus_su3(tmp_su3_0, tmp_su3_0, tmp_su3_2);

        _complex_times_su3(tmp_su3_1, g_f2_smearing[x], g_Q_squared_smearing[x][mu]);
        /* here we compensate for our i factors*/
        tmp_c_0.re = -1.0;
        tmp_c_0.im =  0.0;
        _complex_times_su3(tmp_su3_2, tmp_c_0, tmp_su3_1);
        _su3_plus_su3(tmp_su3_0, tmp_su3_0, tmp_su3_2);

        /*if(x==0 )
          {
          printf("Spiderschwein mu = %d\n", mu);
          print_su3(&(tmp_su3_0));
          printf("\n");
          }*/

        /*
         *  \Sigma' * (above)
         */
        _su3_times_su3(g_stout_force_field[x][mu], g_previous_stout_force_field[x][mu], tmp_su3_0);

        /*if(x==0 )
          {
          printf("Sigma A mu = %d\n", mu);
          print_su3(&(g_stout_force_field[x][mu]));
          printf("\n");
          }*/


        /*
         *  i * (C_\mu)^\dagger * Lambda
         */
        _su3d_times_su3(tmp_su3_0, g_C_smearing[x][mu], g_Lambda_smearing[x][mu]);
        _itimes_su3(tmp_su3_1, tmp_su3_0);
        _su3_plus_su3(g_stout_force_field[x][mu], g_stout_force_field[x][mu], tmp_su3_1);

        /*if(x==0 )
          {*/
        /*printf("Spiderschwein x= %d mu = %d\n", x, mu);
          print_su3(&(g_stout_force_field[x][mu]));
          printf("\n");*/
        /*}*/

        /*
         *  no we do the sum eqtn (75) hep-lat/0311018
         *  since our rho is mu/nu independent we pull in front of the sum
         */
        _su3_zero(tmp_su3_3);
        for(nu=0; nu < dim; nu++)
          if(nu != mu)
          {
            _su3d_times_su3(tmp_su3_0, g_gauge_field_smeared[x][nu], g_Lambda_smearing[x][nu]);
            _su3d_times_su3(tmp_su3_1, g_gauge_field_smeared[g_iup[x][nu]][mu], tmp_su3_0);
            _su3_times_su3(tmp_su3_2, g_gauge_field_smeared[g_iup[x][mu]][nu], tmp_su3_1);


            /*if(x==0 && mu == 0 && nu ==1)
              {
              printf("Sigma D mu = %d\n", mu);*/
            /*print_su3(&(g_stout_force_field[0][0]));*/
            /*print_su3(&(tmp_su3_2));
              printf("\n");
              }*/

            _su3_times_su3(tmp_su3_0, g_Lambda_smearing[g_idn[x][nu]][mu], g_gauge_field_smeared[g_idn[x][nu]][nu]);
            _su3d_times_su3(tmp_su3_1, g_gauge_field_smeared[g_idn[x][nu]][mu], tmp_su3_0);
            _su3d_times_su3(tmp_su3_0, g_gauge_field_smeared[g_iup[g_idn[x][nu]][mu]][nu], tmp_su3_1);
            _su3_plus_su3(tmp_su3_2, tmp_su3_2, tmp_su3_0);


            _su3d_times_su3(tmp_su3_0, g_gauge_field_smeared[g_idn[x][nu]][mu], g_gauge_field_smeared[g_idn[x][nu]][nu]);
            _su3_times_su3(tmp_su3_1, g_Lambda_smearing[g_iup[g_idn[x][nu]][mu]][nu], tmp_su3_0);
            _su3d_times_su3(tmp_su3_0, g_gauge_field_smeared[g_iup[g_idn[x][nu]][mu]][nu], tmp_su3_1);
            _su3_plus_su3(tmp_su3_2, tmp_su3_2, tmp_su3_0);

            _su3_times_su3(tmp_su3_0, g_Lambda_smearing[g_idn[x][nu]][nu], g_gauge_field_smeared[g_idn[x][nu]][nu]);
            _su3d_times_su3(tmp_su3_1, g_gauge_field_smeared[g_idn[x][nu]][mu], tmp_su3_0);
            _su3d_times_su3(tmp_su3_0, g_gauge_field_smeared[g_iup[g_idn[x][nu]][mu]][nu], tmp_su3_1);
            _su3_minus_su3(tmp_su3_2, tmp_su3_2, tmp_su3_0);


            _su3_dagger(tmp_su3_0, g_gauge_field_smeared[x][nu]);
            _su3d_times_su3(tmp_su3_1, g_gauge_field_smeared[g_iup[x][nu]][mu], tmp_su3_0);
            _su3_times_su3(tmp_su3_0, g_gauge_field_smeared[g_iup[x][mu]][nu], tmp_su3_1);
            _su3_times_su3(tmp_su3_1, g_Lambda_smearing[g_iup[x][mu]][nu], tmp_su3_0);
            _su3_minus_su3(tmp_su3_2, tmp_su3_2, tmp_su3_1);


            _su3_dagger(tmp_su3_0, g_gauge_field_smeared[x][nu]);
            _su3_times_su3(tmp_su3_1, g_Lambda_smearing[g_iup[x][nu]][mu], tmp_su3_0);
            _su3d_times_su3(tmp_su3_0, g_gauge_field_smeared[g_iup[x][nu]][mu], tmp_su3_1);
            _su3_times_su3(tmp_su3_1, g_gauge_field_smeared[g_iup[x][mu]][nu], tmp_su3_0);
            _su3_plus_su3(tmp_su3_2, tmp_su3_2, tmp_su3_1);

            _su3_plus_su3(tmp_su3_3, tmp_su3_3, tmp_su3_2);
          }

        /*if(x==0 && mu == 0)
          {
          printf("Sigma Teil C\n");
          print_su3_octave(&(g_stout_force_field[0][0]));
          printf("\n");
          print_su3_octave(&(g_stout_Gamma_field[0][0]));
          printf("\n");
          print_su3_octave(&(g_previous_stout_force_field[0][0]));
          printf("\n");
          }*/

        /*
         *  now we multiply with i * rho and subtract from the first part
         */
        tmp_c_0.re = 0.0;
        tmp_c_0.im = stout_rho;

        /*if(x==0)
          {
          printf("tmp_su3_3\n");
          print_su3(&(tmp_su3_3));
          printf("tmp_c_0.re = %lf   tmp_c_0.im = %lf\n", );
          }*/
        _complex_times_su3(tmp_su3_0, tmp_c_0, tmp_su3_3);
        /*if(x==0)
          {
          printf("tmp_su3_0\n");
          print_su3(&(tmp_su3_0));
          }*/
        /* here we compensate for our additional i factors*/
        tmp_c_0.re = 0.0;
        tmp_c_0.im = 1.0;
        _complex_times_su3(tmp_su3_1, tmp_c_0, tmp_su3_0);
        _su3_minus_su3(g_stout_force_field[x][mu], g_stout_force_field[x][mu], tmp_su3_1);

        /*if(x==0)
          {
          printf("Sigma E mu = %d\n", mu);*/
        /*print_su3(&(g_stout_force_field[0][0]));*/
        /*print_su3(&(g_stout_force_field[x][mu]));
          printf("tmp_su3_1\n");
          print_su3(&(tmp_su3_1));
          printf("tmp_su3_3\n");
          print_su3(&(tmp_su3_3));
          printf("\n");
          }*/
      }

    /*printf("\nENDENDENDENDENDENDENDEND\n");*/
    /*exit(7);*/
    /*for(x = 0; x < VOLUME; x++)
      for(mu = 0; mu < 4; mu++)
      _su3_zero(g_previous_stout_force_field[x][mu]);*/
    /*
     *  finally we save the \Sigma field, because we need it as \Sigma'
     *  in the next iteration
     */
    /*for(x = 0; x < VOLUME; x++)
      for(mu = 0; mu < dim; mu++)
      {

      _su3_assign(g_previous_stout_force_field[x][mu], g_stout_force_field[x][mu]);
      }*/


  }
  /*printf("CHECCJHECKCHECK\n");
    tmp_su3_0.c00.re=1.0;
    tmp_su3_0.c11.re=1.0;
    tmp_su3_0.c22.re=-2.0;
    tmp_su3_0.c00.im=0.0;
    tmp_su3_0.c11.im=0.0;
    tmp_su3_1.c22.im=-0.0;

    tmp_su3_0.c01.re=-2.0;
    tmp_su3_0.c01.im=-10.0;
    tmp_su3_0.c10.re=tmp_su3_0.c01.re;
    tmp_su3_0.c10.im=-tmp_su3_0.c01.im;

    tmp_su3_0.c02.re=7.0;
    tmp_su3_0.c02.im=5.0;
    tmp_su3_0.c20.re=tmp_su3_0.c02.re;
    tmp_su3_0.c20.im=-tmp_su3_0.c02.im;

    tmp_su3_0.c12.re=-4.0;
    tmp_su3_0.c12.im=-10.0;
    tmp_su3_0.c21.re=tmp_su3_0.c12.re;
    tmp_su3_0.c21.im=-tmp_su3_0.c12.im;*/

  /*tmp_c_0.re = 0.0;
    tmp_c_0.im = -1.0;
    printf("1\n");
    print_su3_octave(&(tmp_su3_0));
    printf("2\n");
    _complex_times_su3(tmp_su3_1, tmp_c_0, tmp_su3_0);
    print_su3_octave(&(tmp_su3_1));
    _trace_lambda(tttt, tmp_su3_1);
    printf("%f %f %f %f %f %f %f %f\n", tttt.d1, tttt.d2, tttt.d3, tttt.d4, tttt.d5, tttt.d6, tttt.d7, tttt.d8);
    tttt.d1 /= -2.0;
    tttt.d2 /= -2.0;
    tttt.d3 /= -2.0;
    tttt.d4 /= -2.0;
    tttt.d5 /= -2.0;
    tttt.d6 /= -2.0;
    tttt.d7 /= -2.0;
    tttt.d8 /= -2.0;
    printf("4\n");
    printf("%f %f %f %f %f %f %f %f\n", tttt.d1, tttt.d2, tttt.d3, tttt.d4, tttt.d5, tttt.d6, tttt.d7, tttt.d8);
    _make_su3(tmp_su3_0, tttt);
    printf("5\n");
    print_su3_octave(&(tmp_su3_0));*/

  /*printf("df0 before = %f %f %f %f %f %f %f %f\n", before_force.d1, before_force.d2, before_force.d3, before_force.d4, before_force.d5, before_force.d6, before_force.d7, before_force.d8);*/
  /*_trace_lambda(tttt, g_stout_force_field[0][0]);
    printf("%f %f %f %f %f %f %f %f\n", tttt.d1, tttt.d2, tttt.d3, tttt.d4, tttt.d5, tttt.d6, tttt.d7, tttt.d8);
    tttt.d1 /= -2.0;
    tttt.d2 /= -2.0;
    tttt.d3 /= -2.0;
    tttt.d4 /= -2.0;
    tttt.d5 /= -2.0;
    tttt.d6 /= -2.0;
    tttt.d7 /= -2.0;
    tttt.d8 /= -2.0;
    printf("%f %f %f %f %f %f %f %f\n", tttt.d1, tttt.d2, tttt.d3, tttt.d4, tttt.d5, tttt.d6, tttt.d7, tttt.d8);
    _make_su3(tmp_su3_0, before_force);
    print_su3_octave(&(tmp_su3_0));
    print_su3_octave(&(g_stout_force_field[0][0]));*/
}
