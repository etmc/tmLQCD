#include <complex.h>
#include <math.h>

#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "stout_smear.h"
#include "read_input.h"


/*void stout_smear_force(int x, int mu, su3 force){*/
void stout_smear_force()
{

  extern su3 * gauge_field_smear_iterations;
  extern su3 *** g_gauge_field_smear_iterations;
  extern su3 * stout_force_field;
  extern su3 ** g_stout_force_field;
  extern su3 * previous_stout_force_field;
  extern su3 ** g_previous_stout_force_field;
  extern su3 * Q_smear_iterations;
  extern su3 *** g_Q_smear_iterations;
  extern su3 * C_smear_iterations;
  extern su3 *** g_C_smear_iterations;
  extern su3 * stout_Lambda_field;
  extern su3 ** g_stout_Lambda_field;
  extern su3 * stout_Gamma_field;
  extern su3 ** g_stout_Gamma_field;

  const int dim = 4;

  int x, mu, nu, stout_iter_counter;
  double magnitude, theta, u, u_squared, w, w_squared, xi_0, xi_1;
  double tmp_d_0, tmp_d_1, tmp_d_2;
  complex c0, c1, c0_max, c1_max, tmp_c_0;
  complex h_0, h_1, h_2, f_0, f_1, f_2; 
  complex r_0_1, r_0_2, r_1_1, r_1_2, r_2_1, r_2_2; 
  complex b_1_0, b_2_0, b_1_1, b_2_1, b_1_2, b_2_2; 
  su3 Q_squared, tmp_su3_0, tmp_su3_1, tmp_su3_2; 
  su3 B1, B2; 

  printf("Running stout_smear_force()\n");

  /*
   *  first we save the original force field, so we have \Sigma'
   *  available when iterating thru  eqtn(75) hep-lat/0311018
   */
  for(x = 0; x < VOLUME*dim; x++)
    /*for(mu = 0; mu < dim; mu++)*/
    _su3_assign(previous_stout_force_field[x], stout_force_field[x]);

  for(stout_iter_counter = stout_no_iter-2; stout_iter_counter > -1; stout_iter_counter--)
  {
    for(x = 0; x < VOLUME; x++)
      for(mu = 0; mu < dim; mu++)
      {

        /*
         * first we calculate u and w
         */
        /*
         * c0 from eqtn (14) hep-lat/0311018
         *
         * c0 = (Tr Q^3)/3
         */
        _su3_times_su3(Q_squared, g_Q_smear_iterations[stout_iter_counter][x][mu], g_Q_smear_iterations[stout_iter_counter][x][mu]);

        _trace_su3_times_su3(c0, Q_squared, g_Q_smear_iterations[stout_iter_counter][x][mu]);
        c0.re *= -0.33333333333333333333333333333333333333333;
        c0.im *= -0.33333333333333333333333333333333333333333;

        /*
         * c1 from eqtn (15) hep-lat/0311018
         *
         * c1 = (Tr Q^2)/2
         */
        /*printf("ABABA\n");
          print_su3_octave(&Q_squared);*/
        c1.re = -0.5 * (Q_squared.c00.re + Q_squared.c11.re + Q_squared.c22.re);
        c1.im = -0.5 * (Q_squared.c00.im + Q_squared.c11.im + Q_squared.c22.im);

        
        /*
         * c0_max from eqtn (17) hep-lat/0311018
         *
         * c0_max = 2*(c1/3)^1.5
         *
         * sqrt(2)/(3)^1.5 = 0.27216552697590867757747600830065459910732749785074112538
         */


        /*
         *  as stated in the chroma sources ("stout_utils.cc")
         *  we need to watch out for the case c0 \approx 0
         *  this case is ommitted in hep-lat/0311018
         */
        if (c1.re < 1e-4)
        /*if (c1.re*c1.re + c1.im * c1.im < 1e-4)*/
        {  

          f_0.re = 1-c0.re * c0.re / 720.0;
          f_0.im =  -(c0.re / 6.0) * (1.0 - (c1.re / 20)
              * ( 1.0 - (c1.re / 42.0)));

          f_1.re =  c0.re / 24.0 * (1.0 - c1.re / 15.0 * (1.0 - 3.0 * c1.re / 112.0)) ;
          f_1.im =  1.0 - c1.re / 6.0 * (1.0 - c1.re / 20.0 * (1.0 - c1.re / 42.0))
            - c0.re * c0.re / 5040.0;

          f_2.re = 0.5 * (-1.0 + c1.re / 12.0 * (1.0 - c1.re / 30.0 
                * (1.0 - c1.re / 56.0)) + c0.re * c0.re / 20160.0);
          f_2.im = 0.5 * (c0.re / 60.0 * (1.0 - c1.re / 21.0 * (1.0 - c1.re / 48.0)));

          b_2_0.re = -c0.re / 360.0;
          b_2_0.im =  -(1.0 / 6.0) * (1.0 - (c1.re / 20.0) * (1.0 - c1.re / 42.0));
          
          b_1_0.re = 0.0;
          b_1_0.im = (c0.re / 120.0)*(1.0 - c1.re / 21.0);

          b_2_1.re = (1.0 / 24.0) * (1.0 - c1.re / 15.0 * (1.0 - 3.0 * c1.re / 112.0));
          b_2_1.im = -c0.re / 2520.0;

          b_1_1.re = -c0.re / 360.0 * (1.0 - 3.0 * c1.re / 56.0);
          b_1_1.im = -1.0 / 6.0 * (1.0 - c1.re / 10.0 * (1.0 - c1.re / 28.0));

          b_2_2.re = 0.5 * c0.re / 10080.0;
          b_2_2.im = 0.5 * (1.0 / 60.0 * (1.0 - c1.re / 21.0 * (1.0 - c1.re / 48.0)) );

          b_1_2.re = 0.5 * (1.0 / 12.0 * (1.0 - (2.0 * c1.re / 30.0) * (1.0 - 3.0*c1.re/112.0)) );
          b_1_2.im = 0.5*( -c0.re/1260.0 * (1.0-c1.re/24.0) );

        }
        else
        {
          /*  tmp_c_1 = c1^3 */
          tmp_c_0.re = (-3.0 * c1.im * c1.im + c1.re * c1.re) * c1.re;
          tmp_c_0.im = ( 3.0 * c1.re * c1.re - c1.im * c1.im) * c1.im;

          magnitude = sqrt(tmp_c_0.re * tmp_c_0.re + tmp_c_0.im * tmp_c_0.im);

          /*this is just the square root and not needed*/
          /*c0_max.re = 0.7071067811865475244008443621048490392848359*sqrt(magnitude + tmp.re); 
            if(tmp.im > 0.0)
            c0_max.im = 0.7071067811865475244008443621048490392848359*sqrt(magnitude - tmp.re);
            else
            c0_max.im = -0.7071067811865475244008443621048490392848359*sqrt(magnitude - tmp.re);*/

          /*c0_max.re = 0.27216552697590867757747600830065459910732749785074112538*sqrt(magnitude + tmp.re); 
            if(tmp.im > 0.0)
            c0_max.im = 0.27216552697590867757747600830065459910732749785074112538*sqrt(magnitude - tmp.re);
            else
            c0_max.im = -0.27216552697590867757747600830065459910732749785074112538*sqrt(magnitude - tmp.re);*/

          /*printf("OPOII %f\n", c1.re * c1.re * c1.re);*/
          c0_max.re = 0.3849001794597505841196932900549*sqrt(c1.re * c1.re * c1.re);
          /*c0_max.im = 0.0;*/

          /*
           * theta from eqtn (25) hep-lat/0311018
           *
           * theta = acos()
           */
          theta = acos(c0.re / c0_max.re);

          /*
           * u from eqtn (23) hep-lat/0311018
           *
           * u = sqrt(c1 / 3) * cos(theta / 3)
           */
          u=sqrt(0.333333333333333333333333333333333*c1.re)*cos(0.333333333333333333333333333333333*theta);
          u_squared = u * u;

          /*
           * w from eqtn (24) hep-lat/0311018
           *
           * w = sqrt(c1) * sin(theta / 3)
           *
           * 1/6  = 0.1666666666...
           * 1/42 = 0.0238095238...
           */
          w=sqrt(c1.re)*sin(0.333333333333333333333333333333333*theta);
          w_squared = w * w;

          /*
           * xi_0 from eqtn (33) hep-lat/0311018
           */
          if(fabs(w) < 0.05)
            xi_0=1
              -0.1666666666666666666666666666666666666666666666666667*w_squared*
              (1.0 - 0.05 * w_squared * 
               (1 - 0.0238095238095238095238095238095238095238095238095238 * w_squared));
          else
            xi_0 = sin(w)/w;

          /*
           * xi_1 from eqtn (67) hep-lat/0311018
           */
          xi_1 = (cos(w) - xi_0)/w_squared;

          /*
           * h_0, h_1 and h_2  from eqtn (30), (31) and (32) hep-lat/0311018
           */
          h_0.re = (u_squared - w_squared) * cos(2*u) + cos(-u)*8.0*u_squared*cos(w)
            - sin(-u) * 2.0 * u * (3*u_squared + w_squared)*xi_0;
          h_0.im = (u_squared - w_squared) * sin(2*u) + sin(-u)*8.0*u_squared*cos(w)
            + cos(-u) * 2.0 * u * (3*u_squared + w_squared)*xi_0;

          h_1.re = 2*u*cos(2.0*u) - cos(-u)*2*u*cos(w) - sin(-u)*(3*u_squared-w_squared)*xi_0;
          h_1.im = 2*u*sin(2.0*u) - sin(-u)*2*u*cos(w) + cos(-u)*(3*u_squared-w_squared)*xi_0;

          h_2.re = cos(2*u) - cos(-u) * cos (w) + sin(-u) * 3.0 * u * xi_0;
          h_2.im = sin(2*u) - sin(-u) * cos (w) + cos(-u) * 3.0 * u * xi_0;


          tmp_d_0 = 1.0/(9.0*u_squared-w_squared); 

          f_0.re = tmp_d_0 * h_0.re; 
          f_0.im = tmp_d_0 * h_0.im; 
          f_1.re = tmp_d_0 * h_1.re; 
          f_1.im = tmp_d_0 * h_1.im; 
          f_2.re = tmp_d_0 * h_2.re; 
          f_2.im = tmp_d_0 * h_2.im; 

          r_0_1.re = 2.0*(u*cos(2.0*u) - sin(2.0 * u)*(u_squared -w_squared) 
              + 8.0 * cos(-u) * u * cos(w) + u * cos(-u) * (3*u_squared+w_squared)*xi_0 + 4.0*u_squared*cos(w)*sin(-u)) - sin(-u)*(9.0*u_squared+w_squared)*xi_0;
          r_0_1.im = 2.0*(u*sin(2.0*u) + cos(2.0 * u)*(u_squared -w_squared) 
              + 8.0 * sin(-u) * u * cos(w) + u * sin(-u) * (3*u_squared+w_squared)*xi_0 + 4.0*u_squared*cos(w)*cos(-u) + cos(-u)*(9.0*u_squared+w_squared)*xi_0);

          r_0_2.re = -2.0 * cos(2.0*u) - 2.0 * u * sin(-u) * cos(w) - 2.0 * u * sin(-u)* xi_0 - 6.0 * sin(-u) * u_squared * xi_1;
          r_0_2.im = -2.0 * u * sin(-u) + 2.0 * u * cos(-u) * cos(w) + 2.0 * u * cos(-u) * xi_0 + 6.0 * u * u_squared * cos(-u) * xi_1 - 8.0 * sin(-u) * u * xi_0; 


          r_1_1.re = 2.0 * cos(2*u) - 4.0 * u * sin(2*u) - 2.0 * cos(-u) * cos(w) - cos(-u) * (w_squared - 3.0 * u_squared) * xi_0 - 2.0 * u * sin(-u) * cos(w) - 6.0 * u * sin(-u) * xi_0;
          r_1_1.im = 2.0 * sin(2*u) + 4.0 * u * cos(2*u) + 2.0 * u * cos(-u) * cos(w) + 6.0 * u * cos(-u) * xi_0 - 2.0 * sin(-u) * cos(w) - sin(-u) * (w_squared - 3.0 * u_squared) * xi_0;

          r_1_2.re = cos(-u) *2.0 * u *xi_0 + sin(-u) * cos(w) + sin(-u) * xi_0 - 3.0 * u_squared * sin(-u) * xi_1;
          r_1_2.im = -cos(-u) * cos(w) - cos(-u) * xi_0 + 3.0 * u_squared * cos(-u) * xi_1 + 2.0 * u * sin(-u) * xi_0;

          r_2_1.re = -2.0 * sin(2.0 * u) - 3.0 * cos(-u) * u * xi_0 - sin(-u) * cos(w) - 3.0 * sin(-u) * xi_0;
          r_2_1.im = 2.0 * cos(-u) * cos(w) + 3.0 * cos(-u) * xi_0 - 3.0 * u * sin(-u) * xi_0;

          r_2_2.re = cos(-u) * xi_0 + 3.0 * u * sin(-u) * xi_1;
          r_2_2.im = sin(-u) * xi_0 - 3.0 * cos(-u) * u * xi_1;


          tmp_d_0 = 9.0 * u_squared - w_squared;
          tmp_d_0 = 1.0 / (2.0 * tmp_d_0 * tmp_d_0);

          tmp_d_1 = 3.0 * u_squared - w_squared;
          tmp_d_2 = 15.0 * u_squared + w_squared;


          b_1_0.re = tmp_d_0 * (2.0 * u * r_0_1.re + tmp_d_1 * r_0_2.re - 2.0 * tmp_d_2 * f_0.re);
          b_1_0.im = tmp_d_0 * (2.0 * u * r_0_1.im + tmp_d_1 * r_0_2.im - 2.0 * tmp_d_2 * f_0.im);

          b_2_0.re = tmp_d_0 * (r_0_1.re - 3.0 * u * r_0_2.re - 24.0 * u *tmp_d_2 * f_0.re);
          b_2_0.im = tmp_d_0 * (r_0_1.im - 3.0 * u * r_0_2.im - 24.0 * u *tmp_d_2 * f_0.im);

          b_1_1.re = tmp_d_0 * (2.0 * u * r_1_1.re + tmp_d_1 * r_1_2.re - 2.0 * tmp_d_2 * f_1.re);
          b_1_1.im = tmp_d_0 * (2.0 * u * r_1_1.im + tmp_d_1 * r_1_2.im - 2.0 * tmp_d_2 * f_1.im);

          b_2_1.re = tmp_d_0 * (r_1_1.re - 3.0 * u * r_1_2.re - 24.0 * u *tmp_d_2 * f_1.re);
          b_2_1.im = tmp_d_0 * (r_1_1.im - 3.0 * u * r_1_2.im - 24.0 * u *tmp_d_2 * f_1.im);

          b_1_2.re = tmp_d_0 * (2.0 * u * r_2_1.re + tmp_d_1 * r_2_2.re - 2.0 * tmp_d_2 * f_2.re);
          b_1_2.im = tmp_d_0 * (2.0 * u * r_2_1.im + tmp_d_1 * r_2_2.im - 2.0 * tmp_d_2 * f_2.im);

          b_2_2.re = tmp_d_0 * (r_2_1.re - 3.0 * u * r_2_2.re - 24.0 * u *tmp_d_2 * f_2.re);
          b_2_2.im = tmp_d_0 * (r_2_1.im - 3.0 * u * r_2_2.im - 24.0 * u *tmp_d_2 * f_2.im);
        }

        /*
         * here we calculate B1, B2  eqtn (69) hep-lat/0311018
         */
        _complex_times_su3(tmp_su3_0, b_1_2, g_Q_smear_iterations[stout_iter_counter][x][mu]);
        _su3_one(tmp_su3_1);
        _complex_times_su3(tmp_su3_2, b_1_1, tmp_su3_1);
        _su3_plus_su3(tmp_su3_1 ,tmp_su3_0, tmp_su3_2);
        _su3_times_su3(tmp_su3_0, tmp_su3_1, g_Q_smear_iterations[stout_iter_counter][x][mu]);
        _su3_one(tmp_su3_1);
        _complex_times_su3(tmp_su3_2, b_1_0, tmp_su3_1);
        _su3_plus_su3(B1 ,tmp_su3_0, tmp_su3_2);

        _complex_times_su3(tmp_su3_0, b_2_2, g_Q_smear_iterations[stout_iter_counter][x][mu]);
        _su3_one(tmp_su3_1);
        _complex_times_su3(tmp_su3_2, b_2_1, tmp_su3_1);
        _su3_plus_su3(tmp_su3_1 ,tmp_su3_0, tmp_su3_2);
        _su3_times_su3(tmp_su3_0, tmp_su3_1, g_Q_smear_iterations[stout_iter_counter][x][mu]);
        _su3_one(tmp_su3_1);
        _complex_times_su3(tmp_su3_2, b_2_0, tmp_su3_1);
        _su3_plus_su3(B2, tmp_su3_0, tmp_su3_2);

        /*
         *  now we can calculate \Gamma eqtn (74) hep-lat/0311018
         */

        /*
         *  Tr(\Sigma' * B1 * U ) * Q
         */
        _su3_times_su3(tmp_su3_0, B1, g_gauge_field_smear_iterations[stout_iter_counter][x][mu]);
        _trace_su3_times_su3(tmp_c_0, g_previous_stout_force_field[x][mu], tmp_su3_0);
        _complex_times_su3(g_stout_Gamma_field[x][mu], tmp_c_0, g_Q_smear_iterations[stout_iter_counter][x][mu]);

        /*
         *  Tr(\Sigma' * B2 * U ) * Q^2
         */
        _su3_times_su3(tmp_su3_0, B2, g_gauge_field_smear_iterations[stout_iter_counter][x][mu]);
        _trace_su3_times_su3(tmp_c_0, g_previous_stout_force_field[x][mu], tmp_su3_0);
        _complex_times_su3(tmp_su3_0, tmp_c_0, Q_squared);

        _su3_plus_su3(g_stout_Gamma_field[x][mu], g_stout_Gamma_field[x][mu], tmp_su3_0);

        /*
         *  f_1 * U  *\Sigma' 
         */
        _su3_times_su3(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][x][mu], g_previous_stout_force_field[x][mu]);
        _complex_times_su3(tmp_su3_0, f_1, tmp_su3_0);

        _su3_plus_su3(g_stout_Gamma_field[x][mu], g_stout_Gamma_field[x][mu], tmp_su3_0);

        /*
         *  f_2 * Q * U  *\Sigma' + f_2 * U * \Sigma' * Q
         */
        _su3_times_su3(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][x][mu], g_previous_stout_force_field[x][mu]);
        _su3_times_su3(tmp_su3_1, g_Q_smear_iterations[stout_iter_counter][x][mu], tmp_su3_0);

        _su3_times_su3(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][x][mu],  g_Q_smear_iterations[stout_iter_counter][x][mu]);
        _su3_times_su3(tmp_su3_2, g_previous_stout_force_field[x][mu], tmp_su3_0);
        _su3_plus_su3(tmp_su3_0, tmp_su3_1, tmp_su3_2);

        _complex_times_su3(tmp_su3_0, f_2, tmp_su3_0);

        _su3_plus_su3(g_stout_Gamma_field[x][mu], g_stout_Gamma_field[x][mu], tmp_su3_0);

        /*
         *  now we can calculate \Lambda eqtn (73) hep-lat/0311018
         */

        /*
         *  0.5* (\Gamma + \Gamma^\dagger)
         *
         *  maybe this should be a macro
         */
        g_stout_Lambda_field[x][mu].c00.re = g_stout_Gamma_field[x][mu].c00.re;
        g_stout_Lambda_field[x][mu].c00.im = 0.0;
        g_stout_Lambda_field[x][mu].c11.re = g_stout_Gamma_field[x][mu].c11.re;
        g_stout_Lambda_field[x][mu].c11.im = 0.0;
        g_stout_Lambda_field[x][mu].c22.re = g_stout_Gamma_field[x][mu].c22.re;
        g_stout_Lambda_field[x][mu].c22.im = 0.0;

        g_stout_Lambda_field[x][mu].c01.re =
          0.5*(g_stout_Gamma_field[x][mu].c01.re * g_stout_Gamma_field[x][mu].c10.re
              +g_stout_Gamma_field[x][mu].c01.im * g_stout_Gamma_field[x][mu].c10.im);
        g_stout_Lambda_field[x][mu].c10.re = g_stout_Lambda_field[x][mu].c01.re;
        g_stout_Lambda_field[x][mu].c01.im = 
          0.5*(g_stout_Gamma_field[x][mu].c01.im * g_stout_Gamma_field[x][mu].c10.re
              -g_stout_Gamma_field[x][mu].c01.re * g_stout_Gamma_field[x][mu].c10.im);
        g_stout_Lambda_field[x][mu].c10.im = -g_stout_Lambda_field[x][mu].c01.im;

        g_stout_Lambda_field[x][mu].c02.re = 
          0.5*(g_stout_Gamma_field[x][mu].c02.re * g_stout_Gamma_field[x][mu].c20.re
              +g_stout_Gamma_field[x][mu].c02.im * g_stout_Gamma_field[x][mu].c20.im);
        g_stout_Lambda_field[x][mu].c20.re = g_stout_Lambda_field[x][mu].c02.re;
        g_stout_Lambda_field[x][mu].c02.im = 
          0.5*(g_stout_Gamma_field[x][mu].c02.im * g_stout_Gamma_field[x][mu].c20.re
              -g_stout_Gamma_field[x][mu].c02.re * g_stout_Gamma_field[x][mu].c20.im);
        g_stout_Lambda_field[x][mu].c20.im = -g_stout_Lambda_field[x][mu].c02.im;

        g_stout_Lambda_field[x][mu].c12.re = 
          0.5*(g_stout_Gamma_field[x][mu].c12.re * g_stout_Gamma_field[x][mu].c21.re
              +g_stout_Gamma_field[x][mu].c12.im * g_stout_Gamma_field[x][mu].c21.im);
        g_stout_Lambda_field[x][mu].c21.re = g_stout_Lambda_field[x][mu].c12.re;
        g_stout_Lambda_field[x][mu].c12.im = 
          0.5*(g_stout_Gamma_field[x][mu].c12.im * g_stout_Gamma_field[x][mu].c21.re
              -g_stout_Gamma_field[x][mu].c12.re * g_stout_Gamma_field[x][mu].c21.im);
        g_stout_Lambda_field[x][mu].c21.im = -g_stout_Lambda_field[x][mu].c12.im;

        /*
         *  1/(2*N) * Tr( (\Gamma + \Gamma^\dagger)
         */
        tmp_c_0.re = 0.33333333333333333333333333333333333333*
          (g_stout_Lambda_field[x][mu].c00.re + g_stout_Lambda_field[x][mu].c11.re + g_stout_Lambda_field[x][mu].c22.re);
        tmp_c_0.im = 0.33333333333333333333333333333333333333*
          (g_stout_Lambda_field[x][mu].c00.im + g_stout_Lambda_field[x][mu].c11.im + g_stout_Lambda_field[x][mu].c22.im);
        /*
         *  now we add both contributions
         */
        g_stout_Lambda_field[x][mu].c00.re += tmp_c_0.re; 
        g_stout_Lambda_field[x][mu].c00.im += tmp_c_0.im; 
        g_stout_Lambda_field[x][mu].c11.re += tmp_c_0.re; 
        g_stout_Lambda_field[x][mu].c11.im += tmp_c_0.im; 
        g_stout_Lambda_field[x][mu].c22.re += tmp_c_0.re; 
        g_stout_Lambda_field[x][mu].c22.im += tmp_c_0.im; 





        /*######################*/
        /*int a,b,c;
          printf("0000000\n");

          for(a=0; a< VOLUME*dim; a++)
          _su3_zero(stout_force_field[a]);
          _su3_one(stout_force_field[14]);
          _su3_one(stout_force_field[0]);
          printf("AAAAAA\n");
          print_su3(&(stout_force_field[14]));
          printf("BBBBBB\n");
          print_su3(&(stout_force_field[12]));
          a=5;
          b=3;
          printf("CCCCCC\n");
          printf("\n");
          print_su3(&(g_stout_force_field[a][b]));
          printf("\n");
          print_su3(&(g_stout_force_field[0][0]));
          printf("\n");
          print_su3(&(g_stout_force_field[a][b-1]));
          printf("\n");
          print_su3(&(g_stout_force_field[a+1][b]));

          for(a=0; a< VOLUME*dim; a++)
          _su3_zero(stout_force_field[a]);
          _su3_one(stout_force_field[14]);
          printf("AAAAAA\n");
          print_su3(&(stout_force_field[14]));
          printf("BBBBBB\n");
          print_su3(&(stout_force_field[12]));
          a=4;
          b=3;
          _su3_one(g_previous_stout_force_field[a][b]);
          print_su3(&(g_stout_force_field[a][b]));*/
        /*_su3_one(g_previous_stout_force_field[a][b]);
          print_su3(&(g_previous_stout_force_field[a][b]));*/
        /*print_su3(&(g_previous_stout_force_field[a][b]));

          exit(54);*/
        /*######################*/


        _trace_su3_times_su3(tmp_c_0, tmp_su3_0, tmp_su3_0);
        printf("AAAAAA\n");
        _trace_su3_times_su3(tmp_c_0, g_previous_stout_force_field[x][mu], tmp_su3_0);
        printf("CCCCCC\n");
        printf("tmp_c_0.re = %f tmp_c_0.im = %f\n", tmp_c_0.re, tmp_c_0.im);

        /*printf("tmp    = %f +i* %f \n", tmp_c_1.re, tmp_c_1.im);
          printf("c0_max = %f +i* %f \n", c0_max.re, c0_max.im);

          printf("c0 = %f +i* %f \n", c0.re, c0.im);
          printf("c1 = %f +i* %f \n", c1.re, c1.im);
          printf("c0_max = %f +i* %f \n", c0_max.re, c0_max.im);
          printf("theta  = %f \n", theta);
          printf("u      = %f \n", u);
          printf("w      = %f \n", w);*/
      }


    /*
     *  now we can calculate the forcefield \Signma_\mu(x)
     *  eqtn (75) in hep-lat/0311018
     */
    for(x = 0; x < VOLUME; x++)
      for(mu = 0; mu < dim; mu++)
      {

        /*
         *  f_0 * I + f_1 * Q + f_2 * Q^2
         */
        _su3_times_su3(tmp_su3_0, g_Q_smear_iterations[stout_iter_counter][x][mu], g_Q_smear_iterations[stout_iter_counter][x][mu]);
        _complex_times_su3(tmp_su3_1, f_2, tmp_su3_0);

        _complex_times_su3(tmp_su3_0, f_1, g_Q_smear_iterations[stout_iter_counter][x][mu]);
        _su3_plus_su3(tmp_su3_0, tmp_su3_0, tmp_su3_1);

        tmp_su3_0.c00.re += f_0.re;
        tmp_su3_0.c00.im += f_0.im;
        tmp_su3_0.c11.re += f_0.re;
        tmp_su3_0.c11.im += f_0.im;
        tmp_su3_0.c22.re += f_0.re;
        tmp_su3_0.c22.im += f_0.im;

        /*
         *  \Sigma' * (above)
         */
        _su3_times_su3(g_stout_force_field[x][mu], g_previous_stout_force_field[x][mu], tmp_su3_0);

        /*
         *  i * (C_\mu)^\dagger * Lambda 
         */
        _su3d_times_su3(tmp_su3_0, g_C_smear_iterations[stout_iter_counter][x][mu], g_stout_Lambda_field[x][mu]);
        _itimes_su3(tmp_su3_1, tmp_su3_0);
        _su3_plus_su3(g_stout_force_field[x][mu], g_stout_force_field[x][mu], tmp_su3_1);

        /*
         *  no we do the sum eqtn (75) hep-lat/0311018
         *  since our rho is mu/nu independent we pull in front of the sum
         */
        for(nu=0; nu < dim; nu++)
          if(nu != mu)
          {
            _su3d_times_su3(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][x][nu], g_stout_Lambda_field[x][nu]);
            _su3d_times_su3(tmp_su3_1, g_gauge_field_smear_iterations[stout_iter_counter][g_iup[x][nu]][mu], tmp_su3_0);
            _su3_times_su3(tmp_su3_2, g_gauge_field_smear_iterations[stout_iter_counter][g_iup[x][mu]][nu], tmp_su3_1);


            _su3_times_su3(tmp_su3_0, g_stout_Lambda_field[g_idn[x][nu]][mu], g_gauge_field_smear_iterations[stout_iter_counter][g_idn[x][nu]][nu]);
            _su3d_times_su3(tmp_su3_1, g_gauge_field_smear_iterations[stout_iter_counter][g_idn[x][nu]][mu], tmp_su3_0);
            _su3d_times_su3(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][g_iup[g_idn[x][nu]][mu]][nu], tmp_su3_1);
            _su3_plus_su3(tmp_su3_2, tmp_su3_2, tmp_su3_0);


            _su3d_times_su3(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][g_idn[x][nu]][mu], g_gauge_field_smear_iterations[stout_iter_counter][g_idn[x][nu]][nu]);
            _su3_times_su3(tmp_su3_1, g_stout_Lambda_field[g_iup[g_idn[x][nu]][mu]][nu], tmp_su3_0);
            _su3d_times_su3(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][g_iup[g_idn[x][nu]][mu]][nu], tmp_su3_1);
            _su3_plus_su3(tmp_su3_2, tmp_su3_2, tmp_su3_0);

            _su3_times_su3(tmp_su3_0, g_stout_Lambda_field[g_idn[x][nu]][nu], g_gauge_field_smear_iterations[stout_iter_counter][g_idn[x][nu]][nu]);
            _su3d_times_su3(tmp_su3_1, g_gauge_field_smear_iterations[stout_iter_counter][g_idn[x][nu]][mu], tmp_su3_0);
            _su3d_times_su3(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][g_iup[g_idn[x][nu]][mu]][nu], tmp_su3_1);
            _su3_minus_su3(tmp_su3_2, tmp_su3_2, tmp_su3_0);


            _su3_dagger(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][x][nu]);
            _su3d_times_su3(tmp_su3_1, g_gauge_field_smear_iterations[stout_iter_counter][g_idn[x][nu]][mu], tmp_su3_0);
            _su3_times_su3(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][g_iup[x][mu]][nu], tmp_su3_1);
            _su3_times_su3(tmp_su3_1, g_stout_Lambda_field[g_iup[x][mu]][nu], tmp_su3_0);
            _su3_minus_su3(tmp_su3_2, tmp_su3_2, tmp_su3_1);


            _su3_times_su3d(tmp_su3_0, g_stout_Lambda_field[g_iup[x][nu]][mu], g_gauge_field_smear_iterations[stout_iter_counter][x][nu]);
            _su3d_times_su3(tmp_su3_1, g_gauge_field_smear_iterations[stout_iter_counter][g_iup[x][nu]][mu], tmp_su3_0);
            _su3_times_su3(tmp_su3_0, g_gauge_field_smear_iterations[stout_iter_counter][g_iup[g_iup[x][mu]][nu]][nu], tmp_su3_1);
            _su3_minus_su3(tmp_su3_2, tmp_su3_2, tmp_su3_0);
          }
        
        /*
         *  now we multiply with i * rho and subtract from the first part
         */
        tmp_c_0.re = 0;
        tmp_c_0.im = stout_rho;

        _complex_times_su3(tmp_su3_0, tmp_c_0, tmp_su3_2);
        _su3_minus_su3(g_stout_force_field[x][mu], g_stout_force_field[x][mu], tmp_su3_0);
      }

    /*
     *  finally we save the \Sigma field, because we need it as \Sigma'
     *  in the next iteration
     */
    for(x = 0; x < VOLUME; x++)
      for(mu = 0; mu < dim; mu++)
      {
        _su3_assign(g_previous_stout_force_field[x][mu], g_stout_force_field[x][mu]);
      }
  }
}
