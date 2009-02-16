/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2007 Craig McNeile
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

#ifndef _STOUT_SMEAR_
#define _STOUT_SMEAR_

int stout_smear_gauge_field(const double rho , const int no_iters); 
void scale_su3(su3 *in, double scale);
void  project_anti_herm(su3 *omega);
void  print_su3(su3 *in);
void  print_su3_octave(su3 *in);
void  print_spinor(spinor *in); 
void print_config_to_screen(su3 **in);
void  print_scalar_complex_field_to_screen(complex *in_field);
void  print_scalar_real_field_to_screen(double *in_field); 
void load_config_from_file(su3 **in, char * filename);

#endif
