/* $Id$ */

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
