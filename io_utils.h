/* $Id$ */

#ifndef _IO_UTILS_H
#define _IO_UTILS_H

#ifndef isnan
# define isnan(x)						 \
  (sizeof (x) == sizeof (long double) ? isnan_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isnan_d (x)		 \
   : isnan_f (x))

#endif

int isnan_f  (float       x);
int isnan_d  (double      x);
int isnan_ld (long double x);

void byte_swap(void *ptr, int nmemb);
void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb);
void byte_swap_assign_single2double(void * out_ptr, void * in_ptr, int nmemb);
void single2double(void * out_ptr, void * in_ptr, int nmemb);
void byte_swap_assign_double2single(void * out_ptr, void * in_ptr, int nmemb);
void double2single(void * out_ptr, void * in_ptr, int nmemb);
int big_endian();
int write_ildg_format_xml(char *filename, LimeWriter * limewriter, const int precision);
void single2double_cm(spinor * const R, float * const S);
void double2single_cm(float * const S, spinor * const R);
void zero_spinor(spinor * const R);
int write_xlf_info(const double plaq, const int counter, char * filename);

#endif
