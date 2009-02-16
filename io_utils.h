/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
/* $Id$ */

#ifndef _IO_UTILS_H
#define _IO_UTILS_H

#include"dml.h"

#ifndef isnan
# define isnan(x)						 \
  (sizeof (x) == sizeof (long double) ? isnan_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isnan_d (x)		 \
   : isnan_f (x))

#endif

int isnan_f  (float       x);
int isnan_d  (double      x);
int isnan_ld (long double x);

int write_checksum(char * filename, DML_Checksum * checksum);
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
int write_xlf_info(const double plaq, const int counter, char * filename, const int append);
int write_inverter_info(const double epssq, const int iter, const int heavy, 
			const int append, char * filename);

#endif
