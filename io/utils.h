/***********************************************************************
 * $Id$
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

#ifndef _UTILS_H
#define _UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <dml.h>
#include <io/params.h>
#include <lime.h>

#ifndef isnan
# define isnan(x)						 \
  (sizeof (x) == sizeof (long double) ? isnan_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isnan_d (x)		 \
   : isnan_f (x))

#endif



#ifdef HAVE_LIBLEMON
# include <lemon.h>
#endif /* HAVE_LIBLEMON */

#ifdef HAVE_LIBLEMON
void kill_with_error(MPI_File *fh, int const rank, char const *error);
void read_message_parallel(LemonReader *lemonreader, char **buffer);
void write_message_parallel(LemonWriter *lemonwriter, char const *buffer, uint64_t bytes);
void write_header_parallel(LemonWriter *lemonwriter, int MB, int ME, char *type, uint64_t bytes);

void write_checksum_parallel(LemonWriter *lemonwriter, DML_Checksum const *checksum);
void write_xlf_info_parallel(LemonWriter *lemonwriter, paramsXlfInfo const *info);
#else
void kill_with_error(FILE *fh, int const rank, char const *error);
#endif /* HAVE_LIBLEMON */

int read_message(LimeReader *limereader, char **buffer);
int write_message(LimeWriter * limewriter, char const *buffer, uint64_t bytes);
void write_header(LimeWriter *limewriter, int MB, int ME, char *type, uint64_t bytes);

void write_checksum(LimeWriter *limewriter, DML_Checksum const *checksum);
void write_xlf_info(LimeWriter *limewriter, paramsXlfInfo const *info);
void write_inverter_info(LimeWriter * writer,
			 paramsInverterInfo const *info);

void engineering(char *result, double value, char const *units);
int parse_checksum_xml(char *message, DML_Checksum *checksum);

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


#endif
