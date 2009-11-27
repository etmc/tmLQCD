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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <stdlib.h>
#include <stdio.h>

#include <io/selector.h>
#include <io/params.h>
#include <io/dml.h>


#ifndef isnan
# define isnan(x)						 \
  (sizeof (x) == sizeof (long double) ? isnan_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isnan_d (x)		 \
   : isnan_f (x))

#endif

/* These are factory functions, since the constructors for c-lime and lemon are different
   and they need different ways of opening files. Moving this to utility functions unclutters
   the main code, since we don't need additional #ifdefs anymore.
   Since lemon is collective and c-lime isn't, some care needs to be taken. For now, a writer
   will only be constructed on the node with g_cart_id == 0 for c-lime, while a reader will
   be created everywhere to exploit trivial parallellization. Be careful not to call
   construct_writer and friends from within a "if (node_num == 0)" type statement, because
   it will cause lemon to get deadlocked! */
void construct_writer(WRITER ** writer, char const *filename);
void destruct_writer(WRITER * writer);

void construct_reader(READER ** reader, char const *filename);
void destruct_reader(READER * reader);

void kill_with_error(LIME_FILE *fh, int const rank, char const *error);

int read_message(READER *reader, char **buffer);
int write_message(WRITER * writer, char const *buffer, uint64_t bytes);
void write_header(WRITER * writer, int MB, int ME, char const *type, uint64_t bytes);

void write_checksum(WRITER *writer, DML_Checksum const *checksum, char const *name);
void write_xlf_info(WRITER *writer, paramsXlfInfo const *info);
void write_inverter_info(WRITER * writer, paramsInverterInfo const *info);

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
