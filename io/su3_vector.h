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

#ifndef _SU3VECTOR_H
#define _SU3VECTOR_H

#include <su3.h>

#include <io/selector.h>
#include <io/utils.h>

//int read_su3_vector(su3_vector * const s, su3_vector * const r, char * filename, const int position);
//int read_binary_su3_vector_data(su3_vector * const s, su3_vector * const r, READER * reader, DML_Checksum * checksum);
//int read_binary_su3_vector_data_l(su3_vector * const s, READER * reader, DML_Checksum * checksum);

//int write_su3_vector(WRITER * writer, su3_vector ** const s, su3_vector ** const r, const int flavours, const int prec);
int write_binary_su3_vector_data(su3_vector * const s, su3_vector * const r, WRITER * writer, DML_Checksum *checksum, int const prec, const int t0);

//void write_su3_vector_info(WRITER * writer, const int write_prop_format_flag, paramsInverterInfo * InverterInfo, int append);
//void write_su3_vector_format(WRITER *writer, paramsPropagatorFormat const *format);
//void write_su3_vector_type(WRITER *writer, const int type);

#endif
