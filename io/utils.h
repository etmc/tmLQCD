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

#include <dml.h>

#ifdef HAVE_LIBLEMON
# include <lemon.h>
#endif /* HAVE_LIBLEMON */

#ifdef HAVE_LIBLEMON
void kill_with_error(MPI_File *fh, int const rank, char const *error);
void read_message_parallel(LemonReader * lemonreader, char **buffer);
void write_message_parallel(LemonWriter * lemonwriter, char *buffer, uint64_t bytes);
void write_header_parallel(LemonWriter * lemonwriter, int MB, int ME, char *type, uint64_t bytes);

void write_checksum_parallel(LemonWriter * lemonwriter, DML_Checksum * checksum);
void write_xlf_info_parallel(LemonWriter * lemonwriter, const double plaq, const int counter);
#endif /* HAVE_LIBLEMON */

void engineering(char *result, double value, char const *units);
void parse_checksum_xml(char *message, DML_Checksum * checksum);

#endif
