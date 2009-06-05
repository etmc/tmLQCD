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

#ifndef _GAUGE_H
#define _GAUGE_H

#ifdef HAVE_LIBLEMON
# include <lemon.h>
#endif /* HAVE_LIBLEMON */

#ifdef HAVE_LIBLEMON
void read_lemon_gauge_field_parallel(char *filename, char **scidac_checksum, char **xlf_info, char **ildg_data_lfn);
void read_binary_gauge_data_parallel(LemonReader * lemonreader, DML_Checksum * checksum);

void write_lemon_gauge_field_parallel(char * filename, const double plaq, const int counter, const int prec);
void write_binary_gauge_data_parallel(LemonWriter * lemonwriter, const int prec, DML_Checksum * ans);
void write_ildg_format_parallel(LemonWriter *writer, const int prec);
#endif /* HAVE_LIBLEMON */

#endif
