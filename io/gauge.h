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

#ifndef _GAUGE_H
#define _GAUGE_H

#if HAVE_CONFIG_H
#include<config.h>
#endif

#include <io/selector.h>
#include <io/params.h>
#include <io/utils.h>


int read_gauge_field(char *filename, su3 ** const gf);
int read_binary_gauge_data(READER *reader, DML_Checksum *checksum, paramsIldgFormat * ildgformat, su3 ** const gf);

int write_gauge_field(char * filename, int prec, paramsXlfInfo const *xlfInfo);
int write_binary_gauge_data(WRITER * writer, const int prec, DML_Checksum * checksum);

void write_ildg_format(WRITER *writer, paramsIldgFormat const *format);

#endif
