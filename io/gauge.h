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

#include "lime.h"
#if HAVE_CONFIG_H
#include<config.h>
#endif

#ifdef HAVE_LIBLEMON
# include <lemon.h>
#endif /* HAVE_LIBLEMON */

#include <io/utils.h>

#ifdef HAVE_LIBLEMON
#  define LIME_FILE MPI_File
#  define WRITER LemonWriter
#  define READER LemonReader
#  define CreateReader lemonCreateReader
#  define ReaderNextRecord lemonReaderNextRecord
#  define ReaderType lemonReaderType
#  define ReaderCloseRecord lemonReaderCloseRecord
#  define DestroyReader lemonDestroyReader
#else /* HAVE_LIBLEMON */
#  define LIME_FILE FILE
#  define WRITER LimeWriter
#  define READER LimeReader
#  define CreateReader limeCreateReader
#  define ReaderNextRecord limeReaderNextRecord
#  define ReaderType limeReaderType
#  define ReaderCloseRecord limeReaderCloseRecord
#  define DestroyReader limeDestroyReader
#endif

void read_gauge_field(char *filename, DML_Checksum *scidac_checksum, char **xlf_info, char **ildg_data_lfn);
void read_binary_gauge_data(READER *reader, DML_Checksum *checksum);

void write_gauge_field(char * filename, int prec, paramsXlfInfo const *xlfInfo);
void write_binary_gauge_data(WRITER *writer, const int prec, DML_Checksum * checksum);

void write_ildg_format(WRITER *writer, paramsIldgFormat const *format);
