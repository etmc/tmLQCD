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

#ifndef _PARALLEL_IO_H
#define _PARALLEL_IO_H

#include <c-lemon/include/lime.h>
#include"dml.h"

int read_lime_gauge_field_parallel(char *filename, int *ranlxd_state);
int read_lime_gauge_field_singleprec_parallel(char const * filename);

int read_binary_gauge_data_parallel(LimeReader * limereader, DML_Checksum * checksum);
int read_checksum_parallel(LimeReader * limereader, DML_Checksum * checksum);
int read_rlxd_state_parallel(LimeReader * limereader, int * state);

int write_lime_gauge_field_parallel(char * filename, const double plaq, const int counter, const int prec);

int write_binary_gauge_data_parallel(LimeWriter * limewriter, const int prec, DML_Checksum * ans);
int write_checksum_parallel(LimeWriter * limewriter, DML_Checksum * checksum);
int write_rlxd_state_parallel(LimeWriter * limewriter);
int write_xlf_info_parallel(LimeWriter * limewriter, const double plaq, const int counter);

#endif
