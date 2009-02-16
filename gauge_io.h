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

#ifndef _GAUGE_IO_H
#define _GAUGE_IO_H

#include"dml.h"

/* int write_binary_gauge_data(LimeWriter * limewriter, */
/* 			    const int prec, DML_Checksum * ans); */
/* int read_binary_gauge_data(LimeReader * limereader,  */
/* 			   const double prec, DML_Checksum * ans); */

int write_lime_gauge_field(char * filename, const double plaq, const int counter, const int prec);
int read_lime_gauge_field(char * filename);

#endif
