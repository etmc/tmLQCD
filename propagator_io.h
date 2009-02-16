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

#ifndef _PROPAGATOR_IO_H
#define _PROPAGATOR_IO_H

/* write a one flavour propagator to file */
int write_propagator(spinor * const s, spinor * const r, char * filename, 
		     const int append, const int prec, const int format);

/* write a one flavour source to file */
int write_source(spinor * const s, spinor * const r, char * filename, 
		     const int append, const int prec);

int read_source(spinor * const s, spinor * const r, char * filename, 
		const int format, const int position);

/* write a two flavour propagator to file */
int write_double_propagator(spinor * const s, spinor * const r, 
			    spinor * const p, spinor * const q,
			    char * filename, const int append, const int prec);

int write_propagator_format(char * filename, const int prec, const int no_flavours);
int write_source_format(char * filename, const int prec, const int no_flavours,
			const int t, const int lx, const int ly, const int lz,
			const int is, const int ic);

int write_propagator_type(const int type, char * filename);
int write_source_type(const int type, char * filename);

int get_propagator_type(char * filename);
int get_source_type(char * filename);

/* The binary IO */
int write_lime_spinor(spinor * const s, spinor * const r, char * filename, const int append, const int prec);
int read_lime_spinor(spinor * const s, spinor * const r, char * filename, const int position);

#endif
