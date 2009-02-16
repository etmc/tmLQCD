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

/****************************************************
 * IO routines:
 *
 * write_gauge_field_time_p(char * filename)
 *   writes gauge field configuration to file
 *   with name filename.
 *
 * read_gauge_field_time_p(char * filename)
 *   reads gauge field configuration from file
 *   with name filename.
 *
 * Autor: Ines Wetzorke <wetzorke@ifh.de>
 *        Carsten Urbach <urbach@ifh.de>
 *
 ****************************************************/

#ifndef _IO_H
#define _IO_H

int write_gauge_field_time_p(char * filename);
int read_gauge_field_time_p(char * filename);

int write_spinorfield_eo_time_p(spinor * const s, spinor * const r, char * filename, const int append);
int read_spinorfield_eo_time(spinor * const s, spinor * const r, char * filename);
void write_su3(su3 * up, FILE * f);
int write_lime_gauge_field_old(char * filename, const double plaq, const int counter);
int write_lime_gauge_field_singleprec(char * filename, const double plaq, const int counter);
int read_lime_gauge_field_old(char * filename);
int read_lime_gauge_field_singleprec(char * filename);
int read_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename, 
			       const int ts, const int vol);
int read_spinorfield_cm_swap_single(spinor * const s, spinor * const r, char * filename,
                                    const int ts, const int vol);
int write_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename);

int write_eospinor(spinor * const s, char * filename, 
		   const double evalue, const double prec, const int nstore);
int read_eospinor(spinor * const s, char * filename);

int write_first_messages(FILE * parameterfile, const int integtyp, const int inv);

int write_rlxd_state(char * filename, int * const _state, const int rlxdsize);
int read_rlxd_state(char * filename, int * state, const int rlxdsize);

#endif
