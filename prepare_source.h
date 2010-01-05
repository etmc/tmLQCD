/***********************************************************************
 * $Id$
 *
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

#ifndef _PREPARE_SOURCE_H
#define _PREPARE_SOURCE_H

void prepare_source(const int nstore, const int isample, const int ix, const int op_id, 
		    const int read_source_flag, const int source_splitted, 
		    const int source_location, const int source_time_slice,
		    const int source_type, const int propagator_splitted,
		    char * source_input_filename);

#endif
