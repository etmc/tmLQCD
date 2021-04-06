/***********************************************************************
 *
 * Copyright (C) 2018  Bartosz Kostrzewa
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
 *
 *******************************************************************************/

#ifndef INIT_CRITICAL_GLOBALS_H
#define INIT_CRITICAL_GLOBALS_H

#include "misc_types.h"

/* function to initialise global variables which always need to be defined
 * but which are not initialised by the input file reader, for example, but
 * rather via the command line or explicitly inside some of the 'main' functions 
 * This means that in order to initialise these to some defaults in case tmLQCD
 * is used as a library, it is necessary to initialise them explicitly using
 * this function, which is done in tmLQCD_invert_init (wrapper/lib_wrapper.c) */

void init_critical_globals(const tm_ProgramId_t program_id);

#endif
