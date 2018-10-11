/***********************************************************************
 *  
 * Copyright (C) 2018 Bartosz Kostrzewa
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

#include "global.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

void tm_debug_printf(const int proc_id,
                     const int dbg_level_threshold,
                     const char * format,
                     ...)
{
  if( g_proc_id == proc_id && g_debug_level >= dbg_level_threshold ){
    va_list arglist;
    va_start(arglist, format);
    vprintf(format, arglist);
    va_end(arglist);
  }
}

