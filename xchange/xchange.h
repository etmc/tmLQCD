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
#ifndef _XCHANGE_H
#define _XCHANGE_H

#ifdef HAVE_CONFIG_H
#include "tmlqcd_config.h"
#endif

#include "xchange/xchange_field.h"
#include "xchange/xchange_gauge.h"
#include "xchange/xchange_deri.h"
#include "xchange/xchange_halffield.h"
#include "xchange/xchange_jacobi.h"
#include "xchange/xchange_2fields.h"
#include "xchange/xchange_lexicfield.h"

#  ifdef _USE_TSPLITPAR
#    include "xchange/xchange_field_tslice.h"
#  endif

#endif
