/***********************************************************************
 *
 * Copyright (C) 2026 Roman Gruber
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
 * Allocation utils
 *
 * Author: Roman Gruber
 *         roman.gruber@unibe.ch
 *
 *******************************************************************************/

#ifndef ALLOC_H
#define ALLOC_H


#include <stdbool.h>


#ifdef __cplusplus
extern "C" {
#endif


void *safe_malloc_impl(size_t size, const char *file, int line, const char *func);
void *safe_calloc_impl(size_t size, const char *file, int line, const char *func);
#define safe_malloc(size) safe_malloc_impl((size), __FILE__, __LINE__, __func__)
#define safe_calloc(size) safe_calloc_impl((size), __FILE__, __LINE__, __func__)


#ifdef __cplusplus
}
#endif

#endif
