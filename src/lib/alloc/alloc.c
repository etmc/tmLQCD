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


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>


/**
 * @brief      Safe malloc implementation that checks malloc for NULL. Never use
 *             this function instead use the macro safe_malloc()
 *
 * @param[in]  size  The allocation size in bytes
 * @param[in]  file  __FILE__
 * @param[in]  line  __LINE__
 * @param[in]  func  __func__
 *
 * @return     Pointer returned by malloc()
 */
void *safe_malloc_impl(size_t size, const char *file, int line, const char *func)
{
    if (size <= 0) {
        fprintf(stderr, "safe_malloc: zero-size allocation at %s:%d (%s)\n", file, line, func);
        abort();
    }

    void *p = malloc(size);
    if (p == NULL) {
        fprintf(stderr, "safe_malloc: failed to allocate %zu bytes at %s:%d (%s): %s\n",
                size, file, line, func, strerror(errno));
        abort();
    }
    return p;
}


/**
 * @brief      Identical to safe_malloc_impl above just that is call calloc
 *             instead of malloc. Never use this function instead use the macro
 *             safe_calloc()
 */
void *safe_calloc_impl(size_t size, const char *file, int line, const char *func)
{
   if (size <= 0) {
      fprintf(stderr, "safe_calloc: zero-size allocation at %s:%d (%s)\n", file, line, func);
      abort();
   }

   void *p = calloc(size, 1);
   if (p == NULL) {
      fprintf(stderr, "safe_calloc: failed to allocate %zu bytes at %s:%d (%s): %s\n",
              size, file, line, func, strerror(errno));
      abort();
   }
   return p;
}