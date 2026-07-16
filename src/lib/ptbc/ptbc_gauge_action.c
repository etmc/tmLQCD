/***********************************************************************
 *
 * Copyright (C) 2026 JingJing Li
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
 *******************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include "monomial/monomial.h"
#include "measure_rectangles.h"
#include "measure_gauge_action.h"
#include "global.h"
#include <ptbc.h>

#define err(test, ...) err_impl(test, __func__, __FILE__, __LINE__, __VA_ARGS__)
static void err_impl(const bool test, const char* func, const char* file, const int line, const char* format, ...)
{
    if (test) {
        va_list args;
        char message[1024];
        va_start(args, format);
        vsnprintf(message, 1024, format, args);
        va_end(args);
        char location[1024];
        snprintf(location, 1024, "%s:%d %s", file, line, func);
        fatal_error(message, location);
    }
}

/**
 * @brief find the id of monomial that corresponds to gauge
 * 
 * @return * int 
 */
static int gauge_mnl_id(void) {
  static int id = -1;
  if (id < 0) {
    for (int i = 0; i < no_monomials; i++) {
      if (monomial_list[i].type == GAUGE) { id = i; break; }
    }
    err(id < 0, "No GAUGE monomial found for PTBC swap action evaluation!");
  }
  return id;
}


/**
 * @brief calculate the gauge action
 * 
 * @return * double 
 */
static double ptbc_gauge_action(void) {
  monomial const *mnl = &monomial_list[gauge_mnl_id()];
  double s = g_beta * (mnl->c0 * measure_gauge_action((const su3 **)g_gauge_field, mnl->glambda));
  if (mnl->use_rectangles) {
    s += g_beta * (mnl->c1 * measure_rectangles((const su3 **)g_gauge_field));
  }
  return s;
}


/**
 * @brief   find gauge action difference of with swap partner defect and own defect
 * 
 * @param alt_inst    instance id of swap partner
 * @return * double 
 */
double ptbc_swap_dh(int const alt_inst) {
  int const own_inst = app()->ptbc.instance_id;
  double const s_own = ptbc_gauge_action();
  appm()->ptbc.instance_id = alt_inst;
  double const s_alt = ptbc_gauge_action();
  appm()->ptbc.instance_id = own_inst;
  return s_alt - s_own;
}
