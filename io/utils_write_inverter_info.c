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

#include "spinor.ih"

void write_inverter_info(WRITER * writer, paramsInverterInfo const *info)
{
  char *message;
  n_uint64_t bytes;
  message = (char*)malloc(1024);

  if (info->mms) {
    sprintf(message, "solver = %s\n"
                     "result is for Q^dagger Q!\n"
                     "multiple mass solver\n"
                     "epssq = %e\n"
                     "noiter = %d\n"
                     "kappa = %.12f, inverted mu = %.12f, lowest mu = %.12f\n"
                     "inverter version = %s\n"
                     "date = %s",
                     info->inverter,
                     info->epssq, info->iter, info->kappa,
                     info->cgmms_mass,
                     info->mu, info->package_version,
                     info->date);
  }
  else {
    if (!info->heavy) {
      sprintf(message, "solver = %s\n"
                       "epssq = %e\n"
                       "noiter = %d\n"
                       "kappa = %.12f, mu = %.12f\n"
                       "inverter version = %s\n"
                       "date = %s",
                       info->inverter,
                       info->epssq, info->iter, info->kappa, info->mu,
                       info->package_version, info->date);
    }
    else {
      sprintf(message, "solver = %s\n"
                       "epssq = %e\n"
                       "noiter = %d\n"
                       "kappa = %.12f, mubar = %.12f, epsbar=%.12f\n"
                       "inverter version = %s\n"
                       "date = %s",
                       info->inverter,
                       info->epssq, info->iter, info->kappa, info->mubar,
                       info->epsbar,
                       info->package_version, info->date);
    }
  }
  bytes = strlen(message);
  write_header(writer, 1, 0, "inverter-info", bytes);
  write_message(writer, message, bytes);
  close_writer_record(writer);

  free(message);
  return;
}

