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

#include "utils.ih"

void write_xlf_info(WRITER * writer, paramsXlfInfo const *info)
{
  char *message;
  uint64_t bytes;

  message = (char*)malloc(512);
  if (message == (char*)NULL)
    kill_with_error(writer->fp, g_cart_id, "Memory allocation error in write_xlf_info. Aborting\n");

  if (info->kappa != 0.0) {
    sprintf(message, "plaquette = %14.12f\n"
                     " trajectory nr = %d\n"
                     " beta = %f, kappa = %f, mu = %f, c2_rec = %f\n"
                     " time = %ld\n"
                     " hmcversion = %s\n"
                     " mubar = %f\n"
                     " epsilonbar = %f\n"
                     " date = %s",
                     info->plaq, info->counter, info->beta, info->kappa,
                     info->mu, info->c2_rec, info->time, info->package_version,
                     info->mubar, info->epsilonbar, info->date);
  }
  else {
    sprintf(message, "plaquette = %e\n"
                     " trajectory nr = %d\n"
                     " beta = %f\n"
                     " kappa = %f\n"
                     " 2*kappa*mu = %f\n"
                     " c2_rec = %f\n"
                     " date = %s",
                     info->plaq, info->counter, info->beta, info->kappa,
                     info->mu, info->c2_rec, info->date);
  }
  bytes = strlen(message);

  write_header(writer, 1, 1, "xlf-info", bytes);
  write_message(writer, message, bytes);

  close_writer_record(writer);

  free(message);
  return;
}
