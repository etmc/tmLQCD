/***********************************************************************
* Copyright (C) 2011 Siebren Reker
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

void write_xlf_info_xml(WRITER * writer, paramsXlfInfo const *info)
{
  char *message;
  uint64_t bytes;

  message = (char*)malloc(512);
  if (message == (char*)NULL)
    kill_with_error(writer->fp, g_cart_id, "Memory allocation error in write_xlf_info_xml. Aborting\n");

  if (info->kappa != 0.0) {
    sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<xlf-info>\n"
        "  <plaquette>%14.12f</plaquette>\n"
        "  <trajectory>%d</trajectory>\n"
        "  <beta>%f</beta>\n"
        "  <kappa>%f</kappa>\n"
        "  <mu>%f</mu>\n"
        "  <c2_rec>%f</c2_rec>\n"
        "  <time>%ld</time>\n"
        "  <hmcversion>%s</hmcversion>\n"
        "  <mubar>%f</mubar>\n"
        "  <epsilonbar>%f</epsilonbar>\n"
        "  <date>%s</date>\n"
        "</xlf-info>", info->plaq, info->counter, info->beta, info->kappa,
                       info->mu, info->c2_rec, info->time, info->package_version,
                       info->mubar, info->epsilonbar, info->date);
  bytes = strlen(message);
  }
  else {
    sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<xlf-info>\n"
        "  <plaquette>%e</plaquette>\n"
        "  <trajectory>%d</trajectory>\n"
        "  <beta>%f</beta>\n"
        "  <kappa>%f</kappa>\n"
        "  <2kappamu>%f</2kappamu>\n"
        "  <c2_rec>%f</c2_rec>\n"
        "  <date>%s</date>\n"
        "</xlf-info>", info->plaq, info->counter, info->beta, info->kappa,
                       info->mu, info->c2_rec, info->date);
  }
  bytes = strlen(message);

  write_header(writer, 1, 1, "xlf-info", bytes);
  write_message(writer, message, bytes);

  close_writer_record(writer);

  free(message);
  return;
}
