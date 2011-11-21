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

#include "gauge.ih"

void write_ildg_format(WRITER *writer, paramsIldgFormat const *format)
{
  uint64_t bytes;
  char *buf;

  buf = (char*)malloc(512);
  if (buf == (char*)NULL)
    kill_with_error(writer->fp, g_cart_id, "Memory allocation error in write_ildg_format. Aborting\n");

  sprintf(buf, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
          "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n"
          "            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
          "            xsi:schemaLocation=\"http://www.lqcd.org/ildg/filefmt.xsd\">\n"
          "  <version>1.0</version>\n"
          "  <field>su3gauge</field>\n"
          "  <precision>%d</precision>\n"
          "  <lx>%d</lx>\n"
          "  <ly>%d</ly>\n"
          "  <lz>%d</lz>\n"
          "  <lt>%d</lt>\n"
          "</ildgFormat>",
          format->prec, format->lx, format->ly, format->lz, format->lt);

  bytes = strlen(buf);
  write_header(writer, 1, 0, "ildg-format", bytes); /* ME is 0 because a ildg-binary-data record MUST follow */
  write_message(writer, buf, bytes);
  close_writer_record(writer);
  free(buf);
}
