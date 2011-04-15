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

void write_source_format(WRITER *writer, paramsSourceFormat const *format)
{
  uint64_t bytes;
  char *buf = NULL;
#ifndef HAVE_LIBLEMON
  if(g_cart_id == 0) {
#endif /* ! HAVE_LIBLEMON */
  buf = (char*)malloc(512);
  sprintf(buf, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<etmcFormat>\n"
        "  <field>diracFermion</field>\n"
        "  <precision>%d</precision>\n"
        "  <flavours>%d</flavours>\n"
        "  <lx>%d</lx>\n"
        "  <ly>%d</ly>\n"
        "  <lz>%d</lz>\n"
        "  <lt>%d</lt>\n"
        "  <spin>%d</spin>\n"
        "  <colour>%d</colour>\n"
        "</etmcFormat>",
          format->prec, format->flavours,
          format->lx, format->ly, format->lz, format->lt,
          format->spins, format->colours);
  bytes = strlen(buf);
  /* This message should be preceded by inverter info
   * and followed by propagator format, so MB=ME=0 */
  write_header(writer, 0, 0, "etmc-source-format", bytes);
  write_message(writer, buf, bytes);
  close_writer_record(writer);

  free(buf);
#ifndef HAVE_LIBLEMON
  }
#endif /* ! HAVE_LIBLEMON */
}
