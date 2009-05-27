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
int write_ildg_format_parallel(LemonWriter *writer, const int prec)
{
  uint64_t bytes;
  int status = 0;
  LemonRecordHeader *header;
  char *buf;

  buf = (char*)malloc(512);
  sprintf(buf, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
               "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n"
               "            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
               "            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n"
               "  <version> 1.0 </version>\n"
               "  <field> su3gauge </field>\n"
               "  <precision> %d </precision>\n"
               "  <lx> %d </lx>\n"
               "  <ly> %d </ly>\n"
               "  <lz> %d </lz>\n"
               "  <lt> %d </lt>\n"
               "</ildgFormat>", prec, LX*g_nproc_x, LY*g_nproc_y, LZ*g_nproc_z, T*g_nproc_t);

  bytes = strlen(buf);
  header = lemonCreateHeader(1, 0, "ildg-format", bytes);
  if(header == (LemonRecordHeader*)NULL)
  {
    fprintf(stderr, "LEMON create header ildg-format error\nPanic! Aborting...\n");
    exit(500);
  }
  status = lemonWriteRecordHeader(header, writer);
  if(status < 0 )
  {
    fprintf(stderr, "LEMON write header ildg-format error %d\nPanic! Aborting...\n", status);
    exit(500);
  }
  lemonDestroyHeader(header);
  lemonWriteRecordData(buf, &bytes, writer);
  lemonWriterCloseRecord(writer);
  free(buf);
  return 0;
}
