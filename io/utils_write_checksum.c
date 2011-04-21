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

void write_checksum(WRITER * writer, DML_Checksum const *checksum, char const *name)
{
  char *message;
  uint64_t bytes;
  message = (char*)malloc(512);
  if (message == (char*)NULL) {
    kill_with_error(writer->fp, g_cart_id, "Memory allocation error in write_checksum. Aborting\n");
  }
  sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                   "<scidacChecksum>\n"
                   "  <version>1.0</version>\n"
                   "  <suma>%08x</suma>\n"
                   "  <sumb>%08x</sumb>\n"
                   "</scidacChecksum>", checksum->suma, checksum->sumb);
  bytes = strlen(message);
  /* The message begin bit is 0, because this is written as part of a data message
   * the end bit is 1, since this should be the last record of a message */
  if (name == NULL)
      write_header(writer, 0, 1, "scidac-checksum", bytes);
  else
      write_header(writer, 0, 1, name, bytes);

  write_message(writer, message, bytes);

  close_writer_record(writer);
  free(message);
  return;
}
