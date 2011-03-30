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

int parse_ildgformat_xml(char *message, paramsIldgFormat *ildgformat)
{
  int read_prec = 0, read_lx = 0, read_ly = 0, read_lz = 0, read_lt = 0;
  char *pos = strtok(message, "<> \n\t");

  if (ildgformat == (paramsIldgFormat*)NULL) {
    return 0;
  }

  while (pos)
  {
    if (!strncmp(pos, "precision", 9)) {
      pos = strtok(0, "<> \n\t");
      sscanf(pos, "%d", &ildgformat->prec);
      read_prec = 1;
    }
    if (!strncmp(pos, "lx", 2)) {
      pos = strtok(0, "<> \n\t");
      sscanf(pos, "%d", &ildgformat->lx);
      read_lx = 1;
    }
    if (!strncmp(pos, "ly", 2)) {
      pos = strtok(0, "<> \n\t");
      sscanf(pos, "%d", &ildgformat->ly);
      read_ly = 1;
    }
    if (!strncmp(pos, "lz", 2)) {
      pos = strtok(0, "<> \n\t");
      sscanf(pos, "%d", &ildgformat->lz);
      read_lz = 1;
    }
    if (!strncmp(pos, "lt", 2)) {
      pos = strtok(0, "<> \n\t");
      sscanf(pos, "%d", &ildgformat->lt);
      read_lt = 1;
    }
    pos = strtok(0, "<> \n\t");
  }
  return (read_prec && read_lx && read_ly && read_lz && read_lt);
}
