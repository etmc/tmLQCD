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

int parse_checksum_xml(char *message, DML_Checksum *checksum)
{
  int  read_suma = 0, read_sumb = 0;
  char *pos = strtok(message, "<> \n\t");
 
  if (checksum == (DML_Checksum*)NULL)
    return 0;
  
  while (pos)
  {
    if (!strncmp(pos, "suma", 4))
    {
      pos = strtok(0, "<> \n\t");
      sscanf(pos, "%x", &checksum->suma);
      read_suma = 1;
    }
    if (!strncmp(pos, "sumb", 4))
    {
      pos = strtok(0, "<> \n\t");
      sscanf(pos, "%x", &checksum->sumb);
      read_sumb = 1;
    }
    pos = strtok(0, "<> \n\t");
  }
  return (read_suma && read_sumb);
}
