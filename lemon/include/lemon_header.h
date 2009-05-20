#pragma once

/****************************************************************************
 * LEMON v0.99                                                              *
 *                                                                          *
 * This file is part of the LEMON implementation of the SCIDAC LEMON format. *
 *                                                                          *
 * It is based directly upon the original c-lemon implementation,            *
 * as maintained by C. deTar for the USQCD Collaboration,                   *
 * and inherits its license model and parts of its original code.           *
 *                                                                          *
 * LEMON is free software: you can redistribute it and/or modify            *
 * it under the terms of the GNU General Public License as published by     *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 *                                                                          *
 * LEMON is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 * GNU General Public License for more details. You should have received    *
 * a copy of the GNU General Public License along with LEMON. If not,       *
 * see <http://www.gnu.org/licenses/>.                                      *
 *                                                                          *
 * LEMON was written for the European Twisted Mass Collaboration.           *
 * For support requests or bug reports,                                     *
 * please contact A. Deuzeman (a.deuzeman@rug.nl)                           *
 ****************************************************************************/

typedef struct
{
  unsigned int lemon_version;
  int MB_flag;
  int ME_flag;
  char *type;
  uint64_t data_length;
} LemonRecordHeader;

/* Header manipulators */
LemonRecordHeader *lemonCreateHeader(int MB_flag, int ME_flag, char *type, uint64_t reclen);
void lemonDestroyHeader(LemonRecordHeader *h);
