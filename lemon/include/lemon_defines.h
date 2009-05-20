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

#define LEMON_MAGIC_NO                   0x456789ab
#define LEMON_VERSION                             1
#define LEMON_MAX_BUFSIZE                2147483647
#define HDR_SIZE_BYTES                         144
#define MB_MASK                         ((int)0x80)
#define ME_MASK                         ((int)0x40)
#define MAX_LEMON_HDR_REC_TYPE                  128

enum LEMON_Error_codes
{
  LEMON_SUCCESS = 0,
  LEMON_ERR_LAST_NOT_WRITTEN = -1,
  LEMON_ERR_PARAM = -2,
  LEMON_ERR_HEADER_NEXT = -3,
  LEMON_LAST_REC_WRITTEN = -4,
  LEMON_ERR_WRITE = -5,
  LEMON_EOR = -6,
  LEMON_EOF = -7,
  LEMON_ERR_READ = -8,
  LEMON_ERR_SEEK = -9,
  LEMON_ERR_MBME = -10,
  LEMON_ERR_CLOSE = -11
};
