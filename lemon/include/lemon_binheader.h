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

#include <lemon.h>

static union
{
  uint64_t      int64[HDR_SIZE_BYTES / 8];
  uint32_t      int32[HDR_SIZE_BYTES / 4];
  uint16_t      int16[HDR_SIZE_BYTES / 2];
  unsigned char uchr [HDR_SIZE_BYTES    ];
} lemon_header;

static uint32_t      *lemon_hdr_magic_no = &lemon_header.int32[ 0];
static uint16_t      *lemon_hdr_version  = &lemon_header.int16[ 2];
static uint64_t      *lemon_hdr_data_len = &lemon_header.int64[ 1];
static unsigned char *lemon_hdr_mbme     = &lemon_header.uchr [ 6];
static unsigned char *lemon_hdr_rec_type = &lemon_header.uchr [16];
