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

#include "lemon_header.h"

typedef struct
{
  int is_last;
  int header_nextP;
  MPI_File *fh;
  LemonRecordHeader *curr_header;

  MPI_Offset bytes_total;
  size_t bytes_pad;

  MPI_Comm cartesian;
  int      my_rank;

  MPI_Offset off;
  MPI_Offset pos;
} LemonReader;

/* Reader manipulators */
LemonReader* lemonCreateReader(MPI_File *fp, MPI_Comm cartesian);
void lemonDestroyReader(LemonReader *reader);
int lemonSetReaderPointer(LemonReader *reader, MPI_Offset offset);
MPI_Offset lemonGetReaderPointer(LemonReader *reader);

int lemonReaderNextRecord(LemonReader *reader);
int lemonReaderMBFlag(LemonReader *reader);
int lemonReaderMEFlag(LemonReader *reader);
char *lemonReaderType(LemonReader *reader);
uint64_t lemonReaderBytes(LemonReader *reader);
size_t lemonReaderPadBytes(LemonReader *reader);

int lemonReaderReadData(void *dest, uint64_t *nbytes, LemonReader *reader);
int lemonReaderCloseRecord(LemonReader *reader);
int lemonReaderSeek(LemonReader *reader, MPI_Offset offset, int whence);
int lemonReaderSetState(LemonReader *rdest, LemonReader *rsrc);

/* WAS THE FOLLOWING EVER IMPLEMENTED? */
int lemonEOM(LemonReader *reader);

/* Additions for LEMON follow */
int lemonReadLatticeParallel(LemonReader *reader, void *data,
                            MPI_Offset siteSize, int *latticeDims);
