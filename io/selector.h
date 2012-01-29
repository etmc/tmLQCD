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

#ifndef _IO_SELECTOR_H
#define _IO_SELECTOR_H

#include <lime.h>
#ifdef HAVE_LIBLEMON
# include <lemon.h>
#endif /* HAVE_LIBLEMON */

#ifdef HAVE_LIBLEMON
#  define LIME_FILE MPI_File
#  define WRITER LemonWriter
#  define READER LemonReader
#  define RECORD_HEADER LemonRecordHeader
#  define CreateReader lemonCreateReader
#  define CreateHeader lemonCreateHeader
#  define ReaderBytes lemonReaderBytes
#  define ReaderNextRecord lemonReaderNextRecord
#  define ReaderType lemonReaderType
#  define ReaderCloseRecord lemonReaderCloseRecord
#  define ReaderReadData lemonReaderReadData
#  define WriteRecordHeader lemonWriteRecordHeader
#  define WriteRecordData lemonWriteRecordData
#  define WriterCloseRecord lemonWriterCloseRecord
#  define DestroyReader lemonDestroyReader
#  define DestroyHeader lemonDestroyHeader
#else /* HAVE_LIBLEMON */
#  define LIME_FILE FILE
#  define WRITER LimeWriter
#  define READER LimeReader
#  define RECORD_HEADER LimeRecordHeader
#  define CreateReader limeCreateReader
#  define CreateHeader limeCreateHeader
#  define ReaderBytes limeReaderBytes
#  define ReaderNextRecord limeReaderNextRecord
#  define ReaderType limeReaderType
#  define ReaderCloseRecord limeReaderCloseRecord
#  define ReaderReadData limeReaderReadData
#  define WriteRecordData limeWriteRecordData
#  define WriteRecordHeader limeWriteRecordHeader
#  define WriterCloseRecord limeWriterCloseRecord
#  define DestroyReader limeDestroyReader
#  define DestroyHeader limeDestroyHeader
#endif

#endif
