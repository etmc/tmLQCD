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

paramsGaugeInfo GaugeInfo = { 0., 0, {0,0}, NULL, NULL};

int read_gauge_field(char * filename) {
  int status = 0;
  char *header_type = NULL;
  READER *reader = NULL;

  DML_Checksum checksum_read;
  DML_Checksum checksum_calc;
  int DML_read_flag = 0;
  int gauge_read_flag = 0;
  char *checksum_string = NULL;

  construct_reader(&reader, filename);
  GaugeInfo.gaugeRead = 0;
  while ((status = ReaderNextRecord(reader)) != LIME_EOF) {
    if (status != LIME_SUCCESS) {
      fprintf(stderr, "ReaderNextRecord returned status %d.\n", status);
      break;
    }
    header_type = ReaderType(reader);

    if(g_cart_id == 0 && g_debug_level > 1) {
      fprintf(stderr, "found header %s, will now read the message\n", header_type);
    }

    if (strcmp("ildg-binary-data", header_type) == 0) {
      read_binary_gauge_data(reader, &checksum_calc);
      gauge_read_flag = 1;
      GaugeInfo.gaugeRead = 1;
      GaugeInfo.checksum = checksum_calc;
    }
    else if (strcmp("scidac-checksum", header_type) == 0) {
      read_message(reader, &checksum_string);
      DML_read_flag = parse_checksum_xml(checksum_string, &checksum_read);
      free(checksum_string);
      checksum_string=(char*)NULL;
    }
    else if (strcmp("xlf-info", header_type) == 0) {
      read_message(reader, &GaugeInfo.xlfInfo);
    }
    else if (strcmp("ildg-data-lfn", header_type) == 0) {
      read_message(reader, &GaugeInfo.ildg_data_lfn);
    }
    close_reader_record(reader);
  }

  if (!gauge_read_flag)
    kill_with_error(reader->fp, g_cart_id, "Did not find gauge record. Aborting...\n");


  if (g_debug_level > 0 && g_cart_id == 0)
  {
    printf("# checksum for gaugefield %s\n", filename);
    printf("# calculated: %#x %#x.\n", checksum_calc.suma, checksum_calc.sumb);
    if (DML_read_flag)
      printf("# read:       %#x %#x.\n", checksum_read.suma, checksum_read.sumb);
    else
      printf("# Scidac checksum record not found or malformed.\n");
    fflush(stdout);
  }

  destruct_reader(reader);

  g_update_gauge_copy = 1;
  g_update_gauge_energy = 1;

  return 0;
}
