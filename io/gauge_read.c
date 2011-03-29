/***********************************************************************
 *
 * Copyright (C) 2009-2011 Albert Deuzeman, Siebren Reker, Carsten Urbach
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
  int gauge_binary_status = 0;
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
      if (gauge_read_flag) { /* already read a gauge from this file */
        fprintf(stderr, "In gauge file %s, multiple LIME records with name: \"ildg-binary-data\" found.\n", filename);
        fprintf(stderr, "Unable to verify integrity of the gauge field data.\n");
        return(-1);
      }
      gauge_binary_status = read_binary_gauge_data(reader, &checksum_calc);
      gauge_read_flag = 1;
      GaugeInfo.gaugeRead = 1;
      GaugeInfo.checksum = checksum_calc;
    }
    else if (strcmp("scidac-checksum", header_type) == 0) {
      if(checksum_string == (char*)NULL) {
        read_message(reader, &checksum_string);
        DML_read_flag = parse_checksum_xml(checksum_string, &checksum_read);
        free(checksum_string);
      }
      else { /* checksum_string is not NULL, so a record was already found */
        fprintf(stderr, "In gauge file %s, multiple LIME records with name: \"scidac-checksum\" found.\n", filename);
        fprintf(stderr, "Unable to verify integrity of the gauge field data.\n");
        return(-1);
      }
    }
    else if (strcmp("xlf-info", header_type) == 0) {
      read_message(reader, &GaugeInfo.xlfInfo);
    }
    else if (strcmp("ildg-data-lfn", header_type) == 0) {
      read_message(reader, &GaugeInfo.ildg_data_lfn);
    }
    close_reader_record(reader);
  }

  if (!gauge_read_flag) {
    fprintf(stderr, "Unable to find LIME record with name: \"ildg-binary-data\", while looking in %s.\n", filename);
    fprintf(stderr, "Possible causes: %s is not a LIME file, or not of the appropriate type.\n", filename);
    fprintf(stderr, "No gauge field was read, unable to proceed.\n");
    return(-1);
  }

  if (!DML_read_flag) {
    fprintf(stderr, "LIME record with name: \"scidac-checksum\", in gauge file %s either missing or malformed.\n", filename);
    fprintf(stderr, "Unable to verify integrity of gauge field data.\n");
    return(-1);
  }

  if (g_debug_level >= 0 && g_cart_id == 0)
  {
    printf("# Scidac checksums for gaugefield %s\n", filename);
    printf("# Calculated            : A = %#x B = %#x.\n", checksum_calc.suma, checksum_calc.sumb);
    printf("# Read from LIME headers: A = %#x B = %#x.\n", checksum_read.suma, checksum_read.sumb);
    fflush(stdout);
    if (checksum_calc.suma != checksum_read.suma) {
      fprintf(stderr, "For gauge file %s, calculated and stored values for SciDAC checksum A do not match.\n", filename);
      return(-1);
    }
    if (checksum_calc.sumb != checksum_read.sumb) {
      fprintf(stderr, "For gauge file %s, calculated and stored values for SciDAC checksum B do not match.\n", filename);
      return(-1);
    }
  }

  destruct_reader(reader);

  g_update_gauge_copy = 1;
  g_update_gauge_energy = 1;

  return(0);
}
