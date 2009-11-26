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

int read_lime_gauge_field(char * filename, DML_Checksum *scidac_checksum,
			 char **xlf_info, char **ildg_data_lfn) {
#ifdef HAVE_LIBLEMON
  read_lemon_gauge_field_parallel(filename, scidac_checksum, xlf_info, ildg_data_lfn);
#else
  FILE * ifs;
  int status;
  char * header_type;
  LimeReader * limereader;
  DML_Checksum checksum_read;
  DML_Checksum checksum_calc;
  int DML_read_flag = 0;
  char *checksum_string = NULL;
  int found_ildg_binary_data = 0;

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
#  ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#  endif
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
#  ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#  endif
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(strcmp("ildg-binary-data",header_type) == 0) {
      read_binary_gauge_data(limereader, &checksum_calc);
      found_ildg_binary_data = 1;
    }
    else if(strcmp("xlf-info", header_type) == 0) {
      read_message(limereader, xlf_info);
    }
    else if(strcmp("ildg-data-lfn", header_type) == 0) {
      read_message(limereader, ildg_data_lfn);
    }
    else if (strcmp("scidac-checksum", header_type) == 0) {
      read_message(limereader, &checksum_string);
      DML_read_flag = parse_checksum_xml(checksum_string, &checksum_read);
      if (DML_read_flag && scidac_checksum != (DML_Checksum*)NULL)
	*scidac_checksum = checksum_read;
      free(checksum_string); 
    }
  }
  if(found_ildg_binary_data == 0) {
    if(g_cart_id == 0) {
      fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
      fprintf(stderr, "trying old deprecated file format!\n");
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  if (g_debug_level > 0 && g_cart_id == 0) {
    printf("# checksum for gaugefield %s\n", filename);
    printf("# calculated: %#x %#x.\n", checksum_calc.suma, checksum_calc.sumb);
    if (DML_read_flag)
      printf("# read:       %#x %#x.\n", checksum_read.suma, checksum_read.sumb);
    else
      printf("# Scidac checksum record not found or malformed.\n");
    fflush(stdout);
  }


  limeDestroyReader(limereader);
  fclose(ifs);
  g_update_gauge_copy = 1;
  g_update_gauge_energy = 1;
#endif
  return(0);
}
