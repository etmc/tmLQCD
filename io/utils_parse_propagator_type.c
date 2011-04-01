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

int parse_propagator_type(READER * reader) {
  char *prop_type_string = NULL;
  char *header_type = NULL;
  int prop_type = -1;
  int status = 0;
  int proptypefound = 0, sourcetypefound = 0;

  while ((status = ReaderNextRecord(reader)) != LIME_EOF) {
    if (status != LIME_SUCCESS) {
      fprintf(stderr, "ReaderNextRecord returned status %d.\n", status);
      break;
    }
    header_type = ReaderType(reader);
    if(g_cart_id == 0 && g_debug_level > 1) {
      fprintf(stdout, "found header %s, will now read the message\n", header_type);
      fflush(stdout);
    }
    if (strcmp("propagator-type", header_type) == 0) {
      read_message(reader, &prop_type_string);
      if(strcmp("DiracFermion_Sink", prop_type_string) == 0) 
        prop_type = 0;
      else if(strcmp("DiracFermion_Source_Sink_Pairs", prop_type_string) == 0)
        prop_type = 1;
      else if(strcmp("DiracFermion_ScalarSource_TwelveSink", prop_type_string) == 0)
        prop_type = 2;
      else if(strcmp("DiracFermion_ScalarSource_FourSink", prop_type_string) == 0)
        prop_type = 3;
      else if(strcmp("DiracFermion_Deflation_Field", prop_type_string) == 0)
        prop_type = 4;
      else {
        fprintf(stderr,"Unrecognized propagator-type, found type: %s.\n", prop_type_string);
        break;
      }
      proptypefound = 1;
      if(g_cart_id == 0 && g_debug_level > 0) {
        printf("# file is of type %s for proc %d\n", prop_type_string, g_cart_id);
      }
      free(prop_type_string);
      close_reader_record(reader);
      break;
    }
    if (strcmp("source-type", header_type) == 0) {
      read_message(reader, &prop_type_string);
      if(strcmp("DiracFermion_Source", prop_type_string) == 0)
        prop_type = 10;
      else if(strcmp("DiracFermion_ScalarSource", prop_type_string) == 0)
        prop_type = 11;
      else if(strcmp("DiracFermion_FourScalarSource", prop_type_string) == 0)
        prop_type = 12;
      else if(strcmp("DiracFermion_TwelveScalarSource", prop_type_string) == 0)
        prop_type = 13;
      else {
        fprintf(stderr,"Unrecognized source-type, found type: %s\n", prop_type_string);
        break;
      }
      sourcetypefound = 1;
      if(g_cart_id == 0 && g_debug_level > 0) {
        printf("# file is of type %s", prop_type_string);
      }
      free(prop_type_string);
      close_reader_record(reader);
      break;
    }
    if ((sourcetypefound || proptypefound) == 0) {
      fprintf(stderr, "Unable to find either source-type or propagator-type record.\nWARNING: Continuing in blind faith.\n");
    }
    close_reader_record(reader);
  }
  return(prop_type);
}
