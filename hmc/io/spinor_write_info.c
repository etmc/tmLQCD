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

#include "spinor.ih"

void write_spinor_info(WRITER * writer, const int write_prop_format_flag,
                       paramsInverterInfo * InverterInfo)
{
  write_propagator_type(writer, write_prop_format_flag);
  if(GaugeInfo.xlfInfo != NULL) {
    write_header(writer, 1, 1, "xlf-info", strlen(GaugeInfo.xlfInfo));
    write_message(writer, GaugeInfo.xlfInfo, strlen(GaugeInfo.xlfInfo));
    close_writer_record(writer);
  }
  write_inverter_info(writer, InverterInfo);
  if (GaugeInfo.ildg_data_lfn != NULL)
  {
    write_header(writer, 1, 1, "gauge-ildg-data-lfn-copy", strlen(GaugeInfo.ildg_data_lfn));
    write_message(writer, GaugeInfo.ildg_data_lfn, strlen(GaugeInfo.ildg_data_lfn));
    close_writer_record(writer);
  }
  write_checksum(writer, &GaugeInfo.checksum, "gauge-scidac-checksum-copy");
}
