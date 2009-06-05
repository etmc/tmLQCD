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
#define PACKAGE_VERSION "5.0.1"

void write_xlf_info_parallel(LemonWriter * lemonwriter, const double plaq, const int counter)
{
  char *message;
  uint64_t bytes;
  struct timeval t1;

  message = (char*)malloc(512);
  if (message == (char*)NULL )
  {
    kill_with_error(lemonwriter->fh, lemonwriter->my_rank,
                    "Memory allocation error in write_xlf_info_parallel. Aborting\n");
  }
  gettimeofday(&t1,NULL);
  if(g_kappa > 0. || g_kappa < 0.)
  {
    sprintf(message,"\n plaquette = %e\n trajectory nr = %d\n beta = %f, kappa = %f, mu = %f, c2_rec = %f\n time = %ld\n hmcversion = %s\n mubar = %f\n epsilonbar = %f\n date = %s",
        plaq, counter, g_beta, g_kappa, g_mu/2./g_kappa, g_rgi_C1,t1.tv_sec, PACKAGE_VERSION,
        g_mubar/2./g_kappa, g_epsbar/2./g_kappa, ctime(&t1.tv_sec));
  }
  else
  {
    sprintf(message,"\n plaquette = %e\n trajectory nr = %d\n beta = %f, kappa = %f, 2*kappa*mu = %f, c2_rec = %f\n date = %s",
        plaq, counter, g_beta, g_kappa, g_mu, g_rgi_C1, ctime(&t1.tv_sec));
  }
  bytes = strlen( message );

  write_header_parallel(lemonwriter, 1, 1, "xlf-info", bytes);
  write_message_parallel(lemonwriter, message, bytes);

  lemonWriterCloseRecord(lemonwriter);
  free(message);
}
