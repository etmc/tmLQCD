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

#ifndef _PARAMS_H
#define _PARAMS_H

#include <io/dml.h>

typedef struct
{
  char   date[64];
  char   package_version[32];
  char   inverter[32];

  double epssq;
  double epsbar;
  double kappa;
  double mu;
  double mubar;
  double mu_inverted;
  double mu_lowest;

  double cgmms_mass;

  int    mms;
  int    iter;
  int    heavy;
  int    noflavours;
}
paramsInverterInfo;

typedef struct
{
  int    flavours;
  int    prec;
  int    lx;
  int    ly;
  int    lz;
  int    lt;
}
paramsPropagatorFormat;

typedef struct
{
  int    colours;
  int    flavours;
  int    prec;
  int    lx;
  int    ly;
  int    lz;
  int    lt;
  int    spins;
}
paramsSourceFormat;

typedef struct
{
  char   date[64];
  char   package_version[32];

  double beta;
  double c2_rec;
  double epsilonbar;
  double kappa;
  double mu;
  double mubar;
  double plaq;

  int    counter;

  long int time;
}
paramsXlfInfo;

typedef struct
{
  int    lx;
  int    ly;
  int    lz;
  int    lt;
  int    prec;
}
paramsIldgFormat;

typedef struct {
  double plaquetteEnergy;
  int gaugeRead;
  DML_Checksum checksum;
  char * xlfInfo;
  char * ildg_data_lfn;
} paramsGaugeInfo;

typedef struct {
  int splitted;
  int format;
  int precision;
  char * basename;
} paramsPropInfo;

typedef struct {
  /* later usage for the type of source */
  int type;
  /* splitted or not (really needed?) */
  int splitted;
  /* the IO format (needed?) */
  int format;
  /* the IO precision */
  int precision;
  /* the source location, where it applies*/
  int t, x, y, z;
  /* automatic TS detection */
  int automaticTS;
  /* sample, gauge no and index of source */
  int sample, nstore, ix;
  /* is this a 2 flavour source */
  int no_flavours;
  /* the base filename */
  char * basename;
} paramsSourceInfo;

/* defined in gauge_read.c */
extern paramsGaugeInfo GaugeInfo;
/* defined in spinor_read.c */
extern paramsPropInfo PropInfo;
extern paramsSourceInfo SourceInfo;

paramsIldgFormat       * construct_paramsIldgFormat(int const prec);
paramsPropagatorFormat * construct_paramsPropagatorFormat(int const prec, int const flavours);
paramsSourceFormat     * construct_paramsSourceFormat(int const prec, int const flavours, int const spins, int const sources);
paramsXlfInfo          * construct_paramsXlfInfo(double const plaq, int const counter);
paramsInverterInfo     * construct_paramsInverterInfo(double const epssq, const int iter, 
                                                      const int solver, const int noflavours);
#endif
