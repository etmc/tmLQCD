/* $Id$ */
#ifndef _IO_H
#define _IO_H

int read_lime_gauge_field_doubleprec(double * config, char * filename,
				     const int T, const int LX, const int LY, const int LZ);

int read_lime_gauge_field_singleprec(float * config, char * filename,
				     const int T, const int LX, const int LY, const int LZ);

#endif
