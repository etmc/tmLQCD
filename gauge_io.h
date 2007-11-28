/* $Id$ */

#ifndef _GAUGE_IO_H
#define _GAUGE_IO_H

#include"dml.h"

/* int write_binary_gauge_data(LimeWriter * limewriter, */
/* 			    const int prec, DML_Checksum * ans); */
/* int read_binary_gauge_data(LimeReader * limereader,  */
/* 			   const double prec, DML_Checksum * ans); */

int write_lime_gauge_field(char * filename, const double plaq, const int counter, const int prec);
int read_lime_gauge_field(char * filename);

#endif
