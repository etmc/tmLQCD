/* $Id$ */
#ifndef _DML_H
#define _DML_H

/*
  Header file for the check sum from the QIO package.

*/

#include<stdlib.h>

typedef unsigned int   uint32_t;

/* from qio-2.2.0/include/dml.h  **/
typedef struct {
  uint32_t suma;
  uint32_t sumb;
} DML_Checksum;


typedef uint32_t DML_SiteRank;


/**
    Function prototypes
**/

void DML_checksum_init(DML_Checksum *checksum) ;
int DML_global_xor(uint32_t *x) ;

void DML_checksum_accum(DML_Checksum *checksum, DML_SiteRank rank,
                        char *buf, size_t size) ;

void DML_checksum_combine(DML_Checksum *checksum) ;


void DML_checksum_peq(DML_Checksum *total, DML_Checksum *checksum) ;

uint32_t DML_crc32(uint32_t crc, const unsigned char *buf, size_t len);

#endif



