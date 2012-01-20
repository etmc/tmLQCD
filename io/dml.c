/*
  A subset of the dml library for checksums
  taken from QIO and adapted for tmLQCD by Carsten Urbach
*/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include"dml.h"


/*------------------------------------------------------------------*/
/* Checksum "class" */
/* We do a crc32 sum on the site data -- then do two lexicographic-
   rank-based bit rotations and XORs on the resulting crc32
   checksum */

/* Initialize checksums */
void DML_checksum_init(DML_Checksum *checksum){
  checksum->suma = 0;
  checksum->sumb = 0;
}


#ifdef MPI
int DML_global_xor(uint32_t *x) {
  unsigned long work = (unsigned long)*x;
  unsigned long dest;
  int status;

  status = MPI_Allreduce((void *)&work, (void *)&dest, 1,
                         MPI_UNSIGNED_LONG, MPI_BXOR, MPI_COMM_WORLD);

  if (status == MPI_SUCCESS) {
    *x = (uint32_t)dest;
  }
  return(status);
}
#else
int DML_global_xor(uint32_t *x){return(0);}
#endif


/* Accumulate checksums */
void DML_checksum_accum(DML_Checksum *checksum, DML_SiteRank rank,
                        char *buf, size_t size){

  DML_SiteRank rank29 = rank;
  DML_SiteRank rank31 = rank;
  uint32_t work = DML_crc32(0, (unsigned char*)buf, size);

  rank29 %= 29; rank31 %= 31;

  checksum->suma ^= work<<rank29 | work>>(32-rank29);
  checksum->sumb ^= work<<rank31 | work>>(32-rank31);
}


/* Combine checksums over all nodes */
void DML_checksum_combine(DML_Checksum *checksum){
  DML_global_xor(&checksum->suma);
  DML_global_xor(&checksum->sumb);
}

/* Add single checksum set to the total */
void DML_checksum_peq(DML_Checksum *total, DML_Checksum *checksum){
  total->suma ^= checksum->suma;
  total->sumb ^= checksum->sumb;
}

