#include "global.h"
#include "deflation_block.h"

#define CALLOC_ERROR_CRASH {printf ("calloc errno : %d\n", errno); errno = 0; return 1;}

extern deflation_block *g_deflation_blocks;

int init_deflation_blocks()
{
  int i, j;

  g_deflation_blocks = calloc(2, sizeof(deflation_block));
  for (i = 0; i < 2; ++i)
  {
    if ((void*)(g_deflation_blocks[i].little_basis = calloc(LITTLE_BASIS_SIZE, sizeof(su3spinor *))) == NULL)
      CALLOC_ERROR_CRASH;
    
    for (j = 0; j < LITTLE_BASIS_SIZE; ++j)
      if ((void*)(g_deflation_blocks[i].little_basis[j] = calloc(VOLUME, sizeof(su3spinor))) == NULL)
        CALLOC_ERROR_CRASH;
    
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges = calloc(8, sizeof(su3spinor *))) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[0] = calloc(VOLUME / (2 * LX), sizeof(su3spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[1] = calloc(VOLUME / (2 * LX), sizeof(su3spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[2] = calloc(VOLUME / (2 * LY), sizeof(su3spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[3] = calloc(VOLUME / (2 * LY), sizeof(su3spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[4] = calloc(VOLUME / LZ, sizeof(su3spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[5] = calloc(VOLUME / LZ, sizeof(su3spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[6] = calloc(VOLUME / (2 * T), sizeof(su3sp inor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[7] = calloc(VOLUME / (2 * T), sizeof(su3spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    
    if ((void*)(g_deflation_blocks[i].little_dirac_operator = calloc(9 * LITTLE_BASIS_SIZE * LITTLE_BASIS_SIZE, complex)) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_deflation_blocks[i].local_little_field = calloc(LITTLE_BASIS_SIZE, complex)) == NULL)
      CALLOC_ERROR_CRASH;
  }
  return 0;
}

int free_deflation_blocks()
{
  int i, j;

  for (i = 0; i < 2; ++i)
  {   
    for (j = 0; j < LITTLE_BASIS_SIZE; ++j)
      free(g_deflation_blocks[i].little_basis[j]);
    free(g_deflation_blocks[i].little_basis);
    
    for (j = 0; j < 8; ++j)
      free(g_deflation_blocks[i].little_neighbour_edges[j]);
    free(g_deflation_blocks[i].little_neighbour_edges);
    
    free(g_deflation_blocks[i].little_dirac_operator);
    free(g_deflation_blocks[i].local_little_field);
  }
  
  free(g_deflation_blocks);
  return 0;
}

int add_basis_field(int const index, su3spinor const *field)
{
  /* This should in principle include a transposition step if we want to use LAPACK for all
     the low and dirty linear algebra. Specifically, we would like to use its SVD routines here. */
  int ctr_t;
  int contig_block = VOLUME / (2 * T);
  for (ctr_t = 0; ctr_t < T; ++ctr_t)
  {
    memcpy(g_deflation_block[0].little_basis[index], field + (2 * T) * contig_block, contig_block * sizeof(su3spinor));
    memcpy(g_deflation_block[1].little_basis[index], field + (2 * T + 1) * contig_block, contig_block * sizeof(su3spinor));
  }
  return 0;
}
