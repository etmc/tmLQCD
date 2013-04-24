#ifndef _BUFFERS_ALIGNMENT_H
#define _BUFFERS_ALIGNMENT_H

#include <stdlib.h>

typedef struct
{
  size_t bytes;
  int id;
} meta_data_t;

void *aalloc(size_t const bytes);
static inline void afree(void *ptr);
static inline size_t ameta_bytes(void *ptr);
static inline int ameta_id(void *ptr);

static inline void afree(void *ptr)
{
  free(((void**)ptr)[-1]);
}

static inline size_t ameta_bytes(void *ptr)
{
  return ((meta_data_t*)((void**)ptr)[-1])->bytes;
}

static inline int ameta_id(void *ptr)
{
  return ((meta_data_t*)((void**)ptr)[-1])->id;
}

#endif

