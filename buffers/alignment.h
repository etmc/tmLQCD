#ifndef _BUFFERS_ALIGNMENT_H
#define _BUFFERS_ALIGNMENT_H

#include <stdlib.h>

void *aalloc(size_t const bytes);
static inline void afree(void *ptr);

static inline void afree(void *ptr)
{
  free(((void**)ptr)[-1]);
}

#endif

