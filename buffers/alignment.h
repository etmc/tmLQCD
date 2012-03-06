#pragma once

#include <stdlib.h>

void *aalloc(size_t const bytes);
void afree(void *ptr);

inline void afree(void *ptr)
{
  free(((void**)ptr)[-1]);
}