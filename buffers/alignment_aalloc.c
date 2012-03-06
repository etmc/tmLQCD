#include "alignment.ih"

void *aalloc(size_t const bytes)
{
  void *raw = malloc(sizeof(void*) + ALIGN_BASE + bytes);

  if (raw == NULL)
    return NULL;

  size_t ptr = (size_t)raw + sizeof(void*);
  ptr = ((ptr + ALIGN_BASE) & ~ALIGN_BASE);
  ((void**)ptr)[-1] = raw;
  return (void*)ptr;
}
