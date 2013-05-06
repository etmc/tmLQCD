#ifndef _BUFFERS_ALIGNMENT_H
#define _BUFFERS_ALIGNMENT_H

#include <stdlib.h>
#include <string.h>

typedef struct
{
  int    id;
  size_t bytes;
  char   note[256];
} meta_data_t;

void *aalloc(size_t const bytes);
void *aalloc_annoted(size_t const bytes, char const *note);
static inline void afree(void *ptr);

static inline void ameta_set_note(void *ptr, char const *note);
static inline size_t ameta_bytes(void *ptr);
static inline int ameta_id(void *ptr);
static inline char const *ameta_note(void *ptr);

static inline void afree(void *ptr)
{
  free(((void**)ptr)[-1]);
}

static inline void ameta_set_note(void *ptr, char const *note)
{
  strncpy(((meta_data_t*)((void**)ptr)[-1])->note, note, 255);
  ((meta_data_t*)((void**)ptr)[-1])->note[255] = '\0';
}

static inline size_t ameta_bytes(void *ptr)
{
  return ((meta_data_t*)((void**)ptr)[-1])->bytes;
}

static inline int ameta_id(void *ptr)
{
  return ((meta_data_t*)((void**)ptr)[-1])->id;
}

static inline char const *ameta_note(void *ptr)
{
  return ((meta_data_t*)((void**)ptr)[-1])->note;
}

#endif

