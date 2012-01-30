#define __DEFINE_FREE_UNUSED_BUFFERS(DATATYPE, NAME)									\
void free_unused_ ## NAME ## _buffers()                                                                                 \
{                                                                                                                       \
  for ( ; g_ ## NAME ## _buffers.stack >= 0; --g_ ## NAME ## _buffers.stack, --g_ ## NAME ## _buffers.allocated)        \
  {                                                                                                                     \
    void* ptr = ((void**)g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.stack])[-1];                             \
    free(ptr);                                                                                                          \
  }                                                                                                                     \
}