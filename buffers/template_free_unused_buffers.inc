#define __DEFINE_FREE_UNUSED_BUFFERS(DATATYPE, NAME)									\
void free_unused_ ## NAME ## _buffers()                                                                                 \
{                                                                                                                       \
  for ( ; g_ ## NAME ## _buffers.free > 0; --g_ ## NAME ## _buffers.free, --g_ ## NAME ## _buffers.allocated)           \
  {                                                                                                                     \
    void* ptr = ((void**)g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.free - 1])[-1];                          \
    free(ptr);                                                                                                          \
  }                                                                                                                     \
}
