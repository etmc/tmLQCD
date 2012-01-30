#define __DEFINE_FINALIZE_BUFFERS(DATATYPE, NAME)								\
void finalize_ ## NAME ## _buffers()                                                                            \
{                                                                                                               \
  if (g_ ## NAME ## _buffers.stack != g_ ## NAME ## _buffers.allocated - 1)                                     \
    fatal_error("Finalized g_ ## NAME ## _buffers with unreturned fields!", "finalize_ ## NAME ## _buffers");   \
  for (unsigned int ctr = 0; ctr < (g_ ## NAME ## _buffers.stack); ++ctr)                                       \
  {                                                                                                             \
    void* unaligned_ptr = ((void**)g_ ## NAME ## _buffers.reserve[ctr])[-1];                                    \
    free(unaligned_ptr);                                                                                        \
  }                                                                                                             \
  free(g_ ## NAME ## _buffers.reserve);                                                                         \
}
