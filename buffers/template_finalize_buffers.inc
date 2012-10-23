#define __DEFINE_FINALIZE_BUFFERS(DATATYPE, NAME)							\
void finalize_ ## NAME ## _buffers()									\
{													\
  if (g_ ## NAME ## _buffers.free != g_ ## NAME ## _buffers.allocated)					\
  {													\
    fprintf(stderr, "[DEBUG] Not yet returned: %d - %d fields\n", g_ ## NAME ## _buffers.allocated, g_ ## NAME ## _buffers.free);	\
    fatal_error("Finalized g_" #NAME "_buffers with unreturned fields!", "finalize_" #NAME "_buffers");	\
  }													\
  free_unused_ ## NAME ## _buffers();									\
  free(g_ ## NAME ## _buffers.reserve);									\
  g_ ## NAME ## _buffers.max = 0;									\
}
