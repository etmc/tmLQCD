#define __DEFINE_ALLOCATE_BUFFERS(DATATYPE, NAME) 								\
void allocate_ ## NAME ## _buffers(unsigned int count) 								\
{														\
  if ((g_ ## NAME ## _buffers.allocated + count) > g_ ## NAME ## _buffers.max)					\
    fatal_error("Maximum number of allocated  ## NAME ##  fields exceeded.", "allocate_ ## NAME ## _buffers");	\
														\
  for (unsigned int ctr = 0; ctr < count; ++ctr)								\
  {														\
    void *raw = malloc(sizeof(void*) + ALIGN_BASE + sizeof(DATATYPE) * VOLUMEPLUSRAND + 1);			\
    if (raw == NULL)												\
      fatal_error("Could not allocate the requested amount of memory.", "allocate_ ## NAME ## _buffers");	\
    size_t p = (size_t)raw + sizeof(void*);									\
    p = ((p + ALIGN_BASE) & ~ALIGN_BASE);									\
    ((void**)p)[-1] = raw;											\
														\
    g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.free] = (DATATYPE*)p;	        			\
    ++g_ ## NAME ## _buffers.free;                                                                              \
    ++g_ ## NAME ## _buffers.allocated;										\
  }														\
}
