#define __DEFINE_ALLOCATE_BUFFERS(DATATYPE, NAME) 								\
void allocate_ ## NAME ## _buffers(unsigned int count) 								\
{														\
  if ((g_ ## NAME ## _buffers.allocated + count) > g_ ## NAME ## _buffers.max)					\
    fatal_error("Maximum number of allocated " #NAME " fields exceeded.", "allocate_" #NAME "_buffers");	\
														\
  for (unsigned int ctr = 0; ctr < count; ++ctr)								\
  {														\
    void *p = aalloc(sizeof(DATATYPE) * VOLUMEPLUSRAND + 1);			                                \
    if (p == NULL)												\
      fatal_error("Could not allocate the requested amount of memory.", "allocate_" #NAME "_buffers");		\
														\
    g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.free] = (DATATYPE*)p;	        			\
    ++g_ ## NAME ## _buffers.free;                                                                             \
    ++g_ ## NAME ## _buffers.allocated;										\
  }														\
}
