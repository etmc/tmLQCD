#define __DEFINE_ALLOCATE_BUFFER(DATATYPE, NAME) 					\
void allocate_ ## NAME ## _buffer() 					\
{								\
  if (g_ ## NAME ## _buffers.allocated == g_ ## NAME ## _buffers.max)			\
    fatal_error("Maximum number of allocated " #NAME " fields exceeded.", "allocate_" #NAME "_buffers");		\
								\
  void *p = aalloc(sizeof(DATATYPE) * (VOLUMEPLUSRAND + g_dbw2rand) + 1);			\
  if (p == NULL)							\
    fatal_error("Could not allocate the requested amount of memory.", "allocate_" #NAME "_buffers");		\
								\
  g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.free] = (DATATYPE*)p;	        		\
  ++g_ ## NAME ## _buffers.free;						\
  ++g_ ## NAME ## _buffers.allocated;						\
}
