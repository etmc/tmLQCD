#define __DEFINE_INITIALIZE_BUFFERS(DATATYPE, NAME)			\
void initialize_ ## NAME ## _buffers(unsigned int max)			\
{						\
  g_ ## NAME ## _buffers.max = max;				\
  g_ ## NAME ## _buffers.allocated = 0;				\
  g_ ## NAME ## _buffers.free = 0;				\
  g_ ## NAME ## _buffers.reserve = (DATATYPE**)calloc(max, sizeof(DATATYPE*));	\
}
