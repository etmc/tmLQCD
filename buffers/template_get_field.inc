#define __DEFINE_GET_FIELD(DATATYPE, NAME)					\
NAME ## _field_t get_ ## NAME ## _field()              				\
{										\
  NAME ## _field_t  field;							\
  if (g_ ## NAME ## _buffers.free == 0) /* Need to allocate a new buffer */	\
    allocate_ ## NAME ## _buffers(1);						\
  --g_ ## NAME ## _buffers.free;						\
  field = g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.free];		\
  g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.free] = NULL;		\
  return field;									\
}
