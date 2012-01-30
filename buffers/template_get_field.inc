#define __DEFINE_GET_FIELD(DATATYPE, NAME)							\
NAME ## _field_t get_ ## NAME ## _field()                                                   \
{                                                                                               \
  NAME ## _field_t  NAME ## _field;                                                      \
                                                                                                \
  if (g_ ## NAME ## _buffers.stack < 0) /* Need to allocate a new buffer */                     \
    allocate_ ## NAME ## _buffers(1);                                                           \
                                                                                                \
  NAME ## _field.field = g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.stack];      \
  g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.stack] = NULL;                          \
  --g_ ## NAME ## _buffers.stack;                                                               \
                                                                                                \
  return  NAME ## _field;                                                                    \
}
