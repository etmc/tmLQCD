#define __DEFINE_GET_FIELD_ANNOTATED(DATATYPE, NAME)                            \
NAME ## _field_t get_ ## NAME ## _field_annotated(char const *note)       	\
{					\
  NAME ## _field_t  field;				\
  if (g_ ## NAME ## _buffers.free == 0) /* Need to allocate a new buffer */	\
    allocate_ ## NAME ## _buffer(note);			\
  --g_ ## NAME ## _buffers.free;			\
  field = g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.free];	\
  g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.free] = NULL;	\
  ameta_set_note(field, note);                                                  \
  return field;					\
}
