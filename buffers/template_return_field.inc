#define __DEFINE_RETURN_FIELD(DATATYPE, NAME)                                           \
void return_ ## NAME ## _field(NAME ## _field_t *field)                        \
{                                                                                       \
  g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.free] = field->field;  \
  ++g_ ## NAME ## _buffers.free;                                                        \
  field->field = NULL;                                                         \
}
