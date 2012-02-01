#define __DEFINE_RETURN_FIELD(DATATYPE, NAME)                                           \
void return_ ## NAME ## _field(NAME ## _field_t *NAME ## _field)                        \
{                                                                                       \
  ++g_ ## NAME ## _buffers.stack;                                                       \
  g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.stack] = NAME ## _field->field; \
  NAME ## _field->field = NULL;                                                         \
}
