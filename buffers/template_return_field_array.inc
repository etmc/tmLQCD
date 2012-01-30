#define __DEFINE_RETURN_FIELD_ARRAY(DATATYPE, NAME)									\
void return_ ## NAME ##_field_array(NAME ##_field_array_t * NAME ##_field_array)                                 \
{                                                                                                                       \
  for (unsigned int ctr = 0; ctr < NAME ##_field_array->length; ++ctr)                                              \
  {                                                                                                                     \
    ++g_ ## NAME ##_buffers.stack;                                                                                      \
    g_ ## NAME ##_buffers.reserve[g_ ## NAME ##_buffers.stack] = NAME ##_field_array->field_array[ctr].field;       \
    NAME ##_field_array->field_array[ctr].field = NULL;                                                             \
  }                                                                                                                     \
}
