#define __DEFINE_GET_FIELD_ARRAY(DATATYPE, NAME)					   \
NAME ## _field_array_t get_ ## NAME ## _field_array(unsigned int length)                   \
{                                                                                          \
   NAME ## _field_array_t field_array;                                                     \
   field_array.length = length;                                                   \
   field_array.field_array = (NAME ## _field_t*)calloc(length, sizeof(NAME ## _field_t));  \
                                                                                           \
  for (unsigned int ctr = 0; ctr < length; ++ctr)                                          \
    field_array.field_array[ctr] = get_ ## NAME ## _field();                \
                                                                                           \
  return field_array;                                                            \
}
