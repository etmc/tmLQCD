#define __DEFINE_RETURN_FIELD_ARRAY(DATATYPE, NAME)			 \
void return_ ## NAME ##_field_array(NAME ##_field_array_t * field_array) \
{                                                                        \
  for (unsigned int ctr = 0; ctr < field_array->length; ++ctr)           \
    return_ ## NAME ## _field(&field_array->field_array[ctr]);           \
  free(field_array->field_array);                                        \
}
