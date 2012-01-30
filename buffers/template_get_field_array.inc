#define __DEFINE_GET_FIELD_ARRAY(DATATYPE, NAME)									\
NAME ## _field_array_t get_ ## NAME ## _field_array(unsigned int length)                                             \
{                                                                                                                       \
   NAME ## _field_array_t  NAME ## _field_array;                                                                  	\
   NAME ## _field_array.length = length;                                                                             	\
   NAME ## _field_array.field_array = (NAME ## _field_t*)calloc(length, sizeof(NAME ## _field_t));           		\
                                                                                                                        \
  if (g_ ## NAME ## _buffers.stack < (length - 1)) /* Need to allocate more buffers */                                  \
    allocate_ ## NAME ## _buffers(length - g_ ## NAME ## _buffers.stack - 1);                                           \
                                                                                                                        \
  for (unsigned int ctr = 0; ctr < length; ++ctr)                                                                       \
  {                                                                                                                     \
    NAME ## _field_array.field_array[ctr].field = g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.stack];     	\
    g_ ## NAME ## _buffers.reserve[g_ ## NAME ## _buffers.stack] = NULL;                                                \
    --g_ ## NAME ## _buffers.stack;                                                                                     \
  }                                                                                                                     \
                                                                                                                        \
  return  NAME ## _field_array;                                                                                      	\
}
