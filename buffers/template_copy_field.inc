#define __DEFINE_COPY_FIELD(DATATYPE, NAME)                                           \
void swap_ ## NAME ## _field(NAME ## _field_t *copy, NAME ## _field_t const original) \
{                                                                                     \
  memmove(copy->field, original.field, VOLUMEPLUSRAND * sizeof(DATATYPE));            \
}
