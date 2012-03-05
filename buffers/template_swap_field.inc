#define __DEFINE_SWAP_FIELD(DATATYPE, NAME)					\
void swap_ ## NAME ## _field(NAME ## _field_t *left, NAME ## _field_t *right)   \
{                                                                               \
   DATATYPE *tmp = left->field;                                                 \
   left->field = right->field;                                                  \
   right->field = tmp;                                                          \
}
