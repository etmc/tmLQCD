#define __DECLARE_BUFFER_INTERFACE(DATATYPE, NAME)                                               \
typedef struct                                                                                   \
{                                                                                                \
  DATATYPE **reserve;                                                                            \
  unsigned int max;                                                                              \
  unsigned int allocated;                                                                        \
  unsigned int free;                                                                             \
} NAME ## _buffers_t;                                                                            \
                                                                                                 \
typedef struct                                                                                   \
{                                                                                                \
  DATATYPE *field;                                                                               \
} NAME ## _field_t;                                                                              \
                                                                                                 \
typedef struct                                                                                   \
{                                                                                                \
  NAME ## _field_t *field_array;                                                                 \
  unsigned int length;                                                                           \
} NAME ## _field_array_t;                                                                        \
                                                                                                 \
extern  NAME ## _buffers_t g_ ## NAME ## _buffers;                                               \
                                                                                                 \
void initialize_ ## NAME ## _buffers(unsigned int max);                                          \
void finalize_ ## NAME ## _buffers();                                                            \
                                                                                                 \
void allocate_ ## NAME ## _buffers(unsigned int count);                                          \
void free_unused_ ## NAME ## _buffers();                                                         \
                                                                                                 \
NAME ## _field_t get_ ## NAME ## _field();                                                       \
void return_ ## NAME ## _field(NAME ## _field_t *field);                                         \
                                                                                                 \
NAME ## _field_array_t get_ ## NAME ## _field_array(unsigned int length);                        \
void return_ ## NAME ## _field_array(NAME ## _field_array_t *field_array);                       \
                                                                                                 \
void swap_ ## NAME ## _field(NAME ## _field_t *left, NAME ## _field_t *right);                   \
void copy_ ## NAME ## _field(NAME ## _field_t *copy, NAME ## _field_t const original);
