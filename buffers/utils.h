#ifndef _BUFFERS_UTILS_H
#define _BUFFERS_UTILS_H

#include <string.h>

#include <buffers/adjoint.h>
#include <buffers/gauge.h>
#include <buffers/spinor.h>

extern int VOLUMEPLUSRAND;

void generic_exchange(void *field_in, int bytes_per_site);

static inline void zero_adjoint_field(adjoint_field_t *target);
static inline void zero_gauge_field(gauge_field_t *target);

static inline void exchange_adjoint_field(adjoint_field_t *target);
static inline void exchange_gauge_field(gauge_field_t *target);

static inline void copy_adjoint_field(adjoint_field_t *left, adjoint_field_t const right);
static inline void copy_gauge_field(gauge_field_t *left, gauge_field_t const right);

static inline void swap_adjoint_field(adjoint_field_t *left, adjoint_field_t *right);
static inline void swap_gauge_field(gauge_field_t *left, gauge_field_t *right);

void adjoint_to_gauge(gauge_field_t *out, adjoint_field_t const in);
void gauge_to_adjoint(adjoint_field_t *out, gauge_field_t const in);

/* Inline functions need to be declared inside the header -- hence the following nastiness here... */

#define __DEFINE_BUFFER_INLINES(DATATYPE, NAME)                                                 \
static inline void zero_ ## NAME ##_field(NAME ## _field_t *target)                                     \
{                                                                                               \
  memset(*target, 0.0, VOLUMEPLUSRAND * sizeof(DATATYPE));                                \
}                                                                                               \
                                                                                                \
static inline void exchange_ ## NAME ## _field(NAME ## _field_t *target)                                \
{                                                                                               \
  generic_exchange(*target, sizeof(DATATYPE));                                            \
}                                                                                               \
                                                                                                \
static inline void swap_ ## NAME ## _field(NAME ## _field_t *left, NAME ## _field_t *right)            \
{                                                                                               \
   DATATYPE *tmp = *left;                                                                       \
   *left = *right;					 \
   *right = tmp;                                                                                \
}                                                                                               \
                                                                                                \
static inline void copy_ ## NAME ## _field(NAME ## _field_t *copy, NAME ## _field_t const original)    \
{                                                                                               \
  memmove(*copy, original, VOLUMEPLUSRAND * sizeof(DATATYPE));                                  \
}

__DEFINE_BUFFER_INLINES(su3adj_tuple, adjoint)
__DEFINE_BUFFER_INLINES(su3_tuple, gauge)
// __DEFINE_BUFFER_INLINES(spinor, spinor)

#undef __DEFINE_BUFFER_INLINES

#endif
