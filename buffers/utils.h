#pragma once

#include <buffers/adjoint.h>
#include <buffers/gauge.h>
#include <buffers/spinor.h>

void generic_exchange(void *field_in, int bytes_per_site);

void zero_adjoint_field(adjoint_field_t target);
void zero_gauge_field(gauge_field_t target);
void zero_spinor_field(spinor_field_t target);

void exchange_adjoint_field(adjoint_field_t target);
void exchange_gauge_field(gauge_field_t target);
void exchange_spinor_field(spinor_field_t target);

void copy_adjoint_field(adjoint_field_t *left, adjoint_field_t const right);
void copy_gauge_field(gauge_field_t *left, gauge_field_t const right);
void copy_spinor_field(spinor_field_t *left, spinor_field_t const right);

void swap_adjoint_field(adjoint_field_t *left, adjoint_field_t *right);
void swap_gauge_field(gauge_field_t *left, gauge_field_t *right);
void swap_spinor_field(spinor_field_t *left, spinor_field_t *right);

void zero_adjoint_field(adjoint_field_t target);
void zero_gauge_field(gauge_field_t target);
void zero_spinor_field(spinor_field_t target);

void adjoint_to_gauge(gauge_field_t out, adjoint_field_t const in);
void gauge_to_adjoint(adjoint_field_t out, gauge_field_t const in);

/* Inline functions need to be declared inside the header -- hence the following nastiness here... */

#define __DEFINE_BUFFER_INLINES(DATATYPE, NAME)                                                 \
inline void zero_ ## NAME ##_field(NAME ## _field_t target)                                     \
{                                                                                               \
  memset((void*)target, 0.0, VOLUMEPLUSRAND * sizeof(DATATYPE));                          \
}                                                                                               \
                                                                                                \
inline void exchange_ ## NAME ## _field(NAME ## _field_t target)                                \
{                                                                                               \
  generic_exchange((void*)target, sizeof(DATATYPE));                                      \
}                                                                                               \
                                                                                                \
inline void swap_ ## NAME ## _field(NAME ## _field_t *left, NAME ## _field_t *right)            \
{                                                                                               \
   DATATYPE *tmp = *left;                                                                 \
   *left = *right;									\
   *right = tmp;                                                                          \
}                                                                                               \
                                                                                                \
inline void copy_ ## NAME ## _field(NAME ## _field_t *copy, NAME ## _field_t const original)    \
{                                                                                               \
  memmove(*copy, original, VOLUMEPLUSRAND * sizeof(DATATYPE));                      \
}

__DEFINE_BUFFER_INLINES(su3adj_tuple, adjoint)
__DEFINE_BUFFER_INLINES(su3_tuple, gauge)
__DEFINE_BUFFER_INLINES(spinor, spinor)

#undef __DEFINE_BUFFER_INLINES
