#pragma once

#include <buffers/adjoint.h>
#include <buffers/complex.h>
#include <buffers/gauge.h>
#include <buffers/real.h>
#include <buffers/spinor.h>

void generic_exchange(void *field_in, int bytes_per_site);

void zero_adjoint_field(adjoint_field_t target);
void zero_complex_field(complex_field_t target);
void zero_gauge_field(gauge_field_t target);
void zero_real_field(real_field_t target);
void zero_spinor_field(spinor_field_t target);

void exchange_adjoint_field(adjoint_field_t target);
void exchange_complex_field(complex_field_t target);
void exchange_gauge_field(gauge_field_t target);
void exchange_real_field(real_field_t target);
void exchange_spinor_field(spinor_field_t target);

void copy_adjoint_field(adjoint_field_t *left, adjoint_field_t const right);
void copy_complex_field(complex_field_t *left, complex_field_t const right);
void copy_gauge_field(gauge_field_t *left, copy_gauge_field const right);
void copy_real_field(real_field_t *left,real_field_t const right);
void copy_spinor_field(spinor_field_t *left, spinor_field_t const right);

void swap_adjoint_field(adjoint_field_t *left, adjoint_field_t *right);
void swap_complex_field(complex_field_t *left, complex_field_t *right);
void swap_gauge_field(gauge_field_t *left, swap_gauge_field *right);
void swap_real_field(real_field_t *left,real_field_t *right);
void swap_spinor_field(spinor_field_t *left, spinor_field_t *right);

void zero_adjoint_field(adjoint_field_t target);
void zero_complex_field(complex_field_t target);
void zero_gauge_field(gauge_field_t target);
void zero_real_field(real_field_t target);
void zero_spinor_field(spinor_field_t target);

/* Inline functions need to be declared inside the header -- hence the following nastiness here... */

#define __DEFINE_BUFFER_INLINES(DATATYPE, NAME)                                                 \
inline void copy_ ## NAME ## _field(NAME ## _field_t dest,  NAME ## _field_t orig)              \
{                                                                                               \
  memmove((void*)dest.field, (void*)orig.field, sizeof(DATATYPE) * VOLUMEPLUSRAND + 1);         \
}                                                                                               \
                                                                                                \
inline void zero_ ## NAME ##_field(NAME ## _field_t target)                                     \
{                                                                                               \
  memset((void*)target.field, 0.0, VOLUMEPLUSRAND * sizeof(DATATYPE));                          \
}                                                                                               \
                                                                                                \
inline void exchange_ ## NAME ## _field(NAME ## _field_t target)                                \
{                                                                                               \
  generic_exchange((void*)target.field, sizeof(DATATYPE));                                      \
}                                                                                               \
                                                                                                \
inline void swap_ ## NAME ## _field(NAME ## _field_t *left, NAME ## _field_t *right)            \
{                                                                                               \
   DATATYPE *tmp = left->field;                                                                 \
   left->field = right->field;                                                                  \
   right->field = tmp;                                                                          \
}                                                                                               \
                                                                                                \
inline void copy_ ## NAME ## _field(NAME ## _field_t *copy, NAME ## _field_t const original)    \
{                                                                                               \
  memmove(copy->field, original.field, VOLUMEPLUSRAND * sizeof(DATATYPE));                      \
}

__DEFINE_BUFFER_INLINES(su3adj_tuple, adjoint)
__DEFINE_BUFFER_INLINES(_Complex double, complex)
__DEFINE_BUFFER_INLINES(su3_tuple, gauge)
__DEFINE_BUFFER_INLINES(double, real)
__DEFINE_BUFFER_INLINES(spinor, spinor)

#undef __DEFINE_BUFFER_INLINES
