#pragma once

#include <buffers/adjoint.h>
#include <buffers/complex.h>
#include <buffers/gauge.h>
#include <buffers/real.h>
#include <buffers/spinor.h>

void generic_exchange(void *field_in, int bytes_per_site);

void copy_adjoint_field(adjoint_field_t dest, adjoint_field_t orig);
void copy_complex_field(complex_field_t dest, complex_field_t orig);
void copy_gauge_field(gauge_field_t dest, gauge_field_t orig);
void copy_real_field(real_field_t dest, real_field_t orig);
void copy_spinor_field(spinor_field_t dest, spinor_field_t orig);

void zero_adjoint_field(adjoint_field_t target);
void zero_complex_field(complex_field_t target);
void zero_gauge_field(gauge_field_t target);
void zero_real_field(real_field_t target);
void zero_spinor_field(spinor_field_t target);

void zero_adjoint_field_array(adjoint_field_array_t target);
void zero_complex_field_array(complex_field_array_t target);
void zero_gauge_field_array(gauge_field_array_t target);
void zero_real_field_array(real_field_array_t target);
void zero_spinor_field_array(spinor_field_array_t target);

void exchange_adjoint_field(adjoint_field_t target);
void exchange_adjoint_field_array(adjoint_field_array_t target);

void exchange_complex_field(complex_field_t target);
void exchange_complex_field_array(complex_field_array_t target);

void exchange_gauge_field(gauge_field_t target);
void exchange_gauge_field_array(gauge_field_array_t target);

void exchange_real_field(real_field_t target);
void exchange_real_field_array(real_field_array_t target);

void exchange_spinor_field(spinor_field_t target);
void exchange_spinor_field_array(spinor_field_array_t target);

void convert_gauge_to_adjoint_field(adjoint_field_t dest, gauge_field_t orig);
void convert_adjoint_to_gauge_field(gauge_field_t dest, adjoint_field_t orig);

void convert_gauge_to_adjoint_field_array(adjoint_field_array_t dest, gauge_field_array_t orig);
void convert_adjoint_to_gauge_field_array(gauge_field_array_t dest, adjoint_field_array_t orig);

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
inline void zero_ ## NAME ## _field_array(NAME ## _field_array_t target)                        \
{                                                                                               \
  for (unsigned int idx = 0; idx < target.length; ++idx)                                        \
    zero_ ## NAME ##_field(target.field_array[idx]);	                                        \
}                                                                                               \
                                                                                                \
inline void exchange_ ## NAME ## _field(NAME ## _field_t target)                                \
{                                                                                               \
  generic_exchange((void*)target.field, sizeof(DATATYPE));                                      \
}                                                                                               \
                                                                                                \
inline void exchange_ ## NAME ## _field_array(NAME ## _field_array_t target)                    \
{                                                                                               \
  for (unsigned int idx = 0; idx < target.length; ++idx)                                        \
    exchange_ ## NAME ## _field(target.field_array[idx]);                                       \
}

__DEFINE_BUFFER_INLINES(su3adj_tuple, adjoint)
__DEFINE_BUFFER_INLINES(_Complex double, complex)
__DEFINE_BUFFER_INLINES(su3_tuple, gauge)
__DEFINE_BUFFER_INLINES(double, real)
__DEFINE_BUFFER_INLINES(spinor, spinor)

#undef __DEFINE_BUFFER_INLINES
