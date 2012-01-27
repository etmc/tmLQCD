#pragma once

#include <buffers/gauge.h>

void copy_gauge_field(gauge_field_t dest, gauge_field_t orig)
void generic_exchange(void *field_in, int bytes_per_site);
void exchange_gauge_field(gauge_field_t target);

inline void copy_gauge_field(gauge_field_t dest, gauge_field_t orig)
{
  memmove((void*)dest.field, (void*)orig.field, sizeof(su3_tuple) * VOLUMEPLUSRAND + 1);
}

inline void exchange_gauge_field(gauge_field_t target)
{
  generic_exchange((void*)target.field, sizeof(su3_tuple));
}
