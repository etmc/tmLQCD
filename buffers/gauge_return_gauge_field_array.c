#include "gauge.ih"

void return_gauge_field_array(gauge_field_array_t *gauge_field_array)
{
  for (unsigned int ctr = 0; ctr < gauge_field_array->length; ++ctr)
  {
    ++g_gauge_buffers.stack;
    g_gauge_buffers.reserve[g_gauge_buffers.stack] = gauge_field_array->field_array[ctr];
    gauge_field_array->field_array[ctr] = NULL;
  }
}
