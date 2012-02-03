#include "gauge.ih"

void return_gauge_field(gauge_field_t *gauge_field)
{
  g_gauge_buffers.reserve[g_gauge_buffers.free] = gauge_field->field;
  ++g_gauge_buffers.free;
  gauge_field->field = NULL;
}

