#include "gauge.ih"

gauge_field_array_t get_gauge_field_array(unsigned int length)
{
  gauge_field_array_t gauge_field_array;
  gauge_field_array.length = length;
  gauge_field_array.field_array = (gauge_field_t*)calloc(length, sizeof(gauge_field_t));

  if (g_gauge_buffers.free < length) /* Need to allocate more buffers */
    allocate_gauge_buffers(length - g_gauge_buffers.free);

  for (unsigned int ctr = 0; ctr < length; ++ctr)
  {
    --g_gauge_buffers.free;
    gauge_field_array.field_array[ctr].field = g_gauge_buffers.reserve[g_gauge_buffers.free];
    g_gauge_buffers.reserve[g_gauge_buffers.free] = NULL;
  }

  return gauge_field_array;
}

