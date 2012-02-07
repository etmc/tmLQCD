#include "gauge.ih"

void initialize_gauge_buffers(unsigned int max)
{
  g_gauge_buffers.max = max;
  g_gauge_buffers.allocated = 0;
  g_gauge_buffers.free = 0;
  g_gauge_buffers.reserve = (su3_tuple**)calloc(max, sizeof(su3_tuple*));
}

