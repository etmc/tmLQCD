#include "gauge.ih"

void free_unused_gauge_buffers()
{
  for ( ; g_gauge_buffers.free > 0; --g_gauge_buffers.free, --g_gauge_buffers.allocated)
  {
    void* ptr = ((void**)g_gauge_buffers.reserve[g_gauge_buffers.free - 1])[-1];
    free(ptr);
  }
}
