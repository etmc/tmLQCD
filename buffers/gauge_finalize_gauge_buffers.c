#include "gauge.ih"

void finalize_gauge_buffers()
{
  if (g_gauge_buffers.free != g_gauge_buffers.allocated)
    fatal_error("Finalized g_gauge_buffers with unreturned fields!", "finalize_gauge_buffers");

  free_unused_gauge_buffers();
  free(g_gauge_buffers.reserve);
  g_gauge_buffers.max = 0;
}
