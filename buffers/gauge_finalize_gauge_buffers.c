#include "gauge.ih"

void finalize_gauge_buffers()
{
  if (g_gauge_buffers.stack != g_gauge_buffers.allocated - 1)
    fprintf(stderr, "WARNING: g_gauge_buffers finalized with unreturned fields!"); /* Make error? */
  for (unsigned int ctr = 0; ctr < (g_gauge_buffers.stack); ++ctr)
  {
    /* We need to retrieve the unaligned pointers, stored _before_ the actual fields. */
    void* ptr = ((void**)g_gauge_buffers.reserve[ctr])[-1];
    free(ptr);
  }
  free(g_gauge_buffers.reserve);
}
