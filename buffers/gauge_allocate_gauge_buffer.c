#include "gauge.ih"

void allocate_gauge_buffer(unsigned int count)
{
  if ((g_gauge_buffers.allocated + count) > g_gauge_buffers.max)
  {
    exit(1); /* Hard to handle gracefully here, can improve later... */
    /* fatal_error("Maximum number of allocated gauge fields exceeded."); */
  }
  
  for (unsigned int ctr = 0; ctr < count; ++ctr)
  {
    raw = malloc(sizeof(void*) + ALIGN_BASE + sizeof(su3_tuple) * VOLUMEPLUSRAND + 1);
    /*
     * if (raw == NULL)
     *   fatal_error("Could not allocate the requested amount of memory.");
     */
    p = (size_t)raw + sizeof(void*);
    p = ((p + ALIGN_BASE) & ~ALIGN_BASE);
    ((void**)p)[-1] = raw;

    ++g_gauge_buffers.stack;
    g_gauge_buffers.reserve[g_gauge_buffers.stack] = (su3_tuple*)p;
    ++g_gauge_buffers.allocated;
  }
}