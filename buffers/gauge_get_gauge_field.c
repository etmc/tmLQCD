#include "gauge.ih"

/* This routine not only malloc's a field, but immediately aligns it.
   To keep track of the original address to free the field eventually,
   we store that address _before_ the actual buffer.
   The end user should never have to see the alignment after this. */

gauge_field_t get_gauge_field()
{
  void *raw = NULL;
  size_t p = 0;
  gauge_field_t gauge_field;

  if (g_gauge_buffers.stack < 0) /* Need to allocate a new buffer */
  {
    if (g_gauge_buffers.allocated == g_gauge_buffers.max)
    {
#ifdef MPI
      MPI_Finalize(); /* Hard to handle gracefully here, can improve later... */
#endif
      exit(1);
    }
    ++g_gauge_buffers.stack;
    raw = malloc(sizeof(void*) + ALIGN_BASE + sizeof(su3_tuple) * VOLUMEPLUSRAND + 1);
    p = (size_t)raw + sizeof(void*);
    p = ((p + ALIGN_BASE) & ~ALIGN_BASE);
    ((void**)p)[-1] = raw;
    g_gauge_buffers.reserve[g_gauge_buffers.stack] = (su3_tuple*)p;
    ++g_gauge_buffers.allocated;
  }

  gauge_field.field = g_gauge_buffers.reserve[g_gauge_buffers.stack];
  g_gauge_buffers.reserve[g_gauge_buffers.stack] = NULL;
  --g_gauge_buffers.stack;

  return gauge_field;
}

