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
    allocate_gauge_buffers(1);

  gauge_field.field = g_gauge_buffers.reserve[g_gauge_buffers.stack];
  g_gauge_buffers.reserve[g_gauge_buffers.stack] = NULL;
  --g_gauge_buffers.stack;

  return gauge_field;
}

