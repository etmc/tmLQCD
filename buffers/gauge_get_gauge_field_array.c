#include "gauge.ih"

gauge_field_array_t get_gauge_field_array(unsigned int length)
{
  gauge_field_array_t gauge_field_array;
  gauge_field_array.length = length;
  gauge_field_array.field_array = (su3_tuple**)calloc(length, sizeof(su3_tuple*));

  if (g_gauge_buffers.stack < (length - 1)) /* Need to allocate more buffers */
  {
    if (g_gauge_buffers.allocated >= (g_gauge_buffers.max - length))
    {
#ifdef MPI
      MPI_Finalize();
#endif
      exit(1); /* Hard to handle gracefully here, can improve later... */
    }
    for (unsigned int ctr = 0; ctr < length; ++ctr)
    {
      ++g_gauge_buffers.stack;
      g_gauge_buffers.reserve[g_gauge_buffers.stack] = (su3_tuple*)malloc(sizeof(su3_tuple) * VOLUMEPLUSRAND + 1);
    }
    g_gauge_buffers.allocated += length;
  }

  for (unsigned int ctr = 0; ctr < length; ++ctr)
  {
    gauge_field_array.field_array[ctr] = g_gauge_buffers.reserve[g_gauge_buffers.stack];
    g_gauge_buffers.reserve[g_gauge_buffers.stack] = NULL;
    --g_gauge_buffers.stack;
  }
  return gauge_field_array;
}

