#include <config.h>
#include <lemon.h>

int lemonDestroyWriter(LemonWriter *writer)
{
  if (writer == (LemonWriter*)NULL)
    return LEMON_SUCCESS;

  MPI_Comm_free(&writer->cartesian);

  free(writer);
  return LEMON_SUCCESS;
}
