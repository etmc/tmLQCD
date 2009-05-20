#include <config.h>
#include <lemon.h>

MPI_Offset lemonGetReaderPointer(LemonReader *reader)
{
  return (reader->off + reader->pos);
}
