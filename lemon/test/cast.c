#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  int rank;
  uint32_t uint32;
  uint64_t uint64;
  size_t   size;
  MPI_Offset offset;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    size = 1;
    uint32 = size;
    uint64 = size;
    offset = size;
    fprintf(stderr, "Casting %llu gives:\n\tuint64: %llu\n\tuint32: %llu\n\toffset: %llu\n",
            (unsigned long long)size, (unsigned long long)uint64, (unsigned long long)uint32, (unsigned long long)offset);
  }

  MPI_Finalize();
  return 0;
}
