#include "utils.ih"

void kill_with_error(MPI_File *fh, int const rank, char const *error)
{
  if (rank == 0)
    fprintf(stderr, "%s", error);
  MPI_File_close(fh);
  MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Finalize();
  exit(500);
}
