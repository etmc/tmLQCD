#include "utils.ih"

#ifdef HAVE_LIBLEMON
void kill_with_error(MPI_File *fh, int const rank, char const *error)
{
  if (rank == 0)
    fprintf(stderr, "%s", error);
  MPI_File_close(fh);
  MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Finalize();
  exit(500);
}
#else
void kill_with_error(FILE *fh, int const rank, char const *error) {
  if (rank == 0)
    fprintf(stderr, "%s", error);
  fclose(fh);
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Finalize();
#endif
  exit(500);
}
#endif
