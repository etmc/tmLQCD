#include "utils.ih"

void kill_with_error(LIME_FILE *fh, int const rank, char const *error)
{
  if (rank == 0)
    fprintf(stderr, "%s", error);
#ifdef HAVE_LIBLEMON
  MPI_File_close(fh);
#else
  fclose(fh);
#endif /* HAVE_LIBLEMON */
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Finalize();
#endif
  exit(500);
}
