#include "utils.ih"

void construct_writer(WRITER ** writer, char const *filename)
{
  LIME_FILE *fh = NULL;
  int status = 0;

#ifdef HAVE_LIBLEMON
  fh = (MPI_File*)malloc(sizeof(MPI_File));
  status = MPI_File_open(g_cart_grid, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, fh);
  status = (status == MPI_SUCCESS) ? 0 : 1;
  *writer = lemonCreateWriter(fh, g_cart_grid);
  status = status || (writer == NULL);
#else /* HAVE_LIBLEMON */
  if (g_cart_id == 0)
  {
    fh = fopen(filename, "a");
    status = (fh == NULL);
    *writer = limeCreateWriter(fh);
    status = status || (writer == NULL);
  }
#endif /* HAVE_LIBLEMON */

  if (status)
    kill_with_error(fh, g_cart_id, "Failed to create writer. Aborting...\n");
}
