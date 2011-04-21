#include "utils.ih"

void construct_writer(WRITER ** writer, char * filename, const int append)
{
  LIME_FILE *fh = NULL;
  int status = 0;
  if(g_debug_level > 0 && g_cart_id == 0) {
#ifdef HAVE_LIBLEMON
    printf("# Constructing LEMON writer for file %s for append = %d\n", filename, append);
#else
    printf("# Constructing LIME writer for file %s for append = %d\n", filename, append);
#endif
  }

#ifdef HAVE_LIBLEMON
  fh = (MPI_File*)malloc(sizeof(MPI_File));
  if(append) {
    status = MPI_File_open(g_cart_grid, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, fh);
  }
  else {
    status = MPI_File_open(g_cart_grid, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, fh);
    if(status == MPI_SUCCESS) status = MPI_File_set_size(*fh, 0);
  }
  status = (status == MPI_SUCCESS) ? 0 : 1;
  *writer = lemonCreateWriter(fh, g_cart_grid);
  status = status || (writer == NULL);
#else /* HAVE_LIBLEMON */
  if (g_cart_id == 0)
  {
    if(append) {
      fh = fopen(filename, "a");
    }
    else {
      fh = fopen(filename, "w");
    }
    status = (fh == NULL);
    *writer = limeCreateWriter(fh);
    status = status || (writer == NULL);
  }
#endif /* HAVE_LIBLEMON */

  if (status)
    kill_with_error(fh, g_cart_id, "Failed to create writer. Aborting...\n");
}
