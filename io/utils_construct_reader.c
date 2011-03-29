#include "utils.ih"

void construct_reader(READER ** reader, char * filename)
{
  LIME_FILE *fh = NULL;
  int status = 0;

  if(g_debug_level > 0 && g_cart_id == 0) {
#ifdef HAVE_LIBLEMON
    printf("# Constructing LEMON reader for file %s ...\n", filename);
#else
    printf("# Constructing LIME reader for file %s ...\n", filename);
#endif
  }


#ifdef HAVE_LIBLEMON
  fh = (MPI_File*)malloc(sizeof(MPI_File));
  status = MPI_File_open(g_cart_grid, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh);
  status = (status == MPI_SUCCESS) ? 0 : 1;
#else /* HAVE_LIBLEMON */
  fh = fopen(filename, "r");
  status = (fh == NULL) ? 1 : 0;
  fflush(stderr);
#endif /* HAVE_LIBLEMON */

  if (status) {
    kill_with_error(fh, g_cart_id, "\nUnable to open file for reading.\nPlease verify file existence and access rights.\nUnable to continue.\n");
  }

#ifdef HAVE_LIBLEMON
  *reader = lemonCreateReader(fh, g_cart_grid);
#else /* HAVE_LIBLEMON */
  *reader = limeCreateReader(fh);
#endif /* HAVE_LIBLEMON */

  if (*reader == (READER *)NULL) {
    kill_with_error(fh, g_cart_id, "\nCould not create reader, unable to continue.\n");
  }
}
