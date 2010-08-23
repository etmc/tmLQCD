#include "utils.ih"

void destruct_writer(WRITER * writer)
{
  LIME_FILE *fh = NULL;

#ifdef HAVE_LIBLEMON
  fh = writer->fp;
  lemonDestroyWriter(writer);
  MPI_File_close(fh);
  free(fh); /* NB This assumes construct_writer was used to malloc memory! */
#else /* HAVE_LIBLEMON */
  if (g_cart_id == 0)
  {
    fh = writer->fp;
    limeDestroyWriter(writer);
    fclose(fh);
  }
#endif /* HAVE_LIBLEMON */
}
