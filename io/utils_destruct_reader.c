#include "utils.ih"

void destruct_reader(READER * reader)
{
  LIME_FILE *fh = NULL;

  fh = reader->fp;
  DestroyReader(reader);
#ifdef HAVE_LIBLEMON
  MPI_File_close(fh);
  free(fh); /* NB This assumes construct_writer was used to malloc memory! */
#else /* HAVE_LIBLEMON */
  fclose(fh);
#endif /* HAVE_LIBLEMON */
}
