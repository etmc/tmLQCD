#include "utils.ih"

void destruct_reader(READER *reader) {
  LIME_FILE *fh = NULL;

  fh = reader->fp;
  DestroyReader(reader);
#ifdef TM_USE_LEMON
  MPI_File_close(fh);
  free(fh); /* NB This assumes construct_writer was used to malloc memory! */
#else       /* TM_USE_LEMON */
  fclose(fh);
#endif      /* TM_USE_LEMON */
}
