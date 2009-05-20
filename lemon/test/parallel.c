#include <lemon.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv)
{
  MPI_File fp;

  LemonWriter *w;
  LemonReader *r;
  LemonRecordHeader *h;

  char *data;
  char *data_read;
  int mpisize;
  int rank;
  char *type;

  int ME_flag=1, MB_flag=1, status=0;

  int latDist[] = {0, 0, 0, 0};
  int periods[] = {1, 1, 1, 1};
  /* The following are the local volumes - we extend the lattice as needed. */
  int latSizes[] = {4, 4, 4, 4};
  int latVol = 256;

  MPI_Comm cartesian;
  int i;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  MPI_Dims_create(mpisize, 4, latDist);
  MPI_Cart_create(MPI_COMM_WORLD, 4, latDist, periods, 1, &cartesian);
  for (i = 0; i < 4; ++i)
    latSizes[i] *= latDist[i];
  latVol *= mpisize;
  MPI_Comm_rank(cartesian, &rank);

   /* Start of code - writing */
  MPI_File_open(cartesian, "parallel.test", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);
  w = lemonCreateWriter(&fp, cartesian);

  data = (char*)malloc(257);
  for (i = 0; i < latVol/mpisize; ++i)
    data[i] = (char)(rank + 48);
  data[256] = '\0';

  h = lemonCreateHeader(MB_flag, ME_flag, "parallel-test", latVol);
  status = lemonWriteRecordHeader(h, w);

  lemonDestroyHeader(h);

  lemonWriteLatticeParallel(w, data, sizeof(char), latSizes);

  lemonWriterCloseRecord(w);
  lemonDestroyWriter(w);
  MPI_File_close(&fp);

  /* Reading */
  data_read = (char*)malloc(257);

  MPI_File_open(cartesian, "parallel.test", MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);
  r = lemonCreateReader(&fp, cartesian);

  if (lemonReaderNextRecord(r))
    fprintf(stderr, "Node %d reports: next record failed.\n", rank);

  type = lemonReaderType(r);
  if (strncmp(type, "parallel-test", 13))
    fprintf(stderr, "Node %d reports: wrong type read.\n", rank);

  lemonReadLatticeParallel(r, data_read, sizeof(char), latSizes);
  data_read[256] = '\0';
  if (strncmp(data_read, data, 256))
  {
    fprintf(stderr, "Node %d reports: wrong data read.\n", rank);
    fprintf(stderr, "Node %d wanted: %s\n", rank, data);
    fprintf(stderr, "Node %d got   : %s\n", rank, data_read);
  }
  else
    fprintf(stderr, "Node %d reports data okay.\n", rank);

  lemonDestroyReader(r);
  MPI_File_close(&fp);
  MPI_Finalize();

  free(data);
  free(data_read);

  return(0);
}

