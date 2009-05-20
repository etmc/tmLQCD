#include <mpi.h>
#include <lemon.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

int main(int argc, char **argv)
{
  int rank;
  MPI_File fp;

  LemonWriter *w;
  LemonRecordHeader *h;

  int ME_flag=1, MB_flag=1, status=0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Writing */
  MPI_File_open(MPI_COMM_WORLD, "nfstalker.test", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);

  w = lemonCreateWriter(&fp, MPI_COMM_WORLD);

  h = lemonCreateHeader(MB_flag, ME_flag, "nfstalker-test", 128);

  fprintf(stderr, "Node %d: preparing to write...\n", w->my_rank);
  status = lemonWriteRecordHeader(h, w);

  fprintf(stderr, "Node %d: mopping up...\n", w->my_rank);
  lemonDestroyHeader(h);
  lemonWriterCloseRecord(w);
  lemonDestroyWriter(w);

  MPI_File_close(&fp);
  if (w->my_rank == 0)
    fprintf(stderr, "Finished writing!\n");
  MPI_Finalize();

  return(0);
}
