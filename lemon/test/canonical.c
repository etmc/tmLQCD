#include <stdio.h>
#include <mpi.h>
#include <lemon.h>
#include <string.h>

int main(int argc, char **argv)
{
  int rank;
  MPI_File fp;
  MPI_Info info;
  LemonWriter *w;
  LemonRecordHeader *h;
  char message[] = "LEMON test message";
  uint64_t bytes = strlen(message);

  MPI_Init(&argc, &argv);
  MPI_Info_create(&info);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_File_open(MPI_COMM_WORLD, "canonical.test", MPI_MODE_WRONLY | MPI_MODE_CREATE, info, &fp);
  w = lemonCreateWriter(&fp, MPI_COMM_WORLD);
  h = lemonCreateHeader(1, 1, "lemon-test-text", bytes);

  lemonWriteRecordHeader(h, w);
  lemonWriteRecordData(message, &bytes, w);
  lemonDestroyHeader(h);
  lemonDestroyWriter(w);
  MPI_File_close(&fp);
  MPI_Finalize();
  return 0;
}
