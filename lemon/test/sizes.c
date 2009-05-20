#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  int rank;
  int size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    fprintf(stderr, "Reporting on sizes for different datatypes.\n");

    fprintf(stderr, "\tSize of char..................: %lu.\n", (long unsigned int)sizeof(char));
    fprintf(stderr, "\tSize of int...................: %lu.\n", (long unsigned int)sizeof(int));
    fprintf(stderr, "\tSize of long int..............: %lu.\n", (long unsigned int)sizeof(long int));
    fprintf(stderr, "\tSize of long long int.........: %lu.\n", (long unsigned int)sizeof(long long int));
    fprintf(stderr, "\tSize of size_t................: %lu.\n", (long unsigned int)sizeof(size_t));
    fprintf(stderr, "\tSize of MPI_Offset............: %lu.\n\n", (long unsigned int)sizeof(MPI_Offset));

    fprintf(stderr, "\tSize of unsigned char.........: %lu.\n", (long unsigned int)sizeof(unsigned char));
    fprintf(stderr, "\tSize of unsigned int..........: %lu.\n", (long unsigned int)sizeof(unsigned int));
    fprintf(stderr, "\tSize of unsigned long int.....: %lu.\n", (long unsigned int)sizeof(unsigned long int));
    fprintf(stderr, "\tSize of unsigned long long int: %lu.\n\n", (long unsigned int)sizeof(unsigned long long int));

    fprintf(stderr, "\tSize of int8_t................: %lu.\n", (long unsigned int)sizeof(int8_t));
    fprintf(stderr, "\tSize of int16_t...............: %lu.\n", (long unsigned int)sizeof(int16_t));
    fprintf(stderr, "\tSize of int32_t...............: %lu.\n", (long unsigned int)sizeof(int32_t));
    fprintf(stderr, "\tSize of int64_t...............: %lu.\n\n", (long unsigned int)sizeof(int64_t));

    fprintf(stderr, "\tSize of uint8_t...............: %lu.\n", (long unsigned int)sizeof(uint8_t));
    fprintf(stderr, "\tSize of uint16_t..............: %lu.\n", (long unsigned int)sizeof(uint16_t));
    fprintf(stderr, "\tSize of uint32_t..............: %lu.\n", (long unsigned int)sizeof(uint32_t));
    fprintf(stderr, "\tSize of uint64_t..............: %lu.\n\n", (long unsigned int)sizeof(uint64_t));

    fprintf(stderr, "\tSize of float.................: %lu.\n", (long unsigned int)sizeof(float));
    fprintf(stderr, "\tSize of double................: %lu.\n", (long unsigned int)sizeof(double));
    fprintf(stderr, "\tSize of long double...........: %lu.\n\n", (long unsigned int)sizeof(long double));

    fprintf(stderr, "Reporting on sizes for MPI types.\n");
    MPI_Type_size(MPI_BYTE, &size);
    fprintf(stderr, "\tSize of MPI_BYTE..............: %d.\n", size);
    MPI_Type_size(MPI_PACKED, &size);
    fprintf(stderr, "\tSize of MPI_PACKED............: %d.\n\n", size);

    MPI_Type_size(MPI_CHAR, &size);
    fprintf(stderr, "\tSize of MPI_CHAR..............: %d.\n", size);
    MPI_Type_size(MPI_WCHAR, &size);
    fprintf(stderr, "\tSize of MPI_WCHAR.............: %d.\n\n", size);

    MPI_Type_size(MPI_SIGNED_CHAR, &size);
    fprintf(stderr, "\tSize of MPI_SIGNED_CHAR.......: %d.\n", size);
    MPI_Type_size(MPI_SHORT, &size);
    fprintf(stderr, "\tSize of MPI_SHORT.............: %d.\n", size);
    MPI_Type_size(MPI_INT, &size);
    fprintf(stderr, "\tSize of MPI_INT...............: %d.\n", size);
    MPI_Type_size(MPI_LONG, &size);
    fprintf(stderr, "\tSize of MPI_LONG..............: %d.\n", size);
    MPI_Type_size(MPI_LONG_LONG, &size);
    fprintf(stderr, "\tSize of MPI_LONG_LONG.........: %d.\n\n", size);

    MPI_Type_size(MPI_UNSIGNED_CHAR, &size);
    fprintf(stderr, "\tSize of MPI_UNSIGNED_CHAR.....: %d.\n", size);
    MPI_Type_size(MPI_UNSIGNED_SHORT, &size);
    fprintf(stderr, "\tSize of MPI_UNSIGNED_SHORT....: %d.\n", size);
    MPI_Type_size(MPI_UNSIGNED, &size);
    fprintf(stderr, "\tSize of MPI_UNSIGNED..........: %d.\n", size);
    MPI_Type_size(MPI_UNSIGNED_LONG, &size);
    fprintf(stderr, "\tSize of MPI_UNSIGNED_LONG.....: %d.\n", size);
    MPI_Type_size(MPI_UNSIGNED_LONG_LONG, &size);
    fprintf(stderr, "\tSize of MPI_UNSIGNED_LONG_LONG: %d.\n\n", size);

    MPI_Type_size(MPI_FLOAT, &size);
    fprintf(stderr, "\tSize of MPI_FLOAT.............: %d.\n", size);
    MPI_Type_size(MPI_DOUBLE, &size);
    fprintf(stderr, "\tSize of MPI_DOUBLE............: %d.\n", size);
    MPI_Type_size(MPI_LONG_DOUBLE, &size);
    fprintf(stderr, "\tSize of MPI_LONG_DOUBLE.......: %d.\n\n", size);

    MPI_Type_size(MPI_INTEGER1, &size);
    fprintf(stderr, "\tSize of MPI_INTEGER1..........: %d.\n", size);
    MPI_Type_size(MPI_INTEGER2, &size);
    fprintf(stderr, "\tSize of MPI_INTEGER2..........: %d.\n", size);
    MPI_Type_size(MPI_INTEGER4, &size);
    fprintf(stderr, "\tSize of MPI_INTEGER4..........: %d.\n", size);
    MPI_Type_size(MPI_INTEGER8, &size);
    fprintf(stderr, "\tSize of MPI_INTEGER8..........: %d.\n\n", size);

    MPI_Type_size(MPI_REAL4, &size);
    fprintf(stderr, "\tSize of MPI_REAL4.............: %d.\n", size);
    MPI_Type_size(MPI_REAL8, &size);
    fprintf(stderr, "\tSize of MPI_REAL8.............: %d.\n", size);
  }

  MPI_Finalize();
  return 0;
}
