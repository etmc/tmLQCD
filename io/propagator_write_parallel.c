#include "propagator.ih"

int write_spinor_parallel(spinor * const s, spinor * const r, char * filename,
                          const int append, const int prec)
{
  MPI_File ofs;
  LemonWriter *lemonwriter = NULL;
  LemonRecordHeader *lemonheader = NULL;
  int status = 0;
  int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
  n_uint64_t bytes;
  DML_Checksum checksum;
  int err = 0;

  if(g_cart_id == 0)
  {
    if(append)
      amode |= MPI_MODE_APPEND;

    err = MPI_File_open(g_cart_grid, filename, amode, MPI_INFO_NULL, &ofs);

    if (err != MPI_SUCCESS)
    {
      fprintf(stderr, "Could not open file %s for writing!\n Aborting...\n", filename);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }

    lemonwriter = lemonCreateWriter(&ofs, g_cart_grid);
    if (lemonwriter == (LemonWriter*)NULL)
    {
      fprintf(stderr, "LEMON error in file %s for writing!\n Aborting...\n", filename);
      MPI_File_close(&ofs);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }

    bytes = (n_uint64_t)LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*(n_uint64_t)(sizeof(spinor)*prec/64);
    lemonheader = lemonCreateHeader(1, 0, "scidac-binary-data", bytes);
    status = lemonWriteRecordHeader( lemonheader, lemonwriter);
    if(status != LEMON_SUCCESS )
    {
      fprintf(stderr, "LEMON write header (scidac-binary-data) error %d\n", status);
      MPI_File_close(&ofs);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }
    lemonDestroyHeader(lemonheader);
  }

  status = write_binary_spinor_data_parallel(s, r, lemonwriter, &checksum, prec);
  if(status != LEMON_SUCCESS )
  {
    fprintf(stderr, "LEMON write error %d\n", status);
    MPI_File_close(&ofs);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }

  status = write_checksum_parallel(lemonwriter, &checksum);
  if(status != LEMON_SUCCESS )
  {
    fprintf(stderr, "LEMON write error %d\n", status);
    MPI_File_close(&ofs);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }

  return(0);
}
