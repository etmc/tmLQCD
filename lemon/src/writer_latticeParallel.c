#include <config.h>
#include <lemon.h>

int lemonWriteLatticeParallel(LemonWriter *writer, void *data,
                             MPI_Offset siteSize, int *latticeDims)
{
  MPI_Datatype etype;
  MPI_Datatype ftype;
  MPI_Status status;

  int idx = 0;
  int ndims = 0;
  int localVol = 1;
  int totalVol = 1;
  int *starts;
  int *localDims;
  int *mpiDims;
  int *mpiCoords;
  int *period;

  /* Elementary datatype for reading/writing: a single lattice site of the correct size.
     Note that we always assume we obtained this as raw binary data of the correct type,
     because apparently we cannot trust MPI libraries to deal with conversions internally. */
  MPI_Type_contiguous(siteSize, MPI_BYTE, &etype);
  MPI_Type_commit(&etype);

  /* Gathering of the required MPI data from the cartesian communicator. */
  MPI_Cartdim_get(writer->cartesian, &ndims);
  starts = (int*)malloc(ndims * sizeof(int));
  localDims = (int*)malloc(ndims * sizeof(int));
  mpiDims = (int*)malloc(ndims * sizeof(int));
  mpiCoords = (int*)malloc(ndims * sizeof(int));
  period = (int*)malloc(ndims * sizeof(int));
  MPI_Cart_get(writer->cartesian, ndims, mpiDims, period, mpiCoords);
  free(period);

  /* Calculation of local lattice dimensions from the MPI data we obtained. */
  for (idx = 0; idx < ndims; ++idx)
  {
    localDims[idx] = latticeDims[idx] / mpiDims[idx];
    localVol *= localDims[idx];
    totalVol *= latticeDims[idx];
    starts[idx] = localDims[idx] * mpiCoords[idx];
  }

  /* Build up a filetype that provides the right offsets for the writing of a N-dimensional lattice. */
  MPI_Type_create_subarray(ndims, latticeDims, localDims, starts, MPI_ORDER_C, etype, &ftype);
  MPI_Type_commit(&ftype);

  /* Install the data organization we worked out above on the file as a view */
  /* We want to start writing wherever the latest write occurred (because of the
     no seeking in LemonWriter files rule this should be consistent). We therefore
     take the global maximum of all MPI_long File filepointers in absolute offsets and
     start our writeout operation from there. */
  MPI_Barrier(writer->cartesian);
  MPI_File_set_view(*writer->fh, writer->off + writer->pos, etype, ftype, "native", MPI_INFO_NULL);

  /* Blast away! */
  MPI_File_write_at_all(*writer->fh, writer->pos, data, localVol, etype, &status);
  MPI_File_sync(*writer->fh);

  MPI_Barrier(writer->cartesian);

  writer->pos += totalVol * siteSize;

  /* We should reset the shared file pointer, in an MPI_BYTE based view... */
  MPI_Barrier(writer->cartesian);
  MPI_File_set_view(*writer->fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

  /* Free up the resources we claimed for this operation. */
  MPI_Type_free(&etype);
  MPI_Type_free(&ftype);
  free(starts);
  free(localDims);
  free(mpiDims);
  free(mpiCoords);

  return LEMON_SUCCESS;
}
