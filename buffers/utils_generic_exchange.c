#include "utils.ih"

#ifndef MPI /*Let's deal with this case once and for all*/
void generic_exchange(void *field_in, int bytes_per_site)
{}
#else /* MPI */
void generic_exchange(void *field_in, int bytes_per_site)
{
#if defined _NON_BLOCKING
  int cntr=0;
  MPI_Request request[108];
  MPI_Status  status[108];
#else /* _NON_BLOCKING */
  MPI_Status status;
#endif /* _NON_BLOCKING */
  static int initialized = 0;

  /* We start by defining all the MPI datatypes required */
  static MPI_Datatype site_type;

  static MPI_Datatype slice_X_cont_type, slice_Y_cont_type, slice_Z_cont_type, slice_T_cont_type;
  static MPI_Datatype slice_X_subs_type, slice_Y_subs_type;
  static MPI_Datatype slice_X_gath_type, slice_Y_gath_type, slice_Z_gath_type;

  static MPI_Datatype edge_XY_cont_type, edge_XZ_cont_type, edge_XT_cont_type, edge_YZ_cont_type, edge_YT_cont_type, edge_ZT_cont_type;
  static MPI_Datatype edge_XY_gath_type, edge_XZ_gath_type, edge_XT_gath_type, edge_YZ_gath_type, edge_YT_gath_type, edge_ZT_gath_type;

  unsigned char(*buffer)[bytes_per_site] = field_in; /* To allow for pointer arithmetic */

  // To avoid continuous MPI operations on these local variables, let's declare them static.
  // That means we should only initialize if this is the first use of the function, or if
  // the existing initialization is for the wrong number of bytes per size!
  if (initialized && (initialized != bytes_per_site))
  {
    MPI_Type_free(&site_type);

    MPI_Type_free(&slice_T_cont_type);
    MPI_Type_free(&slice_X_cont_type);
    MPI_Type_free(&slice_Y_cont_type);
    MPI_Type_free(&slice_Z_cont_type);

    MPI_Type_free(&slice_X_subs_type);
    MPI_Type_free(&slice_Y_subs_type);

    MPI_Type_free(&slice_X_gath_type);
    MPI_Type_free(&slice_Y_gath_type);
    MPI_Type_free(&slice_Z_gath_type);

    MPI_Type_free(&edge_XY_cont_type);
    MPI_Type_free(&edge_XZ_cont_type);
    MPI_Type_free(&edge_XT_cont_type);
    MPI_Type_free(&edge_YZ_cont_type);
    MPI_Type_free(&edge_YT_cont_type);
    MPI_Type_free(&edge_ZT_cont_type);

    MPI_Type_free(&edge_XY_gath_type);
    MPI_Type_free(&edge_XZ_gath_type);
    MPI_Type_free(&edge_XT_gath_type);
    MPI_Type_free(&edge_YZ_gath_type);
    MPI_Type_free(&edge_YT_gath_type);
    MPI_Type_free(&edge_ZT_gath_type);

    /* We're ready to reinitialize all these types now... */
    initialized = 0;
  }

  if (!initialized)
  {
    /* Initialization of the datatypes - adapted from mpi_init.c */
    MPI_Type_contiguous(bytes_per_site, MPI_BYTE, &site_type);
    MPI_Type_commit(&site_type);

    MPI_Type_contiguous(LX * LY *LZ, site_type, &slice_T_cont_type);
    MPI_Type_contiguous( T * LY *LZ, site_type, &slice_X_cont_type);
    MPI_Type_contiguous( T * LX *LZ, site_type, &slice_Y_cont_type);
    MPI_Type_contiguous( T * LX *LY, site_type, &slice_Z_cont_type);

    MPI_Type_commit(&slice_T_cont_type);
    MPI_Type_commit(&slice_X_cont_type);
    MPI_Type_commit(&slice_Y_cont_type);
    MPI_Type_commit(&slice_Z_cont_type);

    MPI_Type_contiguous(LY * LZ, site_type, &slice_X_subs_type);
    MPI_Type_contiguous(LZ, site_type, &slice_Y_subs_type);

    MPI_Type_commit(&slice_X_subs_type);
    MPI_Type_commit(&slice_Y_subs_type);

    MPI_Type_vector(T, 1, LX, slice_X_subs_type, &slice_X_gath_type);
    MPI_Type_vector(T * LX, 1, LY, slice_Y_subs_type, &slice_Y_gath_type);
    MPI_Type_vector(T * LX * LY, 1, LZ, site_type,  &slice_Z_gath_type);

    MPI_Type_commit(&slice_X_gath_type);
    MPI_Type_commit(&slice_Y_gath_type);
    MPI_Type_commit(&slice_Z_gath_type);

    MPI_Type_contiguous(2 * T * LZ, site_type, &edge_XY_cont_type);
    MPI_Type_contiguous(2 * T * LY, site_type, &edge_XZ_cont_type);
    MPI_Type_contiguous(2 * LY * LZ, site_type, &edge_XT_cont_type);
    MPI_Type_contiguous(2 * T * LX, site_type, &edge_YZ_cont_type);
    MPI_Type_contiguous(2 * LX * LZ, site_type, &edge_YT_cont_type);

    MPI_Type_contiguous(2 * LX * LY, site_type, &edge_ZT_cont_type);

    MPI_Type_commit(&edge_XY_cont_type);
    MPI_Type_commit(&edge_XZ_cont_type);
    MPI_Type_commit(&edge_XT_cont_type);
    MPI_Type_commit(&edge_YZ_cont_type);
    MPI_Type_commit(&edge_YT_cont_type);
    MPI_Type_commit(&edge_ZT_cont_type);

    MPI_Type_vector(2 * T, LZ, LX * LZ, site_type, &edge_XY_gath_type);
    MPI_Type_vector(2 * T, LY, LY * LX, site_type, &edge_XZ_gath_type);
    MPI_Type_vector(2, 1, T, slice_X_subs_type, &edge_XT_gath_type);
    MPI_Type_vector(2 * T * LX, 1, LY, site_type, &edge_YZ_gath_type);
    MPI_Type_vector(2 * LX, LZ, LY * LZ, site_type, &edge_YT_gath_type);
    MPI_Type_vector(2 * LX * LY, 1, LZ, site_type, &edge_ZT_gath_type);

    MPI_Type_commit(&edge_XY_gath_type);
    MPI_Type_commit(&edge_XZ_gath_type);
    MPI_Type_commit(&edge_XT_gath_type);
    MPI_Type_commit(&edge_YZ_gath_type);
    MPI_Type_commit(&edge_YT_gath_type);
    MPI_Type_commit(&edge_ZT_gath_type);

    initialized = bytes_per_site;
  }

  /* Following are implementations using different compile time flags */
#if defined _NON_BLOCKING
# if defined _INDEX_INDEP_GEOM
#  include "utils_generic_exchange.1.inc"
# else /* _INDEX_INDEP_GEOM */
#  include "utils_generic_exchange.2.inc"
# endif /* _INDEX_INDEP_GEOM */
#else /* _NON_BLOCKING */
# if defined _INDEX_INDEP_GEOM
#  include "utils_generic_exchange.3.inc"
# else /* _INDEX_INDEP_GEOM */
#  include "utils_generic_exchange.4.inc"
# endif /* _INDEX_INDEP_GEOM */
#endif /* _NON_BLOCKING */
}

#endif /* MPI */

