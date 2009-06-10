#include "params.ih"

paramsIldgFormat *construct_paramsIldgFormat(int const prec)
{
  paramsIldgFormat *format = malloc(sizeof(paramsIldgFormat));

  if (format == (paramsIldgFormat*)NULL)
    kill_with_error(NULL, g_cart_id, "Could not allocate paramsIldgFormat.");

  format->prec = prec;
  format->nx = LX * g_nproc_x;
  format->ny = LY * g_nproc_y;
  format->nz = LZ * g_nproc_z;
  format->nt =  T * g_nproc_t;

  return format;
}
