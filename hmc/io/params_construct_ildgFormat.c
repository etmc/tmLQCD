#include "params.ih"

paramsIldgFormat *construct_paramsIldgFormat(int const prec)
{
  paramsIldgFormat *format = malloc(sizeof(paramsIldgFormat));

  if (format == (paramsIldgFormat*)NULL)
    kill_with_error(NULL, g_cart_id, "Could not allocate paramsIldgFormat.");

  format->prec = prec;
  format->nx = L;
  format->ny = L;
  format->nz = L;
  format->nt = T_global;

  return format;
}
