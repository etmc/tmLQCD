#include "params.ih"

paramsIldgFormat *construct_paramsIldgFormat(int const prec)
{
  paramsIldgFormat *format = malloc(sizeof(paramsIldgFormat));

  if (format == (paramsIldgFormat*)NULL)
    kill_with_error(NULL, g_cart_id, "Could not allocate paramsIldgFormat.");

  format->prec = prec;
  format->lx = L;
  format->ly = L;
  format->lz = L;
  format->lt = T_global;

  return format;
}
