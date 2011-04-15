#include "params.ih"

paramsPropagatorFormat *construct_paramsPropagatorFormat(int const prec, int const flavours)
{
  paramsPropagatorFormat *format = malloc(sizeof(paramsPropagatorFormat));

  if (format == (paramsPropagatorFormat*)NULL)
    kill_with_error(NULL, g_cart_id, "Could not allocate paramsPropagatorFormat.");

  format->flavours = flavours;
  format->prec = prec;

  format->lx = LX * g_nproc_x;
  format->ly = LY * g_nproc_y;
  format->lz = LZ * g_nproc_z;
  format->lt =  T * g_nproc_t;

  return format;
}
