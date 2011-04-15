#include "params.ih"

paramsSourceFormat *construct_paramsSourceFormat(int const prec, int const flavours, int const spins, int const colours)
{
  paramsSourceFormat *format = malloc(sizeof(paramsSourceFormat));

  if (format == (paramsSourceFormat*)NULL)
    kill_with_error(NULL, g_cart_id, "Could not allocate paramsSourceFormat.");

  format->prec = prec;
  format->flavours = flavours;

  format->lx = LX * g_nproc_x;
  format->ly = LY * g_nproc_y;
  format->lz = LZ * g_nproc_z;
  format->lt =  T * g_nproc_t;

  format->spins = spins;
  format->colours = colours;

  return format;
}
