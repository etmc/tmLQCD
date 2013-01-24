#include "utils.ih"

void unfold_field(gauge_field_t *target, gauge_field_t const base)
{
  su3 ALIGN tmp;
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
      _su3d_times_su3(tmp, base[x][mu], (*target)[x][mu]);
      memmove(&(*target)[x][mu], &tmp, sizeof(su3));
    }
}
