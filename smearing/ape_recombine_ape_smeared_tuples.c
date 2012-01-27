#include "ape.ih"

void recombine_ape_smeared_tuples(gauge_field_t smeared, gauge_field_t original)
{
  _Complex double *sp = (_Complex double*)&smeared.field[0][0];
  _Complex double *op = (_Complex double*)&original.field[0][0];
  su3 *tp = (su3*)smeared.field[0][0];

  for (int idx = 0; idx < 9 * 4 * VOLUME; ++idx)
    op[idx] = rho_s * sp[idx] + rho_p * op[idx];
  
  for (int idx = 0; idx <     4 * VOLUME; ++idx)
    reunitarize(tp + idx);
}
