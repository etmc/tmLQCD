#include "utils.ih"

void rnd_gauge_trafo(gauge_field_t * target, gauge_field_t const src)
{
  /* Generate a global gauge transformation */
  su3 *gauge_trafo = aalloc(VOLUMEPLUSRAND * sizeof(su3));
  for (unsigned int ix = 0; ix < VOLUME; ++ix)
  {
    random_su3(&gauge_trafo[ix]);
  }
  generic_exchange(gauge_trafo, sizeof(su3));
  
  /* Apply the global gauge transformation */
  su3 ALIGN tmp;
  for (unsigned int ix = 0; ix < VOLUME; ++ix)
  {
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
      _su3_times_su3d(tmp, src[ix][mu], gauge_trafo[g_iup[ix][mu]]);
      _su3_times_su3((*target)[ix][mu], gauge_trafo[ix], tmp);
    }
  }
      
  afree(gauge_trafo);
}

