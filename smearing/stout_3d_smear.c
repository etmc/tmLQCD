#include "stout_3d.ih"

#include "stout_fatten_links.static"

void stout_3d_smear(stout_3d_control *control, gauge_field_t in)
{
  su3 ALIGN staples;
  /* Staple calculation isn't local, so we cannot store changes per site immediately
   * and need a separate input and output. Since the result of one iteration becomes 
   * the input of the next, that means we'll actually need an additional field of 
   * buffer space. */
  gauge_field_t buffer = get_gauge_field_annotated("STOUT 3D BUFFER");

  control->U[0] = in; /* Shallow copy intended */
  
#pragma omp parallel private(staples)
  for(unsigned int iter = 0; iter < control->iterations; ++iter)
  {
#pragma omp for
    for (unsigned int x = 0; x < VOLUME; ++x)
    {
      _su3_assign(buffer[x][0], in[x][0]); // Left untouched, but still needed for future calculations!
      for (unsigned int mu = 1; mu < 4; ++mu)
      {
        generic_staples_3d(&staples, x, mu, in);
        fatten_links(&buffer[x][mu], control->coeff, &staples, &in[x][mu]);
      }
    }
    /* Prepare for the next iteration, swap and exchange fields */
    /* There should be an implicit OMP barrier here */
#pragma omp single
    {
      swap_gauge_field(&control->U[1], &buffer);
      exchange_gauge_field(&control->U[1]); /* The edge terms will have changed by now */
      in = control->U[1]; /* Shallow copy intended */
    }
  }
  
  control->result = control->U[1];
  return_gauge_field(&buffer);
}
