#include "utils.ih"

void  print_su3adj(su3adj const *in)
{
  double const *in_as_array = (double const *)(in);
  printf("[");
  for (int ctr = 0; ctr < 8; ++ctr)
  {
    printf("%c %14.12f", (in_as_array[ctr] > 0 ? '+' : '-'), fabs(in_as_array[ctr]));
    if (ctr < 7)
      printf(", ");
    else
      printf("];\n");
  }
}