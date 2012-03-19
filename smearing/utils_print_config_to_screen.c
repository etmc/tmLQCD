#include "utils.ih"

void  print_config_to_screen(su3 **in) 
{
  for(int x = 0; x < VOLUME; ++x)
    for(int mu = 0; mu < 4; ++mu)
    {
      printf("x = %d  mu = %d\n", x, mu);
      print_su3(&(in[x][mu]));
    }
}
