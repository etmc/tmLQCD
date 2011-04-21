#include "utils.ih"

void  print_config_to_screen(su3 **in) 
{
  int x, mu;
  for(x = 0; x < VOLUME; x++)
    for(mu = 0; mu < 4; mu++)
    {
      printf("x = %d  mu = %d\n", x, mu);
      /*print_su3_full_hex_precision(&(in[x][mu]));*/
      print_su3(&(in[x][mu]));
    }
}
