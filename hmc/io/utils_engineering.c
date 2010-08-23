#include "utils.ih"

static char prefix[] = {'z', 'a', 'f', 'p', 'u', 'm', ' ', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'};

void engineering(char *result, double value, char const *units)
{
  double logval = log10(value);
  int logscale;
  int digits = 2;

  logscale = (int)floor(logval / 3);

  if (logscale > -6 && logscale < 6)
  {
    value /= pow(1E3, (double)logscale);
    if (value > 100)
      digits = 0;
    else
      if (value > 10)
        digits = 1;

    if (logscale)
      sprintf(result, "%.*f %c%s", digits, value, prefix[logscale + 6], units);
    else
      sprintf(result, "%.*f %s", digits, value, units);
  }
  else
  {
    sprintf(result, "%4.2e %s", value, units);
  }
}
