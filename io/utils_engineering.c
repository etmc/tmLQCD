#include "utils.ih"

static char prefix[] = {'z', 'a', 'f', 'p', 'u', 'm', ' ', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'};

void engineering(char *result, double value, char const *units)
{
  double logval = log10(value);
  int logscale;
  int digits;

  logscale = (int)floor(logval/3);

  if (logscale > -6 && logscale < 6)
  {
    value /= pow(1E3, (double)logscale);
    digits = 2 - (int)log10(value);

    sprintf(result, "%.*f %c%s", digits, value, prefix[logscale + 6], units);
  }
  else
  {
    sprintf(result, "%4.2e %s", value, units);
  }
}
