#include "params.ih"

paramsXlfInfo *construct_paramsXlfInfo(double const plaq, int const counter)
{
  struct timeval t1;
  paramsXlfInfo *info = malloc(sizeof(paramsXlfInfo));

  if (info == (paramsXlfInfo*)NULL)
    kill_with_error(NULL, g_cart_id, "Could not allocate paramsXlfInfo.");

  gettimeofday(&t1, NULL);

  info->plaq = plaq;
  info->counter = counter;

  info->beta = g_beta;
  info->kappa = g_kappa;
  info->mu = g_mu / 2. / g_kappa;
  info->c2_rec = g_rgi_C1;
  info->time = t1.tv_sec;

  strcpy(info->package_version, PACKAGE_VERSION);

  info->mubar = g_mubar / 2. / g_kappa;
  info->epsilonbar = g_epsbar / 2. / g_kappa;

  strcpy(info->date, ctime(&t1.tv_sec));
  return(info);
}
