#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

#define LEAPFROG 1
#define SEXTON 2
#define EXTLEAPFROG 3
#define EXTSEXTON 4
#define IMPRLEAPFROG 5
#define MN2 6
#define MN2p 7

typedef void (*integratefk)(const double, const int, const int);

typedef struct{
  int type[10];
  int no_timescales;
  int n_int[10];
  double tau;
  double lambda[10];
  int mnls_per_ts[10][10];
  int no_mnls_per_ts[10];
  integratefk integrate[10];
} integrator;

extern integrator Integrator;


int init_integrator();
void integrate_2mn(const double tau, const int S, const int halfstep);
void integrate_2mnp(const double tau, const int S, const int halfstep);
void integrate_leap_frog(const double tau, const int S, const int halfstep);


#endif
