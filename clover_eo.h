#ifndef _CLOVER_EO_H
#define _CLOVER_EO_H

void Q_psi(int k, int l, double q_off);
void M_psi(int k, int l, double q_off);
void H_eo_psi(int ieo, int l, int k);
void deriv_Sb(int ieo, int l, int k);
void gamma5(int l, int k);
void boundary();

#endif
