/* $Id$ */
#ifndef _INIT_CHI_SPINOR_FIELD_H
#define _INIT_CHI_SPINOR_FIELD_H

int init_chi_up_spinor_field(const int V, const int nr);
void free_chi_up_spinor_field();

int init_chi_dn_spinor_field(const int V, const int nr);
void free_chi_dn_spinor_field();

#endif
