/* $Id$ */
#ifndef _INIT_GAUGE_TMP_H
#define _INIT_GAUGE_TMP_H

extern su3 ** gauge_tmp;

int init_gauge_tmp(const int V);
void free_gauge_tmp();

#endif
