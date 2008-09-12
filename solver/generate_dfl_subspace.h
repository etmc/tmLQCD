/* $Id$ */
#ifndef _GENERATE_DFL_SUBSPACE
#define _GENERATE_DFL_SUBSPACE

int init_dfl_subspace(const int);
int free_dfl_subspace();
int generate_dfl_subspace(const int Ns, const int N);
int generate_dfl_subspace_free(const int Ns, const int N);

extern spinor ** dfl_fields;
#endif
