/***********************************************************************
 *
 * Copyright (C) 2026 JingJing Li
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#ifndef PTBC_H
#define PTBC_H

#include <mpi.h>
#include <app.h>


/* for computing action with PTBC.
 *
 * A link is identified by its start site and its direction. The start site is
 * passed as a LOCAL site ix plus a displacement rather than as the neighbour
 * index, because a neighbour may lie in the MPI halo, where g_coord is not
 * defined. get_ptbc_coeff() resolves the displacement in global coordinates
 * (with periodic wrap), so halo links get the correct coefficient.
 *
 * ix must satisfy ix < VOLUME; the displaced site need not.
 */
bool is_defect(PTBCDefect *def, int const coords[4], int const mu);
double get_ptbc_coeff(int const ix, int const disp[4], int const mu);

/* link (ix, dir) */
static inline double ptbc_coeff0(int const ix, int const dir) {
  int const disp[4] = {0, 0, 0, 0};
  return get_ptbc_coeff(ix, disp, dir);
}

/* link (ix + a*e_alpha, dir) */
static inline double ptbc_coeff1(int const ix, int const alpha, int const a, int const dir) {
  int disp[4] = {0, 0, 0, 0};
  disp[alpha] = a;
  return get_ptbc_coeff(ix, disp, dir);
}

/* link (ix + a*e_alpha + b*e_beta, dir); alpha == beta accumulates */
static inline double ptbc_coeff2(int const ix, int const alpha, int const a, int const beta,
                                 int const b, int const dir) {
  int disp[4] = {0, 0, 0, 0};
  disp[alpha] = a;
  disp[beta] += b;
  return get_ptbc_coeff(ix, disp, dir);
}

/* synchronising the ptbc topology whenever swaps happen */
void ptbc_sync();

/* for swapping rng */
typedef struct{
    int state_s[105];   // state of rlxs
    int state_d[105];   // state of rlxd
} SwapRNG;

void swap_rng(int const dest_inst);

/* Graph of instances */
typedef struct {
  int parent; // parent node
  int children[MAX_N_DEFECTS]; // children nodes
  int n_children; // number of children nodes  
} Node;

typedef struct {
    int root; // root node of the graph 
} Tree;

// setters
void set_tree_root(int const root);
void set_node_parent(int const node_id, int const parent_id);
void add_node_child(int const node_id, int const child_id);
void set_edge(int const parent_id, int const child_id);
void init_node(int const node_id);


// getters
int const get_tree_root();
int const get_node_parent(int const node_id);
int const* get_node_children(int const node_id);
int const get_node_n_children(int const node_id);


// initialiser
void init_ptbc_tree();
void print_ptbc_topo();


// utils
bool if_periodic(int inst_id);
void swap_rate(int const inst_id, int *rate);
const char* ptbc_timer_tag(int const inst_id);

/* Swapping functions (can make static if confirm later not used elsewhere) */
int swap_eo_tent(double const diff_up, double const diff_dn, int eo);
int init_eoswap_pbc(int *tent_order, int eo);
int swap_link(int const partner_inst, double const own_diff);

// gauge action with PTBC
double ptbc_swap_dh(int const alt_inst);

/* The even-odd swap */
void eo_swap(int *Rate, int const eo);


/* up or downstream swap */
void up_swap(int *Rate);
void down_swap(int *Rate);


/* io related */
void ptbc_chdir_instance(void);

#endif