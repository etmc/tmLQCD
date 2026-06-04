#ifndef PTBC_H
#define PTBC_H

#include <mpi.h>
#include <app.h>


/* for computing action with PTBC */
bool is_defect(PTBCDefect *def, int const ix, int const mu);
double get_ptbc_coeff(int const ix, int const mu);


void mpi_gather_base_rank();
void mpi_bcast_base_rank();
void mpi_base_rank_update_fini();

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

#endif