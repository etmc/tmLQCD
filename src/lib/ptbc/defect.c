#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include <ptbc.h>
#include "global.h"
#include "ranlxs.h"
#include "ranlxd.h"

#define err(test, ...) err_impl(test, __func__, __FILE__, __LINE__, __VA_ARGS__)
static void err_impl(const bool test, const char* func, const char* file, const int line, const char* format, ...)
{
    if (test) {
        va_list args;
        char message[1024];
        va_start(args, format);
        vsnprintf(message, 1024, format, args);
        va_end(args);
        char location[1024];
        snprintf(location, 1024, "%s:%d %s", file, line, func);
        fatal_error(message, location);
    }
}

static int base_rank[MAX_N_INSTANCES];  // MPI rank offset of each instance
static MPI_Comm leader_comm;
static bool leader_comm_initialised = false;
static MPI_Request stats[2];

MPI_Datatype swap_info_type= MPI_DATATYPE_NULL; // MPI datatype for swapping info
static Tree tree;
static Node nodes[MAX_N_INSTANCES]; // each node corresponds to an instance

// find distance end-start corrected for peridic bc 
int const static dist(int start, int end, int const g_length) {
  int const d = end - start;
  int const dist = (d>=0) ? d : d + g_length; // assuming periodic bc

  return dist;
}

/**
 * @brief      check if a link lie in defect region, cuts at pos + 1/2 and pos + Ld + 1/2
 *             if either or end of the link is in the defect region, return true 
 *             include both case when link is in the defect and case link is crossing defect boundary
 *
 * @param      def  defect region
 * @param      ix   link start point
 * @param      mu   link direction
 */
bool is_defect(PTBCDefect *def, int const ix, int const mu) {
  if (ix >= VOLUME) return false;
  int const global_dim[4] = {T*g_nproc_t, LX*g_nproc_t, LY*g_nproc_x, LZ*g_nproc_z};
  int* coords = g_coord[ix];

  // check if mu direction is in the defect or crossing the cut
  // i.e. start >= pos[mu] and end <= pos[mu] + Ld[mu] + 1
  // dist >= 0 or dist + 1 <= Ld + 1
  int const dist_mu = dist(def->pos[mu], coords[mu], global_dim[mu]);
  if (dist_mu >= 0 && dist_mu <= def->Ld[mu]) {
    // find three non-mu directions
    int dim3[3];
    int count = 0;
    for (int d=0; d<4; d++) {
      if (d != mu) {
        dim3[count] = d;
        count++;
      }
    }
    // check if the other three directions are within the defect region
    // i.e. start at > pos end at <= pos + Ld
    // dist > 0 or dist + 1 <= Ld
    for (int d=0; d<3; d++) {
      int const dist_d = dist(def->pos[dim3[d]], coords[dim3[d]], global_dim[dim3[d]]);
      if (dist_d > 0 && dist_d < def->Ld[dim3[d]]){
        continue;
      }
      else{
        return false;
      }
    }
    return true;
  }
  else
    return false;  
}

/**
 * @brief      get multiplying factor of parallel tempering locally (assumne no overlapping defects!)
 *
 * @param      ix   starting point of link
 * @param      mu   direction of link
 */
double get_ptbc_coeff(int const ix, int const mu) {
  int const inst = app()->ptbc.instance_id;
  const PTBCInstance *instance = &(app()->ptbc.instances[inst]);

  // if instance is not active, return 1
  if (!instance->active) return 1.;
  
  // loop over defects
  for (int i=0; i<instance->n_coeffs; i++) {
    PTBCDefect *def = instance->defects[i];
    // apply coeff if within defect
    if (is_defect(def, ix, mu)) {
      return instance->coefficients[i];
    }
  }

  return 1.;
}


/**
 * @brief      Swap RNG state of with another instance
 *
 * @param      dest_inst  Destination instance ID
 */
void swap_rng(int const dest_inst) {
  SwapRNG info;
  
  if (dest_inst == app()->ptbc.instance_id) return;

  // find dest rank and src rank
  int local_rank;
  MPI_Comm_rank(app()->mpi.comm, &local_rank);
  int const dest_rank = base_rank[dest_inst] + local_rank;  
  int const src_rank = app()->mpi.world_rank;

  // get states
  rlxs_get(info.state_s);
  rlxd_get(info.state_d);

  // define MPI datatype for swapping info if not defined
  if (swap_info_type == MPI_DATATYPE_NULL) {
    int const blocklengths[2] = {105, 105};
    MPI_Aint const displacements[2] = {0, 105 * sizeof(int)};
    MPI_Datatype const types[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, blocklengths, displacements, types, &swap_info_type);
    MPI_Type_commit(&swap_info_type);
  }

  // swap rng and coeff info
  MPI_Status status;
  MPI_Sendrecv_replace(&info, 1, swap_info_type, dest_rank, 123, src_rank, 123, app()->mpi.world_comm, &status);
  if (status.MPI_ERROR != MPI_SUCCESS) {
    fprintf(stderr, "Error in MPI_Sendrecv_replace in swap_rng_coeff\n");
    MPI_Abort(app()->mpi.world_comm, status.MPI_ERROR);
  }

  // set state
  rlxs_reset(info.state_s);
  rlxd_reset(info.state_d);

  return;
}


/**
 * @brief Gather base rank offset from leader rank of all instances. Must use with the other two!
 * 
 */
void mpi_gather_base_rank() {
  int my_rank;
  MPI_Comm_rank(app()->mpi.comm, &my_rank);

  if (my_rank != 0) return; // only local rank 0 needs to send information, avoid redundancy

  // set communicator if not yet
  if (!leader_comm_initialised) {
    int leader_color = (my_rank == 0) ? 0 : MPI_UNDEFINED;
    
    // leader rank == instance_id
    MPI_Comm_split(app()->mpi.world_comm, leader_color, app()->ptbc.instance_id, &leader_comm);
    leader_comm_initialised = true;
  }

  // collect base rank offsets from leader ranks
  MPI_Iallgather(&(app()->mpi.world_rank), 1, MPI_INT, base_rank, 1, MPI_INT, leader_comm, stats);
}

/**
 * @brief Broadcast base rank offset to all local ranks from leader rank. Must use with the other two!
 * 
 */
void mpi_bcast_base_rank() {
  // check if gather complete
  MPI_Wait(stats, MPI_STATUS_IGNORE);

  // broadcast from leader rank (rank 0) to all other local ranks
  MPI_Ibcast(base_rank, app()->ptbc.n_instances, MPI_INT, 0, app()->mpi.comm, stats+1);
}

/**
 * @brief Finalise base rank update.
 * 
 */
void mpi_base_rank_update_fini(){
  // finalise
  MPI_Wait(stats+1, MPI_STATUS_IGNORE);
}

typedef struct {
  unsigned int n_coeffs;
  double coeffs[MAX_N_INSTANCES];
  unsigned int inst_id;
} InstanceInfo;

int compare_coeff(void const *ia, void const *ib) {
  unsigned int count1 = 0;
  unsigned int count2 = 0;

  const InstanceInfo * a = (const InstanceInfo *) ia;
  const InstanceInfo * b = (const InstanceInfo *) ib;
  
  for (unsigned int i=0; i<a->n_coeffs; i++) {
    if (a->coeffs[i] < b->coeffs[i]) count1++;
    else if (a->coeffs[i] > b->coeffs[i]) count2++;
  }
  if (count1 && count2==0) return 1;  // ib has bigger coeffs
  else if (count2 && count1==0) return -1;  // ia has bigger coeffs
  else if (count1==0 && count2==0) return 0; // two instances with exactly the same coeffs
  
  err((count1 && count2), "Error in compare_coeff: cannot sort instance coefficients in a strictly descending order!");
  return 0;
}

static int if_periodic(PTBCInstance const *instance) {
  for (int i=0; i<instance->n_coeffs; i++) {
    if (instance->coefficients[i] != 1) {
      return false;
    }
  }
  return true;
}


void set_tree_root(int const root) {tree.root = root;};
void set_node_parent(int const node_id, int const parent_id) {nodes[node_id].parent = parent_id;};
void add_node_child(int const node_id, int const child_id) {
    nodes[node_id].children[nodes[node_id].n_children] = child_id;
    nodes[node_id].n_children++;
};
void set_edge(int const parent_id, int const child_id) {
    set_node_parent(child_id, parent_id);
    add_node_child(parent_id, child_id);
}
void init_node(int const node_id) {nodes[node_id].parent=-1; nodes[node_id].n_children=0;};

int const get_tree_root() {return tree.root;};
int const get_node_parent(int const node_id) {return nodes[node_id].parent;};
int const* get_node_children(int const node_id) {return nodes[node_id].children;};
int const get_node_n_children(int const node_id) {return nodes[node_id].n_children;};


/**
 * @brief      Initialise PTBC instance connection graph.
 */
void init_ptbc_tree() {
  PTBCContext const *ptbc_ctx = &(app()->ptbc);
  InstanceInfo info[MAX_N_DEFECTS][MAX_N_INSTANCES]; // instance info n_tentacles x instance per tentacle
  int n_tentacles = 0;
  int tentacle_type[MAX_N_DEFECTS][MAX_N_DEFECTS]; // store defect id for each tentacle
  int tentacle_n_defect[MAX_N_DEFECTS]; // number of defects in each tentacle
  int tentacle_length[MAX_N_DEFECTS] = {0}; // count how many instances in each tentacle
  int n_periodic = 0; // number of periodic instances
  int periodic_id[MAX_N_INSTANCES]; // store periodic instance id, may be multiple periodic instances

  // initialise
  for (int i=0; i<ptbc_ctx->n_instances; i++) {
    init_node(i);
  }

  // log tentacles excluding periodic instances
  for (int id=0; id<ptbc_ctx->n_instances; id++) {
    PTBCInstance const * ins = &(ptbc_ctx->instances[id]);

    if (!if_periodic(ins)) {
      // check if defect already recorded, if not add to list, first exclude periodic instances
      int flag_new_tentacle = 1;
      for (int t=0; t<n_tentacles && flag_new_tentacle; t++) {
        int defect_match = 0;
        for (int nc=0; nc<ins->n_coeffs; nc++) {
          if (ins->defects[nc] - ptbc_ctx->defects == tentacle_type[t][nc]) {
            defect_match++;
          }
        }

        // if match, record instance info in the tentacle
        if (defect_match == ins->n_coeffs) {
          flag_new_tentacle = 0;
          info[t][tentacle_length[t]].n_coeffs = ins->n_coeffs;
          info[t][tentacle_length[t]].inst_id = id;
          for (int nc=0; nc<ins->n_coeffs; nc++) {
            info[t][tentacle_length[t]].coeffs[nc] = ins->coefficients[nc];
          }
          tentacle_length[t]++;
        }
      }

      // add new tentacle if no match
      if (flag_new_tentacle) {
        for (int nc=0; nc<ins->n_coeffs; nc++) {
          tentacle_type[n_tentacles][nc] = ins->defects[nc] - ptbc_ctx->defects;
        }
        tentacle_n_defect[n_tentacles] = ins->n_coeffs;

        // record first instance info in the new tentacle
        info[n_tentacles][0].n_coeffs = ins->n_coeffs;
        info[n_tentacles][0].inst_id = id;
        for (int nc=0; nc<ins->n_coeffs; nc++) {
          info[n_tentacles][0].coeffs[nc] = ins->coefficients[nc];
        }
        tentacle_length[n_tentacles]++;
        n_tentacles++;
      }
    }
    else {
      periodic_id[n_periodic] = id;
      n_periodic++;
    }
  }

  // order instances in each tentacle by coeff 1->0
  for (int t=0; t<n_tentacles; t++) {
    qsort(info[t], tentacle_length[t], sizeof(InstanceInfo), compare_coeff);
  }

  // log tree structure
  for (int t=0; t<n_tentacles; t++) {
    // connect instances in the same tentacle, ordered by coeff 1->0
    for (int n=0; n<tentacle_length[t]-1; n++) {
      set_edge(info[t][n].inst_id, info[t][n+1].inst_id);
    }
  }

  // set tree root
  if (n_periodic) {
    set_tree_root(periodic_id[0]);
  }
  else{
    err(n_periodic, "Error in init_ptbc_tree: no periodic instances found!");
  }

  // Sort periodic instances
  int connected_tentacles=0;
  for (int p=0; p<n_periodic; p++) {
    int periodic_edge_type[MAX_N_DEFECTS]; // store defect id for each periodic edge
    PTBCInstance const *period_instance = &(ptbc_ctx->instances[periodic_id[p]]);
    int n_periodic_edge=0;
    int matches[MAX_N_DEFECTS] = {0}; // count how many matches per tentacle
    
    // attach tentacles to periodic instance(s)
    for (int nc=0; nc<period_instance->n_coeffs; nc++) {
      int flag_match = 0;
      for (int t=0; t<n_tentacles; t++) {
        for (int nc_t=0; nc_t<tentacle_n_defect[t]; nc_t++) {
          if (period_instance->defects[nc] - ptbc_ctx->defects == tentacle_type[t][nc_t]) {
            matches[t]++;
            flag_match = 1;
          }
        }
      }
      if (!flag_match) {
        // must be a periodic-periodic edge
        periodic_edge_type[n_periodic_edge] = period_instance->defects[nc] - ptbc_ctx->defects;
        n_periodic_edge++;
      }
    }

    // connect periodic instance to tentacles
    for (int t=0; t<n_tentacles; t++) {
      if (matches[t] == tentacle_n_defect[t]) {
        // this tentacle is connected to the periodic root
        set_edge(periodic_id[p], info[t][0].inst_id);
        connected_tentacles++;
      }
      else if (matches[t] > 0) {
        err(matches[t], "Error in init_ptbc_tree: incomplete overlap of periodic instance with tentacles!");
      }
    }

    // connect between periodic instances
    // skip the last periodic as it has been checked already
    if (p==n_periodic-1) continue;
    int const cnt = n_periodic_edge;
    // loop over other periodic instances to find match
    for (int p2=p+1; p2<n_periodic; p2++) {
      for (int type=0; type<cnt; type++) {
        PTBCInstance const *instance_p2 = &(ptbc_ctx->instances[periodic_id[p2]]);
        for (int nc=0; nc<instance_p2->n_coeffs; nc++) {
          if (instance_p2->defects[nc] - ptbc_ctx->defects == periodic_edge_type[type]) {
            // connect two periodic instances
            set_edge(periodic_id[p], periodic_id[p2]);
            n_periodic_edge--;
          }
        }
      }
    }
    err(n_periodic_edge!=0, "Error in init_ptbc_tree: invalid connections between periodic instances!");
  }

  // err if tentacles connected to multiple periodic or loose tentacles exist
  err(connected_tentacles>n_tentacles, "Error in init_ptbc_tree: a tentacle is connected to multiple periodic instances!");
  err(connected_tentacles<n_tentacles, "Error in init_ptbc_tree: a tentacle is not connected to any periodic instance!");


  // check for loose connections 
  for (int id=0; id<ptbc_ctx->n_instances; id++) {
    // a lost node that is not connected to anything
    if (get_node_n_children(id) == 0 && get_node_parent(id) == -1) {
      err(1, "Error in init_ptbc_tree: disconnected instance(s) excist(s)!");
    }

    // a non-root periodic node that is not connected to another periodic node
    if (id != get_tree_root() && if_periodic(&(ptbc_ctx->instances[id])) && get_node_parent(id) == -1){
      err(1, "Error in init_ptbc_tree: disconnected PTBC chains!");
    }
  }

  if(g_proc_id == 0) printf("PTBC graph initialised: %d periodic instance(s), %d tentacle(s). \n", n_periodic, n_tentacles);
  return;
}

/**
 * @brief Print information of instances, including instance type upstream IDs and downstream IDs 
 * 
 */
void print_ptbc_topo() {
  // Iterate through the nodes
  const char* typeLabels[] = {
    "open",         // index 0
    "periodic",     // index 1
    "intermediate"  // index 2
  };

  for (int id=0; id<app()->ptbc.n_instances; id++) {
    int type;
    if (if_periodic(&(app()->ptbc.instances[id]))) {
      type = 1; // periodic
    }
    else if (get_node_n_children(id) == 0) {
      type = 0; // open
    }
    else {
      type = 2; // intermediate
    }
    
    printf("Instance ID: %d\tType: %s\t", id, typeLabels[type]);
    printf("Upstream ID: ");
    if (get_node_parent(id) >= 0) printf("%d", get_node_parent(id)); // if parent is -1 the node is the root
    printf("\tDownstream ID: ");
    for (int c=0; c<get_node_n_children(id); c++) {
      printf("%d, ", get_node_children(id)[c]);
    }
    printf("\n\n");
  }
}