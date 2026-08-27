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

#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <ptbc.h>
#include "global.h"
#include "ranlxs.h"
#include "ranlxd.h"
#include "read_input.h"


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

/* Global lattice extents. Fill it on first use. */
static int Ltot[4] = {0};

static void fill_Ltot(void) {
  if (Ltot[0] == 0) {
    Ltot[0] = T * g_nproc_t;
    Ltot[1] = LX * g_nproc_x;
    Ltot[2] = LY * g_nproc_y;
    Ltot[3] = LZ * g_nproc_z;
  }
}

static int base_rank[MAX_N_INSTANCES];  // MPI rank offset of each instance
static MPI_Comm leader_comm = MPI_COMM_NULL;  // communicator for leader ranks of each instance
static bool leader_comm_initialised = false;
static Tree tree;
static Node nodes[MAX_N_INSTANCES]; // each node corresponds to an instance

// PTBC RNG on global rank 0 distribute bundles of random number to leader ranks every step
// Worst-case consumption over one refill period  should is bounded by 
// MAX_N_INSTANCES + MAX_N_DEFECTS - 2,therefore 2 * MAX_N_INSTANCES should be safe
// There are also an out-of-bound checks at runtime see next_random() for details
#define N_RANDOM_PER_INSTANCE (2 * MAX_N_INSTANCES)

static double random_numbers[N_RANDOM_PER_INSTANCE]; // every leader rank rank gets a bundle of random numbers from rank 0
static double *random_numbers_ptr = NULL; // pointer to the next random number in the bundle

/**
 * @brief      Consume the next random number from this leader rank's bundle.
 *             Only leader ranks hold a bundle; calling this elsewhere is a bug.
 *
 * @return     uniform double in [0,1)
 */
double static next_random(void) {
  err(random_numbers_ptr == NULL, "Error in next_random: PTBC random bundle not initialised!");
  err(random_numbers_ptr >= random_numbers + N_RANDOM_PER_INSTANCE,
      "Error in next_random: PTBC random bundle of %d exhausted on instance %d!",
      N_RANDOM_PER_INSTANCE, app()->ptbc.instance_id);

  return *(random_numbers_ptr++);
}

// find distance end-start corrected for peridic bc 
int const static dist(int start, int end, int const g_length) {
  int const d = end - start;
  int const dist = (d>=0) ? d : d + g_length; // assuming periodic bc

  return dist;
}

/**
 * @brief      check if a link lie in defect region
 *             i.e. link is internal to the defect or link is crossing the defect surface
 *
 * @param      def  defect region
 * @param      coords   link start point
 * @param      mu   link direction
 */
bool is_defect(PTBCDefect *def, int const coords[4], int const mu) {

  for (int i = 0; i < 4; i++) {
    int const d = dist(def->pos[i], coords[i], Ltot[i]);
    if (i == mu) {
      if (d > def->Ld[i]) return false;             // link rule: dist in [0, Ld], inclusive
    } else {
      if (d == 0 || d > def->Ld[i]) return false;    // point rule: dist in [1, Ld], strict interior
    }
  }
  return true;
}

/**
 * @brief      get multiplying factor of parallel tempering locally (assumne no overlapping defects!)
 *
 * @param      ix   starting point of link
 * @param      mu   direction of link
 * @param      disp displacement from ix (could extend to halo region)
 */
double get_ptbc_coeff(int const ix, int const disp[4], int const mu) {
  int const inst = app()->ptbc.instance_id;
  const PTBCInstance *instance = &(app()->ptbc.instances[inst]);

  // if instance is not active, return 1
  if (!instance->active) return 1.;
  
  // the base site must be local: g_coord is only defined for ix < VOLUME.
  err(ix < 0 || ix >= VOLUME, "Error in get_ptbc_coeff: base site is not local!");
  int coords[4];
  for (int d=0; d<4; d++) coords[d] = (g_coord[ix][d] + disp[d] + Ltot[d]) % Ltot[d];
  // loop over defects
  for (int i=0; i<instance->n_coeffs; i++) {
    PTBCDefect *def = instance->defects[i];
    // apply coeff if within defect
    if (is_defect(def, coords, mu)) {
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
  // SwapRNG info;
  int state[210];
  if (dest_inst == app()->ptbc.instance_id) return;

  // find dest rank and src rank
  int local_rank;
  MPI_Comm_rank(app()->mpi.comm, &local_rank);
  int const dest_rank = base_rank[dest_inst] + local_rank;
  int const src_rank = app()->mpi.world_rank;

  // get states
  rlxs_get(state);
  rlxd_get(state+105);

  // printf("to %d rng = %d\n", dest_rank, state[100]);

  // swap rng and coeff info
  MPI_Status status;
  MPI_Sendrecv_replace(state, 210, MPI_INT, dest_rank, src_rank, dest_rank, dest_rank, app()->mpi.world_comm, &status);

  // set state
  rlxs_reset(state);
  rlxd_reset(state+105);

  return;
}


static void set_leader_comm() {
  int my_rank;
  MPI_Comm_rank(app()->mpi.comm, &my_rank);
  // set communicator if not yet
  if (!leader_comm_initialised) {
    int leader_color = (my_rank == 0) ? 0 : MPI_UNDEFINED;
    
    // leader rank == instance_id
    MPI_Comm_split(app()->mpi.world_comm, leader_color, app()->ptbc.instance_id, &leader_comm);
    leader_comm_initialised = true;
  }
}

/**
 * @brief Gather base rank offset from leader rank of all instances.
 * 
 */
void static mpi_gather_base_rank() {
  int my_rank;
  MPI_Comm_rank(app()->mpi.comm, &my_rank);

  set_leader_comm();

  // collect base rank offsets from leader ranks
  if (my_rank==0) MPI_Allgather(&(app()->mpi.world_rank), 1, MPI_INT, base_rank, 1, MPI_INT, leader_comm);

  // share with rest of local ranks
  MPI_Bcast(base_rank, app()->ptbc.n_instances, MPI_INT, 0, app()->mpi.comm);

  return;
}

typedef struct {
  unsigned int n_coeffs;
  double coeffs[MAX_N_INSTANCES];
  unsigned int inst_id;
} InstanceInfo;

int static compare_coeff(void const *ia, void const *ib) {
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


/**
 * @brief check if inst_id is a periodic instance
 * 
 * @param inst_id 
 * @return true 
 * @return false 
 */
bool if_periodic(int inst_id) {
  PTBCInstance const *instance = app()->ptbc.instances + inst_id;
  for (int i=0; i<instance->n_coeffs; i++) {
    if (instance->coefficients[i] != 1) {
      return false;
    }
  }
  return true;
}

/**
 * @brief format a small "instN" tag for tm_stopwatch_pop's prefix argument, so PTBC
 *        timer lines show which instance produced them. 
 *
 * @param inst_id captured instance id to tag the line with
 */
const char* ptbc_timer_tag(int const inst_id) {
  static char buf[32];
  snprintf(buf, sizeof(buf), "inst%d", inst_id);
  return buf;
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


#define PTBC_RNG_STATE_SIZE 105
#define PTBC_RNG_STATE_FILE "../ptbc_rng_state"

/**
 * @brief      Write the PTBC swap-decision RNG state to the result directory so that a
 *             restart continues the stream instead of replaying it. Global rank 0 owns
 *             the generator and is the only rank that writes
 */
void write_ptbc_rng_state(void) {
  if (!app()->ptbc.active || app()->mpi.world_rank != 0) return;

  FILE *fp = fopen(PTBC_RNG_STATE_FILE, "w");
  err(fp == NULL, "Error in write_ptbc_rng_state: cannot open %s for writing!", PTBC_RNG_STATE_FILE);

  for (int i = 0; i < PTBC_RNG_STATE_SIZE; i++) {
    fprintf(fp, "%d\n", app()->ptbc.rng_state[i]);
  }
  err(fclose(fp) != 0, "Error in write_ptbc_rng_state: failed to write %s!", PTBC_RNG_STATE_FILE);
}

/**
 * @brief      Restore the PTBC RNG state written by a previous run. Global rank 0 only. 
 *             Abscence of PTBC_RNG_STATE_FILE file means this is a fresh start
 *
 * @return     true if the state was restored, false if there is no state file
 */
bool static read_ptbc_rng_state(void) {
  FILE *fp = fopen(PTBC_RNG_STATE_FILE, "r");
  if (fp == NULL) return false;  // fresh start: seed from ptbc.seed instead

  int state[PTBC_RNG_STATE_SIZE];
  for (int i = 0; i < PTBC_RNG_STATE_SIZE; i++) {
    if (fscanf(fp, "%d", &state[i]) != 1) {
      fclose(fp);
      err(1, "Error in read_ptbc_rng_state: %s is truncated at entry %d!", PTBC_RNG_STATE_FILE, i);
    }
  }
  fclose(fp);

  // rlxd_reset() exits on a bad state, so check what we can first to give a usable message
  err(state[0] != PTBC_RNG_STATE_SIZE,
      "Error in read_ptbc_rng_state: %s holds a state of size %d, expected %d!",
      PTBC_RNG_STATE_FILE, state[0], PTBC_RNG_STATE_SIZE);

  rlxd_reset(state);
  memcpy(appm()->ptbc.rng_state, state, sizeof(state));
  return true;
}

static bool rng_initialised = false;

/**
 * @brief      (Re)fill each leader rank's bundle of swap-decision random numbers.
 *
 *             The PTBC generator is owned by global rank 0 and by nobody else: it is
 *             seeded from ptbc.seed (or restored from PTBC_RNG_STATE_FILE on a restart)
 *             on the first call, and carried forward in ptbc.rng_state. 
 *
 *             Only leader ranks receive a bundle.
 */
void static fill_random_numbers(void) {
  int local_rank;
  MPI_Comm_rank(app()->mpi.comm, &local_rank);

  // find which instance global rank 0 currently resides in
  int rng_inst_id = -1;
  for (int i=0; i<app()->ptbc.n_instances; i++) {
    if (base_rank[i] == 0) {
      rng_inst_id = i;
      break;
    }
  }
  err(rng_inst_id < 0, "Error in fill_random_numbers: cannot find the instance with world rank 0!");

  if (app()->mpi.world_rank == 0) {
    err(rng_inst_id != app()->ptbc.instance_id,
        "Error in fill_random_numbers: mismatch RNG instance and base rank!");

    int tmp[PTBC_RNG_STATE_SIZE];

    // Whether rng is parked in tmp and need to be restored
    bool const parked = rng_initialised;

    // if not initialised, either restore from file or init from seed
    if (!rng_initialised) {
      if (!read_ptbc_rng_state()) {
        err(app()->ptbc.seed <= 0, "Error in fill_random_numbers: PTBC seed is invalid");
        rlxd_init(rlxd_level, app()->ptbc.seed);
      }
      rng_initialised = true;
    }
    else {
      rlxd_get(tmp);                       // park the physics stream
      rlxd_reset(appm()->ptbc.rng_state);  // resume the PTBC stream
    }

    double numbers[N_RANDOM_PER_INSTANCE * app()->ptbc.n_instances];
    ranlxd(numbers, N_RANDOM_PER_INSTANCE * app()->ptbc.n_instances); // uniform doubles in [0,1)
    rlxd_get(appm()->ptbc.rng_state);      // carry the PTBC stream forward

    if (parked) rlxd_reset(tmp);           // restore the physics stream

    MPI_Scatter(numbers, N_RANDOM_PER_INSTANCE, MPI_DOUBLE, random_numbers, N_RANDOM_PER_INSTANCE, MPI_DOUBLE, rng_inst_id, leader_comm);
  }
  else if (local_rank == 0) {
    MPI_Scatter(NULL, N_RANDOM_PER_INSTANCE, MPI_DOUBLE, random_numbers, N_RANDOM_PER_INSTANCE, MPI_DOUBLE, rng_inst_id, leader_comm);
  }

  random_numbers_ptr = random_numbers;  // reset pointer to the start of the array
  return;
}

/**
 * @brief      Initialise PTBC instance connection graph.
 */
void init_ptbc_tree() {
  // modify invalid positions
  fill_Ltot();
  for (int n_def=0; n_def <app()->ptbc.n_defects; n_def++) {
    int *positions = appm()->ptbc.defects[n_def].pos;
    for (int d=0; d<4; d++) {
      while (positions[d] < 0) {
        positions[d] += Ltot[d];
      }
      positions[d] = positions[d] % Ltot[d];
      err(app()->ptbc.defects[n_def].Ld[d] < 0, "Negative defect extent!");
    }
  }

  PTBCContext const *ptbc_ctx = &(app()->ptbc);
  InstanceInfo info[MAX_N_DEFECTS][MAX_N_INSTANCES]; // instance info n_tentacles x instance per tentacle
  int n_tentacles = 0;
  int tentacle_type[MAX_N_DEFECTS][MAX_N_DEFECTS]; // store defect id for each tentacle
  int tentacle_n_defect[MAX_N_DEFECTS]; // number of defects in each tentacle
  int tentacle_length[MAX_N_DEFECTS] = {0}; // count how many instances in each tentacle
  int n_periodic = 0; // number of periodic instances
  int periodic_id[MAX_N_INSTANCES]; // store periodic instance id, may be multiple periodic instances

  // check if defects are valid / no overlap
  for (int i=0; i<ptbc_ctx->n_defects-1; i++) {
    PTBCDefect const* def_ref = &(ptbc_ctx->defects[i]);

    for (int j=i+1; j<ptbc_ctx->n_defects; j++) {
      PTBCDefect const* def = &(ptbc_ctx->defects[j]);

      // two per-axis circular intervals [pos, pos+Ld] overlap iff either
      // defect's start falls within the other's forward span; the 4D boxes
      // overlap only if that holds on every axis simultaneously
      bool def_overlap = true;
      for (int d=0; d<4; d++) {
        bool const axis_overlap =
            dist(def_ref->pos[d], def->pos[d], Ltot[d]) <= def_ref->Ld[d] ||
            dist(def->pos[d], def_ref->pos[d], Ltot[d]) <= def->Ld[d];
        if (!axis_overlap) { def_overlap = false; break; }
      }
      err(def_overlap, "Some defects overlap!");
    }
  }

  // initialise
  for (int i=0; i<ptbc_ctx->n_instances; i++) {
    init_node(i);
  }

  // log tentacles excluding periodic instances
  for (int id=0; id<ptbc_ctx->n_instances; id++) {
    PTBCInstance const * ins = &(ptbc_ctx->instances[id]);

    if (!if_periodic(id)) {
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

  // err if periodic node is not connected to any tentacles
  for (int i=0; i<n_periodic; i++) {
    err(get_node_n_children(periodic_id[i]) == 0, "Error in init_ptbc_tree: certain periodic instance has no tentacles!");
  }

  // check for loose connections 
  for (int id=0; id<ptbc_ctx->n_instances; id++) {
    // a lost node that is not connected to anything
    if (get_node_n_children(id) == 0 && get_node_parent(id) == -1) {
      err(1, "Error in init_ptbc_tree: disconnected instance(s) excist(s)!");
    }

    // a non-root periodic node that is not connected to another periodic node
    if (id != get_tree_root() && if_periodic(id)){
      err(1, "Error in init_ptbc_tree: disconnected PTBC chains!");
    }
  }

  if(g_proc_id == 0) printf("PTBC graph initialised: %d periodic instance(s), %d tentacle(s). \n", n_periodic, n_tentacles);
  

  // set leader communicator for the first time
  set_leader_comm();  
  int size;
  MPI_Comm_size(app()->mpi.world_comm, &size);
  for (int i=0; i<ptbc_ctx->n_instances; i++) base_rank[i] = i * size / ptbc_ctx->n_instances;


  // initialise PTBC rng state if global rank 0
  fill_random_numbers();
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
    if (if_periodic(id)) {
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


/**
 * @brief Function to shuffle the array using the Fisher-Yates algorithm
 * 
 */
void static shuffle(int *arr, int n) {

  for (int i = n - 1; i > 0; i--) {
      // Map double in [0,1) to an index between 0 and i inclusive
      int j = (int)(next_random() * (i + 1));

      // Swap arr[i] with the element at the random index
      int temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
  }

  return;
}

void static swap(int const inst_id) {
  swap_rng(inst_id);  // swap rng
  appm()->ptbc.instance_id = inst_id; // swap instance id
}

/**
 * @brief swap an array between instances
 * 
 * @param inst_id swapping partner
 * @param arr
 * @param size  size of the array 
 * @return void 
 */
void static swap_arr(int const inst_id, int *arr, int const size) {
  int local_rank;
  MPI_Comm_rank(app()->mpi.comm, &local_rank);

  if (local_rank == 0) {
    MPI_Sendrecv_replace(arr, size, MPI_INT, base_rank[inst_id], 0, base_rank[inst_id], 0, app()->mpi.world_comm, MPI_STATUS_IGNORE);
  }

  MPI_Bcast(arr, size, MPI_INT, 0, app()->mpi.comm); // arr the same within an instance

  return;
}
/**
 * @brief check for swap if for even connections on the tentacles
 * 
 * @param d_up action difference with upstream
 * @param d_dn action difference with downstream
 * @param eo  even or odd swaps, 0 for even 1 for odd
 * @return int 1 if need to further swap with pbc 0 if done
 */

int static swap_eo_tent(double const d_up, double const d_dn, int eo) {
  double diff_up[MAX_N_INSTANCES], diff_dn[MAX_N_INSTANCES];
  int local_info[2], local_rank;
  MPI_Comm_rank(app()->mpi.comm, &local_rank);
  int local_inst = app()->ptbc.instance_id;

  if (local_rank == 0) {
    MPI_Gather(&d_up, 1, MPI_DOUBLE, diff_up, 1, MPI_DOUBLE, 0, leader_comm);
    MPI_Gather(&d_dn, 1, MPI_DOUBLE, diff_dn, 1, MPI_DOUBLE, 0, leader_comm);
  }

  // rank 0 of instance 0 is the coordinator for swap decisions
  if (local_inst == 0 && local_rank == 0) {
    // count number of tentacles and find open bc positions
    int pos[MAX_N_DEFECTS];
    int n_tent = 0;
    for (int i=0; i<app()->ptbc.n_instances; i++) {
      if (nodes[i].n_children == 0) {
        pos[n_tent] = i;
        n_tent ++;
      }
    }

    // make swapping decisions
    int swap_info[MAX_N_INSTANCES * 2]={0}; // store decision and partner
    // check for all tentacles
    for (int t=0; t<n_tent; t++) {
      int upstream = get_node_parent(pos[t]);
      if (eo) { // if odd start from link 1
        pos[t] = upstream;
        upstream = get_node_parent(pos[t]);
      }
      // links connected to periodic nodes are not swapped
      while(!if_periodic(pos[t]) && !if_periodic(upstream) && local_rank == 0) {
        double const total_diff = diff_up[pos[t]] + diff_dn[upstream];
        double expmdh = exp(-total_diff);

        double const u = next_random();
        int const accept = expmdh > u;
        printf("[swap_eo_link] inst %d <-> %d: own_diff=%.6e partner_diff=%.6e total_diff=%.6e expmdh=%.6e random number=%.6e accept=%d\n",
          pos[t], upstream, diff_up[pos[t]], diff_dn[upstream], total_diff, expmdh, u, accept);

        
        swap_info[pos[t] * 2] = accept;
        swap_info[pos[t] * 2 + 1] = upstream;
        swap_info[upstream * 2] = accept;
        swap_info[upstream * 2 + 1] = pos[t];

        pos[t] = get_node_parent(upstream);
        upstream = get_node_parent(pos[t]);

      }
    }

    // notify swap decisions to all instances
    MPI_Scatter(swap_info, 2, MPI_INT, local_info, 2, MPI_INT, 0, leader_comm);
  } else if (local_rank == 0) {
    // receive swap decisions from coordinator
    MPI_Scatter(NULL, 2, MPI_INT, local_info, 2, MPI_INT, 0, leader_comm);
  }
  
  // bcast to rest of the local ranks in the instance
  MPI_Bcast(local_info, 2, MPI_INT, 0, app()->mpi.comm);

  // execute swap if accept
  if (local_info[0]) {
    swap(local_info[1]);  // swap instance id
  }

  return local_info[0];  // return 1 if swap accepted, 0 if not
}

/**
 * @brief synchronise the base ranks of all instances after swap(s)
 * 
 */
void static ptbc_sync() {
  int local_rank;
  MPI_Comm_rank(app()->mpi.comm, & local_rank);
  // reset leader comm to be reinitialised for new instance.
  //Only local_rank==0 actually holds a real leader_comm handle to free.
  leader_comm_initialised = false;
  if (local_rank == 0) {
    MPI_Comm_free(&leader_comm);
  }

  mpi_gather_base_rank();   // update offsets
}

/**
 * @brief Get the number of periodic instance
 * 
 * @return int number of periodic instances
 */
int get_n_periodic() {
  int n_periodic = 0;
  for (int i=0; i<app()->ptbc.n_instances; i++) {
    if (if_periodic(i)) {
      n_periodic++;
    }
  }
  return n_periodic;
}

/**
 * @brief Get the periodic nodes
 * 
 * @param periodic_id the array of periodic IDs
 * @return int number of periodic instances
 */
int get_periodic(int *periodic_id) {
  int n_periodic = 0;
  for (int i=0; i<app()->ptbc.n_instances; i++) {
    if (if_periodic(i)) {
      periodic_id[n_periodic] = i;
      n_periodic++;
    }
  }
  return n_periodic;
}


/**
 * @brief initialise periodic-instance swap order; for a periodic instance, returns the
 *        number of tentacles swapping with it this pass; for a tentacle instance, returns
 *        the round index (matching eo_swap's loop) at which it swaps with its periodic
 *        neighbour, or -1 if it is not swapping this pass
 *
 * @param tent_order  array of tentacle indices that are involved in the even/odd swap (periodic instances only)
 * @param eo even or odd swaps, 0 for even 1 for odd
 * @return int  periodic instance: number of tentacles involved; tentacle instance: assigned round index, or -1
 */
int static init_eoswap_pbc(int *tent_order, int eo){
  int local_inst = app()->ptbc.instance_id;
  int local_rank;
  MPI_Comm_rank(app()->mpi.comm, &local_rank);
  int my_round = -1; // round index at which a tentacle swaps with its periodic neighbour, -1 if not swapping

  // if periodic instance, may need swap with multiple tentacles
  if (if_periodic(local_inst) && local_rank == 0) {
    // find distances from open bc, or return 1 for dummy links between periodic BCs
    int *dist = (int *)malloc(sizeof(int) * get_node_n_children(local_inst));
    for (int t=0; t<get_node_n_children(local_inst); t++) {
      int current_node = get_node_children(local_inst)[t];
      int tent_length = 1;
      while (get_node_n_children(current_node) != 0 && !if_periodic(current_node)) {
        current_node = get_node_children(current_node)[0];
        tent_length++;
      }
      dist[t] = tent_length;
    }

    // determine tentalces that are in the even swap
    int *swapping_tents = (int *)malloc(sizeof(int) * get_node_n_children(local_inst));
    int n_swapping_tents = 0;
    for (int t=0; t<get_node_n_children(local_inst); t++) {
      // find tentacles involved in even swap, exclude dummy
      if ((dist[t]+eo)%2 == 1 && !if_periodic(get_node_children(local_inst)[t])) {
        swapping_tents[n_swapping_tents] = t;
        n_swapping_tents++;
      }
    }


    // shuffle index, randomise which direction is swapped first
    shuffle(swapping_tents, n_swapping_tents);

    // return random order tentacle idices 
    for (int i=0; i<n_swapping_tents; i++) {
      tent_order[i] = swapping_tents[i];
    }
    
    // notify other ranks
    int * decision = (int *)malloc(sizeof(int) * app()->ptbc.n_instances);
    for (int i=0; i<app()->ptbc.n_instances; i++)  decision[i] = -1;
    for (int i=0; i<n_swapping_tents; i++)  
      decision[get_node_children(local_inst)[swapping_tents[i]]] = i;

    int periodic_id[MAX_N_DEFECTS];
    int const n_periodic = get_periodic(periodic_id);
    
    MPI_Request* request = (MPI_Request*)malloc(n_periodic * sizeof(MPI_Request));
    int *decisions = (int *) malloc(sizeof(int) * n_periodic); // dummy receives
    for (int i = 0; i < n_periodic; i++) {
      if (periodic_id[i] == local_inst) MPI_Iscatter(decision, 1, MPI_INT, &decisions[i], 1, MPI_INT, local_inst, leader_comm, &request[i]);
      else                              MPI_Iscatter(NULL, 1, MPI_INT, &decisions[i], 1, MPI_INT, periodic_id[i], leader_comm, &request[i]);
    }

    
    MPI_Bcast(&n_swapping_tents, 1, MPI_INT, local_rank, app()->mpi.comm);
    MPI_Bcast(tent_order, n_swapping_tents, MPI_INT, local_rank, app()->mpi.comm);

    MPI_Waitall(n_periodic, request, MPI_STATUSES_IGNORE);

    free(swapping_tents);
    free(dist);
    free(request);
    free(decision);
    free(decisions);

    return n_swapping_tents; // return number of swaps for the periodic instance
  }
  else if (local_rank == 0 && !if_periodic(local_inst)) {
    // count number of periodic instances
    int periodic_id[MAX_N_DEFECTS];
    int const n_periodic = get_periodic(periodic_id);

    MPI_Request* request = (MPI_Request*)malloc(n_periodic * sizeof(MPI_Request));
    // receive from all periodic ranks
    int decisions[MAX_N_DEFECTS] = {0};
    for (int i=0; i<n_periodic; i++) {
      MPI_Iscatter(NULL, 1, MPI_INT, &decisions[i], 1, MPI_INT, periodic_id[i], leader_comm, request+i);
    }

    MPI_Waitall(n_periodic, request, MPI_STATUSES_IGNORE);
    free(request);

    // if selected, record my round index
    for (int i=0; i<n_periodic; i++) {
      if (decisions[i] != -1) my_round = decisions[i];
    }

    MPI_Bcast(&my_round, 1, MPI_INT, 0, app()->mpi.comm);
    return my_round;
  }
  else if (if_periodic(local_inst) && local_rank!=0) {
    int n_swapping_tents;
    MPI_Bcast(&n_swapping_tents, 1, MPI_INT, 0, app()->mpi.comm);
    MPI_Bcast(tent_order, n_swapping_tents, MPI_INT, 0, app()->mpi.comm);
    return n_swapping_tents; // return number of swaps for the periodic instance
  }
  else {
    MPI_Bcast(&my_round, 1, MPI_INT, 0, app()->mpi.comm);
    return my_round;
  }
}

/**
 * @brief Try swapping between periodic instance and another neighbour. Decision made by pbc
 *
 * @param partner_inst
 * @param own_diff
 * @return int 1 accept; 0 reject
 */
int static swap_link(int const partner_inst, double const own_diff) {
  int local_rank;
  MPI_Comm_rank(app()->mpi.comm, &local_rank);
  int const local_inst = app()->ptbc.instance_id;

  double partner_diff;
  int accept;
  if (local_rank==0){ 
    MPI_Sendrecv(&own_diff, 1, MPI_DOUBLE, partner_inst, 0,
                  &partner_diff, 1, MPI_DOUBLE, partner_inst, 0,
                  leader_comm, MPI_STATUS_IGNORE);

    // the upstream decides
    if (get_node_parent(partner_inst)==local_inst && local_rank==0) {
      double const expmdh = exp(-(own_diff + partner_diff));
      double const u = next_random();
      accept = expmdh > u;
      printf("[swap_link] inst %d <-> %d: own_diff=%.6e partner_diff=%.6e total_diff=%.6e expmdh=%.6e random number=%.6e accept=%d\n",
            local_inst, partner_inst, own_diff, partner_diff, own_diff + partner_diff, expmdh, u, accept);
      MPI_Send(&accept, 1, MPI_INT, partner_inst, partner_inst, leader_comm);
    } 
    else if (get_node_parent(local_inst)==partner_inst && local_rank==0) {
      MPI_Recv(&accept, 1, MPI_INT, partner_inst, local_inst, leader_comm, MPI_STATUS_IGNORE);
    } 
    else {
      err(get_node_parent(local_inst)!=partner_inst && get_node_parent(partner_inst)!=local_inst, "Error in swap_link: the two instances are not neighbours!");
    }
  }

  MPI_Bcast(&accept, 1, MPI_INT, 0, app()->mpi.comm);
  if (accept) swap(partner_inst);
  return accept; // return 1 if swap accepted, 0 if not
}

/**
 * @brief Perform dummy swap
 * @param Rate
 */
void static swap_dummy(int* Rate) {
  int local_rank;
  MPI_Comm_rank(app()->mpi.comm, &local_rank);
  int local_inst = app()->ptbc.instance_id;

  // find the pbc
  int periodic_id[MAX_N_DEFECTS] = {0};
  int const n_periodic = get_periodic(periodic_id);
  
  int *link_id = (int *) malloc(sizeof(int) * (n_periodic - 1)); // dummy links

  if (local_inst == 0 && local_rank == 0) {
    // exclude root, left with dummy links
    int cnt = 0;
    int const root = get_tree_root();
    for (int i=0; i<n_periodic; i++) {
      if (periodic_id[i] != root) {
        link_id[cnt] = periodic_id[i]; 
        cnt++;
      }
    }

    // randomise link order
    shuffle(link_id, n_periodic - 1);
  }

  // broadcast shuffle result
  if (local_rank == 0) MPI_Bcast(link_id, n_periodic - 1, MPI_INT, 0, leader_comm);
  MPI_Bcast(link_id, n_periodic-1, MPI_INT, 0, app()->mpi.comm);

  // loop over all dummy links
  for (int i=0; i<n_periodic-1; i++) {
    int accepted = 0;
    int partner = -1;
    if (if_periodic(local_inst)) {
      if (local_inst == link_id[i]) {
        partner = get_node_parent(local_inst);
        double const diff = ptbc_swap_dh(partner);
        accepted = swap_link(partner, diff);
      }
      else if (local_inst == get_node_parent(link_id[i])) {
        partner = link_id[i];
        double const diff = ptbc_swap_dh(partner);
        accepted = swap_link(partner, diff);
      }
    }
    if (accepted) {
      // addressed via partner (not local_inst/own id), same convention as swap_rng(): see the
      // note in eo_swap's pbc loop.
      swap_arr(partner, Rate, 1);
      local_inst = app()->ptbc.instance_id;
    }
    ptbc_sync();
  }

  free(link_id);
  return;
}

/**
 * @brief Perform even or odd swap
 * 
 * @param Rate the swap rate follows instance id
 * @param eo   even 0, odd 1
 */
void static eo_swap(int *Rate, int eo){
  // swap tentacle
  int const inst_id = app()->ptbc.instance_id;
  double const diff_up = get_node_parent(inst_id) != -1? ptbc_swap_dh(get_node_parent(inst_id)): 0.;
  double const diff_dn = get_node_n_children(inst_id) > 0 ? ptbc_swap_dh(get_node_children(inst_id)[0]) : 0.;
  
  int const swap_accepted = swap_eo_tent(diff_up, diff_dn, eo); // even-odd swap
  ptbc_sync();
  if (swap_accepted){ 
    swap_arr(inst_id, Rate, 1);  // if swap happened, swap acceptance rate
  }
  int inst_id_new = app()->ptbc.instance_id;

  // initialise pbc swaps
  int tent_order[MAX_N_DEFECTS];
  int n_swaps = init_eoswap_pbc(tent_order, eo);

  // find number of swaps because all ranks need to synchronise after any swap
  int max_swaps;
  MPI_Allreduce(&n_swaps, &max_swaps, 1, MPI_INT, MPI_MAX, app()->mpi.world_comm);

  // swaps with periodic
  // periodic instances swap on every round i < n_swaps (n_swaps = its count of swapping tentacles,
  // each round picking the next tentacle from tent_order); a tentacle instance swaps when it is
  // my_round. All these information swaps when a swap happens
  for (int i = 0; i < max_swaps; i++) {
    int const inst_is_periodic = if_periodic(inst_id_new);
    int const my_turn = inst_is_periodic ? (i < n_swaps) : (i == n_swaps);
    int accepted = 0;
    int partner_inst = -1;
    if (my_turn) {
      partner_inst = inst_is_periodic ? get_node_children(inst_id_new)[tent_order[i]] : get_node_parent(inst_id_new);
      double const own_diff = ptbc_swap_dh(partner_inst);
      accepted = swap_link(partner_inst, own_diff);
    }
    if (accepted) {
      // properties follow the instance id
      swap_arr(partner_inst, Rate, 1);
      swap_arr(partner_inst, &n_swaps, 1);
      swap_arr(partner_inst, tent_order, MAX_N_DEFECTS);
      inst_id_new = app()->ptbc.instance_id;
    }
    ptbc_sync();
  }

  if (eo==0)  swap_dummy(Rate); // dummies are length 1, even swap

  return;
}


/**
 * @brief Perform even then odd swaps
 * 
 * @param Rate 
 * @return * void 
 */
void even_odd_swap(int *Rate) {
  eo_swap(Rate, 0); // even
  eo_swap(Rate, 1); // odd
  fill_random_numbers(); // refill random numbers for next round
  return;
}

/**
 * @brief Perform odd then even swaps
 * 
 * @param Rate 
 * @return * void 
 */
void odd_even_swap(int *Rate) {
  eo_swap(Rate, 1); // odd
  eo_swap(Rate, 0); // even
  fill_random_numbers(); // refill random numbers for next round
  return;
}

/**
 * @brief swap upstream (open to periodic) on the tentacles
 * 
 * @param Rate  current acceptance rate
 * @param ud    upstream 0 downstream 1
 * 
 */
void static swap_updown_tent(int *Rate, int const ud){
  int local_inst = app()->ptbc.instance_id;
  int local_rank;
  MPI_Comm_rank(app()->mpi.comm, &local_rank);

  // count number of syncs
  int n_max_links;
  // count number of tentacles and find open bc positions
  int pos[MAX_N_DEFECTS];
  int n_tent = 0;
  for (int i=0; i<app()->ptbc.n_instances; i++) {
    if (nodes[i].n_children == 0) {
      pos[n_tent] = i;
      n_tent ++;
    }
  }

  // find lengths of tentacles
  int lengths[MAX_N_DEFECTS] = {0};
  for (int t=0; t<n_tent; t++) {
    int upstream = get_node_parent(pos[t]);
    while(!if_periodic(upstream)) {
      lengths[t]++;
      pos[t] = upstream;
      upstream = get_node_parent(pos[t]);
    }
  }

  // find max number of links
  n_max_links = lengths[0];
  for (int i=1; i<n_tent; i++) {
    if (lengths[i] > n_max_links) n_max_links = lengths[i];
  }
    

  // count own distance to open 
  int d = 0;
  if (!if_periodic(local_inst)) {
    int id = local_inst;
    while(get_node_n_children(id)!=0) {
      d++;
      id = get_node_children(id)[0];
    }
  }
  else {
    d = -1;
  }

  // loop over links, the nodes on a link will launch swaps
  for (int i=0; i<n_max_links; i++) {
    int accepted = 0;
    int const up_cond = ud ? n_max_links-i-1 : i;
    int const dn_cond = ud ? n_max_links-i  : i+1;
    if (d == up_cond && !if_periodic(get_node_parent(local_inst))) { // my turn to swap up
      int const upstream = get_node_parent(local_inst);
      double const own_diff = ptbc_swap_dh(upstream);
      accepted = swap_link(upstream, own_diff);
    }
    else if (d == dn_cond && !if_periodic(get_node_children(local_inst)[0])){ // my turn to swap dn
      int const downstream = get_node_children(local_inst)[0];
      double const own_diff = ptbc_swap_dh(downstream);
      accepted = swap_link(downstream, own_diff);
    }
    ptbc_sync();
    if (accepted) {
      swap_arr(local_inst, Rate, 1);
      swap_arr(local_inst, &d, 1); // swap all info
      local_inst = app()->ptbc.instance_id;
    }
  }

  return;
}


/**
 * @brief attempt pbc swap with all attached tentacles
 * 
 * @param rate 
 * @return * void 
 */
void static swap_pbc(int *rate) {
  int local_inst = app()->ptbc.instance_id;
  int local_rank;
  MPI_Comm_rank(app()->mpi.comm, &local_rank);
  int my_round = -1;

  // need to know the number of synchronisations
  int periodic_id[MAX_N_DEFECTS] = {0};
  int const n_periodic = get_periodic(periodic_id);
  int n_swaps = 0;
  int max_swaps = 0;
  int order[MAX_N_DEFECTS];
  
  //periodic instance
  if (if_periodic(local_inst)) {
    int *children = get_node_children(local_inst);
    int const n_children = get_node_n_children(local_inst);

    // local rank 0 decide order
    if (local_rank == 0) {
      // find links that are not dummy
      for (int i=0; i<n_children; i++) {
        if (!if_periodic(children[i])) {
          order[n_swaps] = children[i];
          n_swaps++;
        }
      }

      // shuffle order
      shuffle(order, n_swaps);

      // notify deicision, -1 if no swap, round index if swap
      int * decision = (int *) malloc(sizeof(int) * app()->ptbc.n_instances);
      for (int i=0; i<app()->ptbc.n_instances; i++) decision[i] = -1;
      for (int i=0; i<n_swaps; i++) decision[order[i]] = i;
      
      MPI_Request* request = (MPI_Request*)malloc(n_periodic * sizeof(MPI_Request));
      int * decisions = (int *) malloc(sizeof(int) * n_periodic); // dummy receive    
      for (int i = 0; i < n_periodic; i++) {
        if (periodic_id[i] == local_inst) MPI_Iscatter(decision, 1, MPI_INT, &decisions[i], 1, MPI_INT, local_inst, leader_comm, &request[i]);
        else                              MPI_Iscatter(NULL, 1, MPI_INT, &decisions[i], 1, MPI_INT, periodic_id[i], leader_comm, &request[i]);
      }
      
      MPI_Waitall(n_periodic, request, MPI_STATUSES_IGNORE);
      free(decision);
      free(decisions);
      free(request);
    }

    // broadcast information
    MPI_Bcast(&n_swaps, 1, MPI_INT, 0, app()->mpi.comm);
    MPI_Bcast(order, n_swaps, MPI_INT, 0, app()->mpi.comm);

    // count max swaps 
    MPI_Allreduce(&n_swaps, &max_swaps, 1, MPI_INT, MPI_MAX, app()->mpi.world_comm);
  }
  // non-periodic ranks
  else {
    if (local_rank == 0) {
      // get my round number if I swap
      MPI_Request* request = (MPI_Request*)malloc(n_periodic * sizeof(MPI_Request));
      int * decisions = (int *) malloc(sizeof(int) * n_periodic);
      // receive from all periodic ranks
      for (int i=0; i<n_periodic; i++) {
        MPI_Iscatter(NULL, 1, MPI_INT, &decisions[i], 1, MPI_INT, periodic_id[i], leader_comm, request+i);
      }

      MPI_Waitall(n_periodic, request, MPI_STATUSES_IGNORE);
      free(request);

      // if selected, record my round index
      for (int i=0; i<n_periodic; i++) {
        if (decisions[i] != -1) my_round = decisions[i];
      }
      free(decisions);
    }
    // bcast info to local instance
    MPI_Bcast(&my_round, 1, MPI_INT, 0, app()->mpi.comm);
    
    // count max swaps 
    MPI_Allreduce(&n_swaps, &max_swaps, 1, MPI_INT, MPI_MAX, app()->mpi.world_comm);
  }

  // loop over max swaps, for periodic, swap if not finished. for non-periodic, swap if my_round 
  int inst= local_inst;
  for (int i=0; i<max_swaps; i++) {
    int accepted = 0;
    int partner = -1;
    if (if_periodic(inst) && i<n_swaps) {
      partner = order[i];
      double diff = ptbc_swap_dh(partner);
      accepted = swap_link(partner, diff);
    }
    else if (!if_periodic(inst) && i==my_round) {
      partner = get_node_parent(inst);
      double diff = ptbc_swap_dh(partner);
      accepted = swap_link(partner, diff);
    }
    if (accepted) {
      swap_arr(partner, rate, 1);
      swap_arr(partner, &n_swaps, 1);
      swap_arr(partner, order, MAX_N_DEFECTS);
      inst = app()->ptbc.instance_id;
    }
    ptbc_sync();
  }

  return;
}

/**
 * @brief Upstream swap, open (leaf) -> periodic (root)
 * 
 * @param Rate 
 * @return * void 
 */
void up_swap(int *Rate) {
  swap_updown_tent(Rate, 0);
  swap_pbc(Rate);
  swap_dummy(Rate);
  fill_random_numbers(); // refill random numbers for next round
}


/**
 * @brief Downstream swap, periodic (root) -> open (leaf)
 * 
 * @param Rate 
 * @return * void 
 */
void down_swap(int *Rate) {
  swap_dummy(Rate);
  swap_pbc(Rate);
  swap_updown_tent(Rate, 1);
  fill_random_numbers(); // refill random numbers for next round
}