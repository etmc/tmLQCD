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

/* Minimal test of the PTBC rng:
 *   1. every rank of every instance starts from a distinct rng state, also
 *      distinct from the PTBC stream on world rank 0
 *   2. swap_rng() hands the rng state over to the partner instance
 *
 * Only the rng-relevant part of the app context is filled, no input file is read.
 * Run with a multiple of N_INSTANCES ranks in a fresh directory, e.g.
 *   mpirun -n 4 ./ptbc_rng_test
 * Results go to stderr, since stdout is redirected to instance_xx/hmc.log.
 */

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "global.h"
#include "read_input.h"
#include "start.h"
#include "ranlxs.h"
#include "ranlxd.h"
#include "mpi_init.h"
#include "ptbc.h"

#define N_INSTANCES 4
#define PTBC_SEED 4242
#define SEED 127782
#define RLX_SIZE 105                 // rlxs_size() == rlxd_size()
#define STATE_SIZE (2 * RLX_SIZE)    // rlxs state followed by rlxd state

static void get_state(int state[STATE_SIZE]) {
  rlxs_get(state);
  rlxd_get(state + RLX_SIZE);
}

// combine the result of all ranks, report on world rank 0
static bool report(char const *name, bool const ok) {
  int local = ok, global;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LAND, app()->mpi.world_comm);
  if (app()->mpi.world_rank == 0) fprintf(stderr, "[ptbc_rng_test] %-16s %s\n", name, global ? "PASS" : "FAIL");
  return global;
}

/**
 * @brief      Every rank holds distinct rlxs and rlxd states, and no rlxd state equals
 *             the PTBC stream on world rank 0.
 */
static bool test_distinct_states(void) {
  int n_ranks;
  MPI_Comm_size(app()->mpi.world_comm, &n_ranks);

  int state[STATE_SIZE];
  get_state(state);

  int *all = NULL;
  if (app()->mpi.world_rank == 0) all = malloc(sizeof(int) * STATE_SIZE * n_ranks);
  MPI_Gather(state, STATE_SIZE, MPI_INT, all, STATE_SIZE, MPI_INT, 0, app()->mpi.world_comm);

  bool ok = true;
  if (app()->mpi.world_rank == 0) {
    for (int i = 0; i < n_ranks; i++) {
      int const *si = all + i * STATE_SIZE;
      if (memcmp(si + RLX_SIZE, app()->ptbc.rng_state, RLX_SIZE * sizeof(int)) == 0) {
        fprintf(stderr, "  world rank %d: rlxd state equals the PTBC stream\n", i);
        ok = false;
      }
      for (int j = i + 1; j < n_ranks; j++) {
        int const *sj = all + j * STATE_SIZE;
        if (memcmp(si, sj, RLX_SIZE * sizeof(int)) == 0) {
          fprintf(stderr, "  world ranks %d and %d: same rlxs state\n", i, j);
          ok = false;
        }
        if (memcmp(si + RLX_SIZE, sj + RLX_SIZE, RLX_SIZE * sizeof(int)) == 0) {
          fprintf(stderr, "  world ranks %d and %d: same rlxd state\n", i, j);
          ok = false;
        }
      }
    }
    free(all);
  }

  return report("distinct states", ok);
}

/**
 * @brief      Swap rng states with instances on a ring, periodic in the index.
 *             Round 0 pairs (0,1),(2,3),...; round 1 pairs (1,2),(3,4),...,(n-1,0).
 *             After swap_rng() each rank must hold the state its partner had before,
 *             and the PTBC stream on world rank 0 must stay untouched.
 */
static bool test_swap(void) {
  int const n = app()->ptbc.n_instances;
  int const id = app()->ptbc.instance_id;
  int n_ranks, local_rank;
  MPI_Comm_size(app()->mpi.world_comm, &n_ranks);
  MPI_Comm_rank(app()->mpi.comm, &local_rank);

  bool ok = true;
  for (int round = 0; round < 2; round++) {
    int const partner = (id % 2 == round % 2) ? (id + 1) % n : (id - 1 + n) % n;
    int const partner_rank = partner * (n_ranks / n) + local_rank;  // instances are block-contiguous

    // reference exchange, independent of swap_rng
    int mine[STATE_SIZE], expected[STATE_SIZE], got[STATE_SIZE];
    get_state(mine);
    MPI_Sendrecv(mine, STATE_SIZE, MPI_INT, partner_rank, 0, expected, STATE_SIZE, MPI_INT,
                 partner_rank, 0, app()->mpi.world_comm, MPI_STATUS_IGNORE);

    int ptbc_before[RLX_SIZE];
    memcpy(ptbc_before, app()->ptbc.rng_state, sizeof(ptbc_before));

    swap_rng(partner);
    get_state(got);

    if (memcmp(got, expected, sizeof(got)) != 0) {
      fprintf(stderr, "  round %d: world rank %d did not receive the state of instance %d\n",
              round, app()->mpi.world_rank, partner);
      ok = false;
    }
    if (app()->mpi.world_rank == 0 && memcmp(ptbc_before, app()->ptbc.rng_state, sizeof(ptbc_before)) != 0) {
      fprintf(stderr, "  round %d: PTBC stream on world rank 0 changed by swap_rng\n", round);
      ok = false;
    }
  }

  return report("periodic swap", ok);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // only the rng-relevant part of the app context
  appm()->ptbc.n_instances = N_INSTANCES;
  appm()->ptbc.seed = PTBC_SEED;
  app()->ptbc.initialize();  // split into instances, set instance_id, chdir to instance_xx

  // globals otherwise set by read_input and tmlqcd_mpi_init
  MPI_Comm_size(app()->mpi.comm, &g_nproc);
  MPI_Comm_rank(app()->mpi.comm, &g_proc_id);
  rlxd_level = 1;

  init_ptbc_rng();  // base ranks and PTBC stream on world rank 0
  start_ranlux(rlxd_level, SEED);

  bool ok = test_distinct_states();
  ok = test_swap() && ok;

  MPI_Finalize();
  return ok ? 0 : 1;
}
