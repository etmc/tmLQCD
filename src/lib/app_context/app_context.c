/***********************************************************************
 *
 * Copyright (C) 2026 Roman Gruber
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
 *
 * App context module
 *
 * Author: Roman Gruber
 *         roman.gruber@unibe.ch
 *
 *******************************************************************************/

#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include "fatal_error.h"


/**
 * @brief      The rank topology struct
 */
typedef struct {
    int number_of_nodes;    // total number of nodes in the job
    int number_of_ranks;    // total number of processes in the job, i.e. size of the world communicator
    int ranks_per_node;     // number of ranks per node
    int node_index;         // index enumerating the node (unique per node)
    int node_rank;          // rank number inside the node
} RankTopology;


/**
 * @brief      Gracefully error with function name, file and line number along
 *             the error message
 *
 * @param      test  The test
 * @param      ...   Format parameters
 */
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


static void initialize(void);
static AppContext app_instance = {
    .mpi = {
        .comm = MPI_COMM_WORLD, // default communicator
        .world_comm = MPI_COMM_WORLD,
    },
    .ptbc = {
        .instance_id = 0,
        .n_instances = 1,
        .n_defects = 0,
        .active = false,
        .initialize = initialize,
        .instances = {{.active = false}},
        .defects = {{.active = false}}
    }
};


/**
 * @brief      Return the global *immutable* application context struct. To be
 *             used when reading parameters.
 *
 * @return     Global application context struct
 */
const AppContext* app(void)
{
    return &app_instance;
}


/**
 * @brief      Return the global *mutable* application context struct. To be
 *             used when initializing/setting parameters.
 *
 * @return     Global application context struct
 */
AppContext* appm(void)
{
    return &app_instance;
}


/**
 * @brief      Return rank topology.
 *
 * @return     The topology.
 */
static RankTopology get_topology(void)
{
    int world_rank;
    RankTopology topo;
    MPI_Comm node_comm, leader_comm;

    MPI_Comm_rank(app_instance.mpi.world_comm, &world_rank);
    MPI_Comm_size(app_instance.mpi.world_comm, &topo.number_of_ranks);
    MPI_Comm_split_type(app_instance.mpi.world_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
    MPI_Comm_rank(node_comm, &topo.node_rank);
    bool is_leader = topo.node_rank == 0;
    MPI_Comm_split(app_instance.mpi.world_comm, is_leader ? 0 : MPI_UNDEFINED, world_rank, &leader_comm);
    if (is_leader) MPI_Comm_size(leader_comm, &topo.number_of_nodes);
    if (is_leader) MPI_Comm_rank(leader_comm, &topo.node_index);
    MPI_Bcast(&topo.number_of_nodes, 1, MPI_INT, 0, app_instance.mpi.world_comm);
    MPI_Bcast(&topo.node_index, 1, MPI_INT, 0, node_comm);

    topo.ranks_per_node = topo.number_of_ranks/topo.number_of_nodes;

    return topo;
}


/**
 * @brief      Initializes the difference application instances.
 */
static void initialize(void)
{
    printf("\033[0;31m[PTBC] Number of chains = %d\033[0m\n", app_instance.ptbc.n_instances);

    int flag;
    MPI_Initialized(&flag);
    err(!flag, "Initialize has to be called *after* MPI_Init().");

    MPI_Comm_rank(app_instance.mpi.world_comm, &app_instance.mpi.world_rank);

    // do nothing in case of a single chain
    if (app_instance.ptbc.n_instances == 1) return;

    app_instance.ptbc.active = true;
    RankTopology topo = get_topology();

    int instance_size = topo.number_of_ranks / app_instance.ptbc.n_instances;

    err(topo.number_of_ranks % app_instance.ptbc.n_instances != 0,
        "PTBC_NCHAINS = %d must divide total number of ranks = %d",
        app_instance.ptbc.n_instances, topo.number_of_ranks);
    err(instance_size % topo.ranks_per_node != 0 && topo.ranks_per_node % instance_size != 0,
        "The number of processes per node = %d and instance_size = %d: one must be divisible by the other",
        topo.ranks_per_node, instance_size);

    // We perform a topology-aware splitting of processes into instances using
    // MPI_Comm_split. Processes within the same node should preferably be
    // associated to the same instance. We have 3 cases:
    //
    // Case 1: If we have one instance per node, there is nothing special to
    // consider.
    //
    // Case 2: If instances span multiple nodes, they should span over the
    // minimal number of nodes possible. Instances only cover whole nodes. We
    // have no notion of nodes being "close" to each other, alhtough we group
    // nodes together with adjacent node indices. Node indices are inherited
    // from rank numbers. If ranks with adjacent indices are "close", then nodes
    // with adjacent node indices are "close".
    //
    // Case 3: If we have multiple instances per node, ranks in the same
    // instance should have adjacent world rank numbers, i.e. they should be
    // "close" to each other cache-wise. No instance covers more than one node.
    int color;
    int key = app_instance.mpi.world_rank; // order of the ranks is kept
    if (instance_size == topo.ranks_per_node) { // case 1: one instance per node
        color = topo.node_index;
    } else if (instance_size % topo.ranks_per_node == 0) { // case 2: one instance spans multiple nodes
        int nodes_per_instance = instance_size / topo.ranks_per_node;
        int remainder = topo.node_index % nodes_per_instance;
        color = (topo.node_index - remainder) / nodes_per_instance; // group nodes with adjacent node indices
    } else if (topo.ranks_per_node % instance_size == 0) { // case 3: multiple instances per node
        int instances_per_node = topo.ranks_per_node / instance_size;
        int remainder = topo.node_rank % instance_size;
        int per_node_instance_index = (topo.node_rank - remainder) / instance_size;
        color = topo.node_index*instances_per_node + per_node_instance_index;
    }

    MPI_Comm_split(app_instance.mpi.world_comm, color, key, &app_instance.mpi.comm);

    int n;
    MPI_Comm_size(app_instance.mpi.comm, &n);
    err(instance_size != n, "Rank topology is not uniform");

    int instance_rank;
    MPI_Comm_rank(app_instance.mpi.comm, &instance_rank);
    app_instance.ptbc.instance_id = color;

    printf("\033[0;31m[PTBC] world rank = %d/%d in instance_id = %d/%d, as instance_rank = %d/%d\033[0m\n",
        app_instance.mpi.world_rank, topo.number_of_ranks,
        app_instance.ptbc.instance_id, app_instance.ptbc.n_instances,
        instance_rank, instance_size);


    //err(true, "bailing out");
    /*if (app_instance.ptbc.instance_id != 0) {
        char logfile[1024];
        snprintf(logfile, 1024, "logfile_%.2d.log", app_instance.ptbc.instance_id);
        freopen(logfile, "w", stdout);
    }*/

    // Every instance just changes into a subdirectory "instance_xx". Relative
    // paths work, absolute paths not.
    struct stat st = {0};
    char subdir[1024];
    snprintf(subdir, 1024, "instance_%.2d", app_instance.ptbc.instance_id);
    if (stat(subdir, &st) == -1)
        mkdir(subdir, 0700);

    chdir(subdir);
}
