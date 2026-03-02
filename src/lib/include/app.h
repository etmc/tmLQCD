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
 * Simple MPI header wrapper
 *
 * Author: Roman Gruber
 *         roman.gruber@unibe.ch
 *
 *******************************************************************************/

#ifndef APP_H
#define APP_H


#include <stdbool.h>


#if defined(TM_USE_MPI)
#include <mpi.h>
#endif


#ifndef MPI_VERSION
typedef int MPI_Comm;
#endif


#ifdef __cplusplus
extern "C" {
#endif


#define MAX_N_DEFECTS 10
#define MAX_N_INSTANCES 10


typedef enum direction_t {
   DIRECTION_T = 0,
   DIRECTION_X = 1,
   DIRECTION_Y = 2,
   DIRECTION_Z = 3
} direction_t;


typedef struct {
    MPI_Comm comm;          // MPI instance communicator
    MPI_Comm world_comm;    // MPI world communicator
    int world_rank;         // MPI world rank
} MPIContext;


typedef struct {
    bool active;        // Whether the defect is active or not
    int Ld[3];          // Extents of the defect
    direction_t along;  // Along which dimension
} PTBCDefect;


typedef struct {
    bool active;          // Whether the instance is active or not
    int n_coeffs;         // Number of coefficients / defect this instance is associated to
    PTBCDefect** defects; // List of defects where this instance is associated to
    double* coefficients; // List of coefficients for the defects
} PTBCInstance;


typedef struct {
    bool active;                                // Whether PTBC mode is active or not
    int instance_id;                            // Instance ID
    int n_instances;                            // Number of instances
    int n_defects;                              // Number of defects
    PTBCInstance instances[MAX_N_INSTANCES];    // List of all instances
    PTBCDefect defects[MAX_N_DEFECTS];          // List of all defects
    void (*initialize)(void);                   // PTBC algorithm initializer
} PTBCContext;


typedef struct {
    MPIContext mpi;     // MPI context
    PTBCContext ptbc;   // PTBC context
} AppContext;


const AppContext* app(void);
AppContext* appm(void);


#ifdef __cplusplus
}
#endif

#endif
