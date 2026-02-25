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
#include <mpi.h>
#include "fatal_error.h"


static AppContext app_instance = {
    .mpi = {
        .comm = MPI_COMM_WORLD // default communicator
    }
};


const AppContext* app(void)
{
    return &app_instance;
}


void app_context_init(const MPI_Comm comm)
{
    static bool initialized = false;

    if (initialized) fatal_error("Application context already initialized", __func__);

    app_instance.mpi.comm = comm;
    initialized = true;
}


void app_context_finalize(void)
{
    
}
