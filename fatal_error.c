#include <stdio.h>
#include <global.h>

#ifdef MPI
#include <mpi.h>
#endif

void fatal_error(char const *error, char const *function)
{
  if (error != NULL)
  {
    fprintf(stderr, "FATAL ERROR\n");
    if (function != NULL)
    {
#ifdef MPI
      fprintf(stderr, "  Within %s (reported by node %d):\n", function, g_proc_id);
#else
      fprintf(stderr, "  Within %s:\n", function);
#endif
    }
    fprintf(stderr, "    %s\n", error);
    fflush(stderr);
  }
  
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Finalize();
#endif

  exit(500);
}
