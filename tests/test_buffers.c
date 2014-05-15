
#include <global.h>
#include <config.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "test_buffers_gauge.h"

TEST_SUITES {
  TEST_SUITE_ADD(BUFFERS_GAUGE),
  TEST_SUITES_CLOSURE
};

int main(int argc,char *argv[]){
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
#else
  g_proc_id = 0;
#endif

  CU_SET_OUT_PREFIX("regressions/");
  CU_RUN(argc,argv);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}


