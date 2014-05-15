
#if HAVE_CONFIG_H
#include<config.h>
#endif
#include "../global.h"
#include "test_rat_init.h"

TEST_SUITES {
  TEST_SUITE_ADD(RAT),
  TEST_SUITES_CLOSURE
};

int main(int argc,char *argv[]){
  CU_SET_OUT_PREFIX("regressions/");
  CU_RUN(argc,argv);
  return 0;
}



