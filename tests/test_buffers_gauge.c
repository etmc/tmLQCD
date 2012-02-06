#include <config.h>
#include <global.h>

#include <cu/cu.h>

#include <buffers/gauge.h>

#define EPS 5e-16

int VOLUMEPLUSRAND = 10000;

/* g_gauge_buffers is defined in one of the includes! */
extern gauge_buffers_t g_gauge_buffers;

TEST(buffers_gauge_allocate_finalize) {
  int test = 0;
  const unsigned int max = 10;

  initialize_gauge_buffers(max);

  assertEqualsM(g_gauge_buffers.max,max,"Buffers were not initialized correctly! max != 10 \n");

  finalize_gauge_buffers();
  assertFalseM(test,"TODO: No good test condition for a failed finalize.\n");

  assertEqualsM(g_gauge_buffers.allocated,0,"Finalize error, allocated != 0 \n");
  assertEqualsM(g_gauge_buffers.free,0,"Finalize error, free != 0 \n");
}

/* TODO: add test for reaching max, but since this terminates the program currently,
 *       need to wait for a cleaner failure condition to be implemented */

TEST(buffers_gauge_get_return) {
  const int max = 10;

  initialize_gauge_buffers(max);
  assertEqualsM(g_gauge_buffers.max,max,"Buffers were not initialized correctly! max != 10 \n");


  gauge_field_t gauge_field[max];

  gauge_field[0] = get_gauge_field();

  assertEqualsM(g_gauge_buffers.allocated,1,"Get error, allocated != 1 \n");
  assertEqualsM(g_gauge_buffers.free,0,"Get error, free != 0 \n");

  for(int i = 1; i < 6; ++i) {
    gauge_field[i] = get_gauge_field();
  }

  assertEqualsM(g_gauge_buffers.allocated,6,"Get error, allocated != 6 \n");

  return_gauge_field(&gauge_field[5]);

  assertEqualsM(gauge_field[5].field,NULL,"Return error, field pointer not NULL \n");
  assertEqualsM(g_gauge_buffers.allocated,6,"Return error, allocated != 6 \n");
  assertEqualsM(g_gauge_buffers.free,1,"Return error, free != 1 \n");

  gauge_field[5] = get_gauge_field();

  assertNotEqualsM(gauge_field[5].field,NULL,"Get error, field pointer still NULL \n");
  assertEqualsM(g_gauge_buffers.allocated,6,"Get error, allocated != 6 \n");
  assertEqualsM(g_gauge_buffers.free,0,"Get error, free != 0 \n");

  allocate_gauge_buffers(2);

  assertEqualsM(g_gauge_buffers.allocated,8,"Allocate error, allocated != 8 \n");
  assertEqualsM(g_gauge_buffers.free,2,"Allocate error, free != 2 \n");

  gauge_field[6] = get_gauge_field();

  assertEqualsM(g_gauge_buffers.allocated,8,"Get error, allocated != 8 \n");
  assertEqualsM(g_gauge_buffers.free,1,"Get error, free != 1 \n");

  for(int i = 4; i <= 6; ++i) {
    return_gauge_field(&gauge_field[i]);
  }

  assertEqualsM(g_gauge_buffers.free,4,"Return error, free != 4 \n");

  free_unused_gauge_buffers();

  assertEqualsM(g_gauge_buffers.allocated,4,"Free error, allocated != 4 \n");
  assertEqualsM(g_gauge_buffers.free,0,"Free error, free != 0 \n");

  for(int i = 0; i < 4; ++i) {
    return_gauge_field(&gauge_field[i]);
  }
  
  finalize_gauge_buffers();

  assertEqualsM(g_gauge_buffers.allocated,0,"Finalize error, allocated != 0 \n");
  assertEqualsM(g_gauge_buffers.free,0,"Finalize error, free != 0 \n");
}

