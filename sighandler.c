/* $Id$ */

/************************************************************
 *
 * Routines to handle system signals
 *
 * void catch_ill_inst(int s)
 *
 * catches illegal instructions signal
 * and writes an error indication to
 * stdout.
 *
 * input:
 *  int s: signal number (not needed)
 *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "io.h"
#include "ranlxd.h"

/* quark propagators defined in main.c */
extern spinor **** qprop;

/* During critical regions one does not want */
/* the configuration to be dumped */
/* in this case set dontdump to 1 while */
/* the program is in the region */
/* don't forget to reset this value... */
int dontdump=0;

/* If a signal is catched while dontdump==1 */
/* forcedump is set to 1 */
/* This can be used to dump data to disk and */
/* exit savely after the critical region has finished */
int forcedump=0;


/* Catch an illegal instruction in order */
/* to give the user a hint what was wrong */
void catch_ill_inst(int s){
  printf("An illegal instruction occured!\n");
#ifdef SSE
  printf("Your code was compiled to use SSE1 instructions.\n");
#endif
#ifdef SSE2
  printf("Your code was compiled to use SSE2 instructions.\n");
#endif
  printf("Probably this caused the exception.\n");
  printf("Please check whether you processor understands SSE (1 or 2) instructions!\n");
  printf("Aborting...\n");
  fflush(stdout);
  exit(0);
}

extern su3 gauge_tmp[VOLUME][4] ALIGN;

/* catch some signals as SIGUSR1|2 and SIGTERM */
/* to save the current configuration and */
/* random number state */
/* This might help to save computing time */
void catch_del_sig(int s){
  int ix, mu;
  int rlxd_state[105];
  FILE * rlxdfile = NULL;
  su3 *v,*w;

  if(dontdump==0){
    printf("We got a signal SIGUSR1|2 or SIGTERM!\n");
    printf("Dumping configuration and Rlxd state to disk!\n");
    printf("Exit savely!\n");
    fflush(stdout);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    /* Save all the stuff needed for restarting */
    if(ranlxd_init == 1){
      if(g_proc_id==0) {
	rlxd_get(rlxd_state);
	rlxdfile=fopen("last_state","w");
	fwrite(rlxd_state,sizeof(rlxd_state),1,rlxdfile);
	fclose(rlxdfile);
      }
    }
    else{
      printf("ranlxd was not initialized, thus nothing to dump...\n");
    }
    for(ix=0;ix<VOLUME;ix++) {
      for(mu=0;mu<4;mu++){
	/* Auch MIST */
	v=&g_gauge_field[ix][mu];
	w=&gauge_tmp[ix][mu];
	_su3_assign(*v,*w);
      }
    }
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    write_gauge_field_time_p("last_configuration");
    fflush(stdout);
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(0);
  }
  else{
    printf("\nWaiting to dump until critical stuff is done...!\n\n");
    fflush(stdout);
    forcedump=1;
  }
}
