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
 *  int s: signal number (not used)
 *
 *
 * void catch_del_sig(int s)
 *
 * catches some user defined signals
 * and saves configuration and 
 * random number status to disk
 *
 * input:
 *  int s: signal number (not used)
 ************************************************************/

#ifndef _SIGHANDLER_H
#define _SIGHANDLER_H
/* During critical regions one does not want */
/* the configuration to be dumped */
/* in this case set dontdump to 1 while */
/* the program is in the region */
/* don't forget to reset this value... */
extern int dontdump;

/* If a signal is catched while dontdump==1 */
/* forcedump is set to 1 */
/* This can be used to dump data to disk and */
/* exit savely after the critical region has finished */
extern int forcedump;

/* Catch an illegal instruction in order */
/* to give the user a hint what was wrong */
void catch_ill_inst(int);

/* catch some signals as SIGUSR1|2 and SIGTERM */
/* to save the current configuration and */
/* random number state */
/* This might help to save computing time */
void catch_del_sig(int);

#endif
