/**********************************************************
 * 
 * exchange routines for the borders of a timeslice of spinor fields
 *
 * Author: Luigi Scorzato
 *
 **********************************************************/

#ifndef _XCHANGE_FIELDTS_H
#define _XCHANGE_FIELDTS_H

#define EVEN 1 
#define  ODD 0 

#ifdef MPI
void xchange_field_open(spinor * const , const int , const int , MPI_Request * , MPI_Status *);  
void xchange_field_close(MPI_Request * , MPI_Status * , int );
void xchange_field_slice(spinor * const , const int , const int );
#endif

#endif
