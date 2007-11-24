/* $Id$ */

#ifndef _PROPAGATOR_IO_H
#define _PROPAGATOR_IO_H

int write_propagator(spinor * const s, spinor * const r, char * filename, 
		     const int append, const int prec);

int write_double_propagator(spinor * const s, spinor * const r, 
			    spinor * const p, spinor * const q,
			    char * filename, const int append, const int prec);

int write_propagator_format(char * filename, const int prec, const int no_flavours);
int write_propagator_type(const int type, char * filename);
int get_propagator_type(char * filename);
int write_lime_spinor(spinor * const s, spinor * const r, char * filename, const int append, const int prec);
int read_lime_spinor(spinor * const s, spinor * const r, char * filename, const int position);

#endif
