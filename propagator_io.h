/* $Id$ */

#ifndef _PROPAGATOR_IO_H
#define _PROPAGATOR_IO_H

/* write a one flavour propagator to file */
int write_propagator(spinor * const s, spinor * const r, char * filename, 
		     const int append, const int prec);

/* write a one flavour source to file */
int write_source(spinor * const s, spinor * const r, char * filename, 
		     const int append, const int prec);

/* write a two flavour propagator to file */
int write_double_propagator(spinor * const s, spinor * const r, 
			    spinor * const p, spinor * const q,
			    char * filename, const int append, const int prec);

int write_propagator_format(char * filename, const int prec, const int no_flavours);
int write_source_format(char * filename, const int prec, const int no_flavours,
			const int t, const int lx, const int ly, const int lz,
			const int is, const int ic);

int write_propagator_type(const int type, char * filename);
int write_source_type(const int type, char * filename);

int get_propagator_type(char * filename);
int get_source_type(char * filename);

/* The binary IO */
int write_lime_spinor(spinor * const s, spinor * const r, char * filename, const int append, const int prec);
int read_lime_spinor(spinor * const s, spinor * const r, char * filename, const int position);

#endif
