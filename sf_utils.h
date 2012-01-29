
#ifndef _SF_UTILS_H
#define _SF_UTILS_H

/* This is the standard square plaquette term. */
/* It calculates the sum over all squares on the lattice */
/* and then divides by 3.  There is no coupling dependence. */
double calc_sq_plaq( void );

/* This is the same as above, but the plaquettes that involve */
/* the time-like links from sites at t=0 to t=1, t=T-2 to t=T-1 */
/* and t=T-1 to t=0 are ignored.  Also all space-space squares */
/* for t=0 and t=T-1 are ignored. */
double calc_bulk_sq_plaq( void );

double calc_boundary_space_space_sq_plaq( void );

double calc_boundary_space_time_sq_plaq( void );

/* It is designed so that bulk + boundary + wrapped = total square contribution. */
double calc_wrapped_sq_plaq( void );

/* This is the standard rectangle term. */
/* It calculates the sum over all rectangles on the lattice */
/* and then divides by 3. */
double calc_rect_plaq( void );

#endif
