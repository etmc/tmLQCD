/* $Id$ */
#ifndef _FORTRAN_MY_H
#define _FORTRAN_MY_H

#if (defined NOF77UNDERSCORE || defined NOF77_)
#define _FT(s) s
#else
#define _FT(s) s ## _
#endif

#endif
