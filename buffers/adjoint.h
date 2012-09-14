#ifndef _BUFFERS_ADJOINT_H
#define _BUFFERS_ADJOINT_H

#include <buffers/alignment.h>
#include <su3adj.h>

typedef su3adj su3adj_tuple[4];

#include "template_header.inc"

__DECLARE_BUFFER_INTERFACE(su3adj_tuple, adjoint)

#undef __DECLARE_BUFFER_INTERFACE

#endif

