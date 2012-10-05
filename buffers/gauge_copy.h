#ifndef _BUFFERS_GAUGE_H
#define _BUFFERS_GAUGE_H

#include <buffers/alignment.h>
#include <su3.h>

#ifndef ALIGN_BASE
#  define ALIGN_BASE 0x0f
#endif

typedef su3 su3_c_tuple[8];

#include "template_header.inc"

__DECLARE_BUFFER_INTERFACE(su3_c_tuple, gauge_copy)

#undef __DECLARE_BUFFER_INTERFACE

