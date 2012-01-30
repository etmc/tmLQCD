#pragma once

#include <su3.h>

#ifndef ALIGN_BASE
#  define ALIGN_BASE 0x0f
#endif

typedef su3 su3_tuple[4];

#include "template_header.inc"

__DECLARE_BUFFER_INTERFACE(su3_tuple, gauge)

#undef __DECLARE_BUFFER_INTERFACE
