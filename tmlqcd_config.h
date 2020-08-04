#ifndef TMLQCD_CONFIG_H
#define TMLQCD_CONFIG_H

/* This file exists to deal with the issue that we include packages
 * which do some rather dubious things with autoconf-generated variables.
 * In particular, this is for compatibility with c-lime, which
 * #undef's PACKAGE PACKAGE_BUGREPORT PACKAGE_NAME PACKAGE_STRING PACKAGE_TARNAME PACKAGE_VERSION VERSION
 * Depending on the include order, either our definitions in tmlqcd_config.h survive,
 * or they do not...
 */

 #undef PACKAGE
 #undef PACKAGE_BUGREPORT
 #undef PACKAGE_NAME
 #undef PACKAGE_STRING
 #undef PACKAGE_TARNAME
 #undef PACKAGE_VERSION

 // include the autoconf-generated header
 #include "tmlqcd_config_internal.h"
 static const char* const TMLQCD_PACKAGE_BUGREPORT = PACKAGE_BUGREPORT;
 static const char* const TMLQCD_PACKAGE_NAME = PACKAGE_NAME;
 static const char* const TMLQCD_PACKAGE_STRING = PACKAGE_STRING;
 static const char* const TMLQCD_PACKAGE_TARNAME = PACKAGE_TARNAME;
 static const char* const TMLQCD_PACKAGE_VERSION = PACKAGE_VERSION;

 // unlike lime, we don't undef anything

#endif
