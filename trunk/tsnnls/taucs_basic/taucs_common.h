/* Taucs_common.

   This file added to the taucs_basic distribution in TSNNLS to
   include the dmalloc header in a TAUCS compile cycle, hoping to pick
   up a subtle porting bug. 

   Jason Cantarella. 4/1/2007 */

#ifndef TAUCS_COMMON
#define TAUCS_COMMON

#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

#endif
