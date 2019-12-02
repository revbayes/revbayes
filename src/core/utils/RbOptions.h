#ifndef RbOptions_H
#define RbOptions_H

/* Debug switches */
/* It is useful to list the switches here but it is preferable to switch
   the defines on in the IDE rather than by uncommenting them here, so
   that accidental commits do not disturb other developers. Beware! */
//#define ASSERTIONS_ALL
//#define ASSERTIONS_TREE
//#define ASSERTIONS_DISTRIBUTIONS
//#define DEBUG_ALL
//#define DEBUG_BISON_FLEX
//#define RB_MPI        // Allows use of MPI (mpi.h) features

//#define TESTING

/* Feature enabling switches */
#ifndef RB_ARM
#define SSE_ENABLED
#endif
//#define AVX_ENABLED


/* Test whether we should use linenoise */
#if !defined (NO_LINENOISE)
#define USE_LIB_LINENOISE
#endif

/* Test whether we need to debug everything. */
#if defined (DEBUG_ALL)

    // switch all assertions on
    #ifndef ASSERTIONS_ALL
    #define ASSERTIONS_ALL
    #endif

    // switch debugging parser on
    //#ifndef DEBUG_BISON_FLEX
    //#define DEBUG_BISON_FLEX
    //#endif


#endif




/* Test whether we need to debug everything. */
#if defined (ASSERTIONS_ALL)

    // switch all assertions on
    #ifndef ASSERTIONS_DISTRIBUTIONS
    #define ASSERTIONS_DISTRIBUTIONS
    #endif

    #ifndef ASSERTIONS_TREE
    #define ASSERTIONS_TREE
    #endif

#endif


//#endif


// AdmixtureGraph depends on armadillo for linear algebra
// Uncomment the first line to enable the armadillo library
//#define USE_LIB_ARMADILLO
#ifdef USE_LIB_ARMADILLO
#include <armadillo>
#endif

#endif
