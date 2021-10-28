/**
 * @file
 * This file contains declarations of helper functions for BEAGLE.
 *
 * @brief Declaration of helper functions for the BEAGLE library.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author Daniel L. Ayres
 * @license GPL version 3
 *
 * $Id$
 */

#ifndef BeagleUtilities_h
#define BeagleUtilities_h

#include <string>

namespace RevBayesCore {
    
    namespace BeagleUtilities {
        
        std::string printBeagleFlags(long inFlags); //!< Print BEAGLE capability flags

    }
}

#endif /* BeagleUtilities_h */

