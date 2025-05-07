/**
 * @file
 * This file contains the declaration of BinaryState, which is
 * the class for the Standard data types in RevBayes.
 *
 * @brief Declaration of BinaryState
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-05-24 09:58:04 +0200 (Thu, 24 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: BinaryState.h 1568 2012-05-24 07:58:04Z hoehna $
 */

#ifndef BinaryState_H
#define BinaryState_H

#include <cstddef>
#include <iosfwd>

#include "StandardState.h"

namespace RevBayesCore {

    class BinaryState : public StandardState {
    
    public:
                                        BinaryState(size_t n=2);                                        //!< Default constructor
                                        BinaryState(const std::string &s);
    };

}

#endif

