/**
 * @file RbMathVector
 * This file contains the math on vectors.
 *
 * @brief Implementation of vector algebra.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */


#ifndef RbMathVector_H
#define RbMathVector_H


#include <vector>

namespace RevBayesCore {
    
    namespace RbMath {
    
        void                        normalize(std::vector<double>& x, double sum);                                      //!< Normalize a vector so that its sum is sum
        void                        normalizeToMin(std::vector<double>& x, double sum);                                      //!< Normalize a vector so that its sum is sum

    }
    
}

#endif
