/**
 * @file Assignable.h
 *
 *
 * @brief Declaration of template cloneable. The only method this interface defines
 * is T* clone(void) const. The clone function is required in e.g. DAG nodes if the
 * type is abstract (e.g. RateMatrix, Function, ...)
 *
 * (c) Copyright 2009-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 */

#ifndef Assignable_H
#define Assignable_H

#include <string>

namespace RevBayesCore {
    
    class Assignable {
        
    public:
        virtual                         ~Assignable(void) {}
        
        virtual Assignable&              assign( const Assignable &a ) = 0;         //!< Assign the object and prevent slicing
    };
    
}

#endif

