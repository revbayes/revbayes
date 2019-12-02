/**
 * @file Cloneable.h
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

#ifndef Cloneable_H
#define Cloneable_H

#include <string>

namespace RevBayesCore {
    
    class Cloneable {
        
    public:
        virtual                         ~Cloneable(void) {}
        
        virtual Cloneable*              clone( void ) const = 0;                                    //!< Create a clone/copy of the object
    };
    
}

#endif

