/**
 * @file Assign.h
 *
 *
 * @brief Declaration of Assign. This class clones objects by either calling clone, if the
 * object is derived from Cloneable, or calling the copy constructor. Hence, the Assign can be
 * used in templates when it is unknown if either the copy constructor (basic classes)
 * or the clone function (abstract classes) exists.
 *
 * (c) Copyright 2009-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 */

#ifndef Assign_H
#define Assign_H

#include <string>

namespace RevBayesCore {
    
    template <typename objType, int>
    // general case: T is not derived from Assignable
    // calls copy constructor
    class Assign {
        
    public:
        //!< Assign the given object.
        static void                     doAssign( objType &a, const objType &b ) { a = b; }
    };
    
    template <typename objType>
    // T is derived from Assignable
    // calls clone
    class Assign<objType,1> {
        
    public:
        //!< Assign the given object.
        static void                     doAssign( objType &a, const objType &b ) { a.assign(b); }
    };
    
}



#endif

