#ifndef MemberObject_H
#define MemberObject_H

#include <string>
#include <vector>

namespace RevBayesCore {
    
    class DagNode;
    
    /**
     * \brief Interface for member-objects
     *
     * A member object is an object that provides convenient access to member methods.
     * These member methods can be called via a MemberFunction.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2014-07-17
     *
     */
    template <class valueType>
    class MemberObject {
        
    public:
        virtual                         ~MemberObject(void) {}
        
        virtual void                    executeMethod( const std::string &n, const std::vector<const DagNode*> &args, valueType &retValue) const = 0;       //!< Execute the member-method
    };
    
}

#endif

