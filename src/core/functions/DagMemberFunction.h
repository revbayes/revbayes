#ifndef DagMemberFunction_H
#define DagMemberFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
    
    template<class valueType>
    class MemberObject;
    
    
    /**
     * \brief Member function that maps a function call to a method of a member object.
     *
     * A member function is a function that is used inside deterministic nodes and provides
     * access to a member method of a member object.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2014-07-17
     *
     */
    template<class valueType>
    class DagMemberFunction : public TypedFunction<valueType> {
        
    public:
        // constructors and destructor
        DagMemberFunction(const std::string &n, const DagNode *o, const std::vector<const DagNode* > &a);
        virtual                                        ~DagMemberFunction(void);
        
        // pure virtual public methors
        DagMemberFunction<valueType>*                   clone(void) const;                                                              //!< Clone the function
        void                                            update(void);                                                                   //!< Update the value of the function
        
    protected:
        
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        std::string                                     method_name;
        const DagNode*                                  the_member_variable;
        std::vector<const DagNode* >                    args;
    };
    
}

#include "Cloneable.h"
#include "Cloner.h"
#include "DagNode.h"
#include "IsDerivedFrom.h"


template <class valueType>
RevBayesCore::DagMemberFunction<valueType>::DagMemberFunction(const std::string &n, const DagNode *o, const std::vector<const DagNode* > &a) : TypedFunction<valueType>( new valueType() ),
    method_name( n ),
    the_member_variable( o ),
    args( a )
{
    
    this->addParameter( the_member_variable );
    typename std::vector<const DagNode* >::iterator it;
    for (it = args.begin(); it != args.end(); ++it)
    {
        this->addParameter( *it );
    }
    
}


template <class valueType>
RevBayesCore::DagMemberFunction<valueType>::~DagMemberFunction( void )
{
    
}


template <class valueType>
RevBayesCore::DagMemberFunction<valueType>* RevBayesCore::DagMemberFunction<valueType>::clone( void ) const
{
    
    return new DagMemberFunction<valueType>(*this);
}



template <class valueType>
void RevBayesCore::DagMemberFunction<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if ( the_member_variable == oldP )
    {
        the_member_variable = newP;
    }
    else
    {
        for (size_t i = 0; i < args.size(); ++i)
        {
            if (oldP == args[i])
            {
                args[i] = newP;
                // we can jump out of the loop now
                break;
            }
        }
    }
    
}



template <class valueType>
void RevBayesCore::DagMemberFunction<valueType>::update( void )
{
    
    const MemberObject<valueType> *member_object = dynamic_cast<const MemberObject<valueType>* >( the_member_variable );
    member_object->executeMethod(method_name,args,*this->value);
    
}


#endif
