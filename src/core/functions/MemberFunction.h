#ifndef MemberFunction_H
#define MemberFunction_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

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
    template<class memberObjectType, class valueType>
    class MemberFunction : public TypedFunction<valueType> {
        
    public:
        // constructors and destructor
        MemberFunction(const std::string &n, const TypedDagNode<memberObjectType> *o, const std::vector<const DagNode* > &a);
        virtual                                    ~MemberFunction(void);
        
        // public methods
        MemberFunction<memberObjectType,valueType>* clone(void) const;                                                              //!< Clone the function
        void                                        update(void);                                                                   //!< Update the value of the function
        
    protected:
        
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Implementation of swaping parameters

    private:
        
        std::string                                 method_name;
        const TypedDagNode<memberObjectType>*       the_member_variable;
        std::vector<const DagNode* >                args;
    };
    
}

#include "Cloneable.h"
#include "Cloner.h"
#include "IsDerivedFrom.h"


template <class memberObjectType, class valueType>
RevBayesCore::MemberFunction<memberObjectType,valueType>::MemberFunction(const std::string &n, const TypedDagNode<memberObjectType> *o, const std::vector<const DagNode* > &a) : TypedFunction<valueType>( nullptr ),
    method_name( n ),
    the_member_variable( o ),
    args( a )
{
    if constexpr (not IsAbstract<valueType>::Is)
	this->value = new valueType;

    this->addParameter( the_member_variable );
    typename std::vector<const DagNode* >::iterator it;
    for (it = args.begin(); it != args.end(); ++it)
    {
        this->addParameter( *it );
    }
    
}

template <class memberObjectType, class valueType>
RevBayesCore::MemberFunction<memberObjectType,valueType>::~MemberFunction( void )
{
    
}



template <class memberObjectType, class valueType>
RevBayesCore::MemberFunction<memberObjectType,valueType>* RevBayesCore::MemberFunction<memberObjectType, valueType>::clone( void ) const
{
    
    return new MemberFunction<memberObjectType,valueType>(*this);
}



template <class memberObjectType, class valueType>
void RevBayesCore::MemberFunction<memberObjectType,valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( the_member_variable == oldP )
    {
        the_member_variable = static_cast< const TypedDagNode<memberObjectType>* > (newP);
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



template <class memberObjectType, class valueType>
void RevBayesCore::MemberFunction<memberObjectType,valueType>::update( void )
{
    const memberObjectType& the_object = dynamic_cast<const memberObjectType& >( the_member_variable->getValue() );
    
    if constexpr (IsAbstract<valueType>::Is)
    {
	delete this->value;
	the_object.executeMethodAbstract(method_name,args,this->value);
    }
    else
    {
	// The memberObjectType might be Core::MatrixBoolean, Core::RbVector,
	// But Core::MatrixBoolean and Core::RbVector don't have the right update types.
	the_object.executeMethod(method_name,args,*this->value);

//	const MemberObject<valueType>& the_member_object = dynamic_cast<const MemberObject<valueType>& >( the_member_variable->getValue() );
//	the_member_object.executeMethod(method_name,args,*this->value);
    }
    
}


#endif
