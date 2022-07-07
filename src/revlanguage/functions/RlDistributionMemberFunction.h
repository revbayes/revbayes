#ifndef RlDistributionMemberFunction_H
#define RlDistributionMemberFunction_H

#include "RlMemberMethod.h"
#include "RlTypedFunction.h"

namespace RevLanguage {
    
    template <typename distributionType, typename retType>
    class DistributionMemberFunction : public TypedFunction<retType>, public MemberMethod {
        
    public:
        DistributionMemberFunction(const std::string &n, const typename distributionType::rlValueType* o, ArgumentRules* argRules, bool f=false);                                 //!< Constructor
        virtual ~DistributionMemberFunction(){};
        
        // Basic utility functions
        DistributionMemberFunction*                                     clone(void) const;                              //!< Clone the object
        static const std::string&                                       getClassType(void);                             //!< Get class name
        static const TypeSpec&                                          getClassTypeSpec(void);                         //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                    //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                        //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< typename retType::valueType>*      createFunction(void) const;                     //!< Execute function
        void                                                            setMemberObject(const RevPtr<RevVariable> &obj);//!< Set the member object to which this function belongs
        const ArgumentRules&                                            getArgumentRules(void) const;                   //!< Get argument rules
        void                                                            setMethodName(const std::string& name);         //!< Set name of member method
        
        
    private:
        
        const ArgumentRules*                                            argument_rules;                                  //!< Argument rules (different for different member functions)
        std::string                                                     func_name;                                       //!< Name of member method
        RevPtr<RevVariable>                                             object;
        const typename distributionType::rlValueType*                   the_member_object;
        bool                                                            force_updates;
        
    };
    
    
    
}


#include "DistributionMemberFunction.h"
#include "RbException.h"

/** default constructor */
template <typename distributionType, typename retType>
RevLanguage::DistributionMemberFunction<distributionType, retType>::DistributionMemberFunction( const std::string &n, const typename distributionType::rlValueType *o, ArgumentRules* ar, bool f ) : TypedFunction<retType>(  ),
    argument_rules( ar ),
    func_name( n ),
    object( NULL ),
    the_member_object( o ),
    force_updates( f )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename distributionType, typename retType>
RevLanguage::DistributionMemberFunction<distributionType, retType>* RevLanguage::DistributionMemberFunction<distributionType, retType>::clone( void ) const
{
    
    return new DistributionMemberFunction<distributionType, retType>( *this );
}


template <typename distributionType, typename retType>
RevBayesCore::TypedFunction< typename retType::valueType >* RevLanguage::DistributionMemberFunction<distributionType, retType>::createFunction( void ) const
{
    
    std::vector<const RevBayesCore::DagNode*> arg_nodes;
    for ( size_t i=0; i<this->args.size(); ++i)
    {
        arg_nodes.push_back( this->args[i].getVariable()->getRevObject().getDagNode() );
    }
    
    const RevBayesCore::TypedDagNode<typename distributionType::rbValueType>* tmp = the_member_object->getDagNode();
    const RevBayesCore::StochasticNode<typename distributionType::rbValueType>* o = dynamic_cast<const RevBayesCore::StochasticNode<typename distributionType::rbValueType>* >( tmp );
    if ( o == NULL )
    {
        throw RbException("Could not cast the member object.");
    }
    
    RevBayesCore::DistributionMemberFunction<typename distributionType::rbValueType, typename retType::valueType> *func = new RevBayesCore::DistributionMemberFunction<typename distributionType::rbValueType, typename retType::valueType>(this->func_name, o, arg_nodes);
    func->setForceUpdates( force_updates );
    
    return func;
}


/** Get argument rules */
template <typename memberObjectType, typename retType>
const RevLanguage::ArgumentRules& RevLanguage::DistributionMemberFunction<memberObjectType, retType>::getArgumentRules(void) const
{
    
    return *argument_rules;
}



template <typename distributionType, typename retType>
const std::string& RevLanguage::DistributionMemberFunction<distributionType, retType>::getClassType(void)
{
    
    static std::string revClassType = "DistributionMemberFunction<" + distributionType::getClassType() + "," + retType::getClassType() + ">";
    
	return revClassType;
}

/* Get class type spec describing type of object */
template <typename distributionType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::DistributionMemberFunction<distributionType, retType>::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return revClassTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename distributionType, typename retType>
std::string RevLanguage::DistributionMemberFunction<distributionType, retType>::getFunctionName( void ) const
{
    
    return func_name;
}


template <typename distributionType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::DistributionMemberFunction<distributionType, retType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


template <typename distributionType, typename retType>
void RevLanguage::DistributionMemberFunction<distributionType, retType>::setMemberObject( const RevPtr<RevVariable> &obj)
{
    
    // we do not own the object itself because one object can have multiple member functions
    object = obj;
    the_member_object = static_cast< const typename distributionType::rlValueType *>( &(obj->getRevObject()) );
    
}


#endif

