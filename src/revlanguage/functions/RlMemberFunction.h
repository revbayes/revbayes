#ifndef RlMemberFunction_H
#define RlMemberFunction_H

#include "RlMemberMethod.h"
#include "RlTypedFunction.h"


namespace RevLanguage {
    
    template <typename memberObjectType, typename retType>
    class MemberFunction : public TypedFunction<retType>, public MemberMethod {
        
    public:
        MemberFunction(const std::string &name, const memberObjectType* o, ArgumentRules* argRules);                                             //!< Constructor
        MemberFunction(const MemberFunction &m);                                                                        //!< Constructor
        virtual ~MemberFunction(void);
        
        MemberFunction&                                                 operator=(const MemberFunction &m);             //!< Assignment operator
        
        // Basic utility functions
        MemberFunction*                                                 clone(void) const;                              //!< Clone the object
        static const std::string&                                       getClassType(void);                             //!< Get class name
        static const TypeSpec&                                          getClassTypeSpec(void);                         //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                        //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< typename retType::valueType>*      createFunction(void) const;                     //!< Execute function
        void                                                            setMemberObject(const RevPtr<RevVariable> &obj);   //!< Set the member object to which this function belongs
        const ArgumentRules&                                            getArgumentRules(void) const;                   //!< Get argument rules
        
        
    protected:
        
        
    private:
        
        const ArgumentRules*                                            argument_rules;                                  //!< Argument rules (different for different member functions)
        std::string                                                     method_name;                                       //!< Name of member method
        RevPtr<RevVariable>                                             object;
        const memberObjectType*                                         the_member_object;
    };
    
}


#include "MemberFunction.h"
#include "RbException.h"

/** default constructor */
template <typename memberObjectType, typename retType>
RevLanguage::MemberFunction<memberObjectType, retType>::MemberFunction( const std::string &name, const memberObjectType* o, ArgumentRules* ar ) : TypedFunction<retType>(),
    argument_rules( ar ),
    method_name( name ),
    object( NULL ),
    the_member_object( o )
{
    
}


/** copy constructor */
template <typename memberObjectType, typename retType>
RevLanguage::MemberFunction<memberObjectType, retType>::MemberFunction( const MemberFunction<memberObjectType, retType> &m) : TypedFunction<retType>(m),
    argument_rules( new ArgumentRules(*m.argument_rules) ),
    method_name( m.method_name ),
    object( m.object ),
    the_member_object( m.the_member_object )
{
    
}


/** default constructor */
template <typename memberObjectType, typename retType>
RevLanguage::MemberFunction<memberObjectType, retType>::~MemberFunction( void )
{
    delete argument_rules;
}


/** copy constructor */
template <typename memberObjectType, typename retType>
RevLanguage::MemberFunction<memberObjectType, retType>& RevLanguage::MemberFunction<memberObjectType, retType>::operator=( const MemberFunction<memberObjectType, retType> &m)
{

    // check for self-assignment
    if ( this != &m )
    {
        TypedFunction<retType>::operator=(m);
        
        delete argument_rules;
        
        argument_rules      = new ArgumentRules(*m.argument_rules);
        method_name         = m.method_name;
        object              = m.object;
        the_member_object   = m.the_member_object;

    }
    
    return *this;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename memberObjectType, typename retType>
RevLanguage::MemberFunction<memberObjectType, retType>* RevLanguage::MemberFunction<memberObjectType, retType>::clone( void ) const
{
    
    return new MemberFunction<memberObjectType, retType>( *this );
}


template <typename memberObjectType, typename retType>
RevBayesCore::TypedFunction< typename retType::valueType >* RevLanguage::MemberFunction<memberObjectType, retType>::createFunction( void ) const
{
        
    std::vector<const RevBayesCore::DagNode*> argNodes;
    for ( size_t i=0; i<this->args.size(); ++i)
    {
        argNodes.push_back( this->args[i].getVariable()->getRevObject().getDagNode() );
    }
    
    const RevBayesCore::TypedDagNode<typename memberObjectType::valueType>* o = the_member_object->getDagNode();
    if ( o == NULL )
    {
        throw RbException("Could not cast the member object.");
    }
    
    RevBayesCore::MemberFunction<typename memberObjectType::valueType, typename retType::valueType> *func = new RevBayesCore::MemberFunction<typename memberObjectType::valueType, typename retType::valueType>(this->method_name, o, argNodes);
    
    return func;
}


/** Get argument rules */
template <typename memberObjectType, typename retType>
const RevLanguage::ArgumentRules& RevLanguage::MemberFunction<memberObjectType, retType>::getArgumentRules(void) const
{
    
    return *argument_rules;
}


template <typename memberObjectType, typename retType>
const std::string& RevLanguage::MemberFunction<memberObjectType, retType>::getClassType(void)
{
    
    static std::string revClassType = "MemberFunction<" + memberObjectType::getClassType() + "," + retType::getClassType() + ">";
    
	return revClassType;
}

/* Get class type spec describing type of object */
template <typename memberObjectType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::MemberFunction<memberObjectType, retType>::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return revClassTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename memberObjectType, typename retType>
std::string RevLanguage::MemberFunction<memberObjectType, retType>::getFunctionName( void ) const
{
    
    return method_name;
}


template <typename memberObjectType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::MemberFunction<memberObjectType, retType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



template <typename memberObjectType, typename retType>
void RevLanguage::MemberFunction<memberObjectType, retType>::setMemberObject( const RevPtr<RevVariable> &obj) {
    
    // we do not own the object itself because one object can have multiple member functions
    object = obj;
    the_member_object = static_cast< const memberObjectType *>( &(obj->getRevObject()) );
}


#endif

