#ifndef RlDagMemberFunction_H
#define RlDagMemberFunction_H

#include "DagNode.h"
#include "RlMemberMethod.h"
#include "RlTypedFunction.h"

namespace RevLanguage {
    
    template <typename retType>
    class DagMemberFunction : public TypedFunction<retType>, public MemberMethod {
        
    public:
        DagMemberFunction(const std::string &n, const RevBayesCore::DagNode* o, ArgumentRules* argRules, bool f=true);                                 //!< Constructor
        virtual ~DagMemberFunction(){};
        
        // Basic utility functions
        DagMemberFunction*                                              clone(void) const;                              //!< Clone the object
        static const std::string&                                       getClassType(void);                             //!< Get class name
        static const TypeSpec&                                          getClassTypeSpec(void);                         //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                        //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< typename retType::valueType>*      createFunction(void) const;                     //!< Execute function
        void                                                            setMemberObject(const RevPtr<RevVariable> &obj);   //!< Set the member object to which this function belongs
        const ArgumentRules&                                            getArgumentRules(void) const;                   //!< Get argument rules
        void                                                            setMethodName(const std::string& name);         //!< Set name of member method
        
        
    private:
        
        const ArgumentRules*                                            argument_rules;                                  //!< Argument rules (different for different member functions)
        std::string                                                     func_name;                                       //!< Name of member method
        RevPtr<RevVariable>                                             object;
        const RevBayesCore::DagNode*                                    the_member_object;
        bool                                                            force_updates;
        
    };
    
    
    
}


#include "DagMemberFunction.h"
#include "RbException.h"

/** default constructor */
template <typename retType>
RevLanguage::DagMemberFunction<retType>::DagMemberFunction( const std::string &n, const RevBayesCore::DagNode *o, ArgumentRules* ar, bool f ) : TypedFunction<retType>(  ),
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
template <typename retType>
RevLanguage::DagMemberFunction<retType>* RevLanguage::DagMemberFunction< retType>::clone( void ) const
{
    
    return new DagMemberFunction< retType>( *this );
}


template <typename retType>
RevBayesCore::TypedFunction<typename retType::valueType >* RevLanguage::DagMemberFunction< retType>::createFunction( void ) const
{
    
    std::vector<const RevBayesCore::DagNode*> argNodes;
    for ( size_t i=0; i<this->args.size(); ++i)
    {
        argNodes.push_back( this->args[i].getVariable()->getRevObject().getDagNode() );
    }
    
    RevBayesCore::DagMemberFunction<typename retType::valueType> *func = new RevBayesCore::DagMemberFunction<typename retType::valueType>(this->func_name, the_member_object, argNodes);
    func->setForceUpdates( force_updates );
    
    return func;
}


/** Get argument rules */
template <typename retType>
const RevLanguage::ArgumentRules& RevLanguage::DagMemberFunction<retType>::getArgumentRules(void) const
{
    
    return *argument_rules;
}



template <typename retType>
const std::string& RevLanguage::DagMemberFunction<retType>::getClassType(void)
{
    
    static std::string revClassType = "DagMemberFunction<" + retType::getClassType() + ">";
    
    return revClassType;
}

/* Get class type spec describing type of object */
template <typename retType>
const RevLanguage::TypeSpec& RevLanguage::DagMemberFunction< retType>::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return revClassTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename retType>
std::string RevLanguage::DagMemberFunction< retType>::getFunctionName( void ) const
{
    
    return func_name;
}


template <typename retType>
const RevLanguage::TypeSpec& RevLanguage::DagMemberFunction< retType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


template <typename retType>
void RevLanguage::DagMemberFunction< retType>::setMemberObject( const RevPtr<RevVariable> &obj)
{
    
    // we do not own the object itself because one object can have multiple member functions
    object = obj;
    the_member_object = (obj->getRevObject().getDagNode());
}


#endif

