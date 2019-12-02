#ifndef Func_workspaceVector_H
#define Func_workspaceVector_H

#include "Procedure.h"

#include <string>

namespace RevLanguage {
    
    /**
     * @brief Func_workspaceVector: function creating workspace vectors
     *
     * This templated function constructs vectors and is used for language
     * constructs such as "v( x1, x2, ..., xn)" and "[ x1, x2, ..., xn ]" when
     * the elements are workspace objects.
     */
    template <typename valType>
    class Func_workspaceVector : public Procedure {
        
    public:
        Func_workspaceVector(void);                                                                 //!< Default constructor
        
        // Basic utility functions
        Func_workspaceVector*           clone(void) const;                                          //!< Clone the object
        static const std::string&       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                 getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        const ArgumentRules&            getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                 getReturnType(void) const;                                  //!< Get type of return value
        
        // Execute function
        RevPtr<RevVariable>             execute(void);                                              //!< Execute function
        
    };
    
}


#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "RbUtil.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"


/** Default constructor */
template <typename valType>
RevLanguage::Func_workspaceVector<valType>::Func_workspaceVector() : Procedure()
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename valType>
RevLanguage::Func_workspaceVector<valType>* RevLanguage::Func_workspaceVector<valType>::clone( void ) const
{
    return new Func_workspaceVector( *this );
}


/** Execute function: assemble a workspace vector */
template <typename valType>
RevLanguage::RevPtr<RevLanguage::RevVariable> RevLanguage::Func_workspaceVector<valType>::execute( void )
{
    RevBayesCore::RbVector<valType> params;
    for ( size_t i = 0; i < args.size(); i++ )
    {
        const valType &val = static_cast<const valType &>( args[i].getVariable()->getRevObject() );
        params.push_back( val );
    }
    
    WorkspaceVector<valType> *theVector = new WorkspaceVector<valType>( params );
    
    return new RevVariable( theVector );
}


/** Get argument rules */
template <typename valType>
const RevLanguage::ArgumentRules& RevLanguage::Func_workspaceVector<valType>::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "", valType::getClassTypeSpec(), "first value", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new Ellipsis ( "more values", valType::getClassTypeSpec() ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object (static version) */
template <typename valType>
const std::string& RevLanguage::Func_workspaceVector<valType>::getClassType( void )
{
    static std::string rev_type = "Func_workspaceVector<" + valType::getClassType() + ">";
    
	return rev_type;
}


/** Get Rev type spec of object (static version) */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_workspaceVector<valType>::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename valType>
std::string RevLanguage::Func_workspaceVector<valType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "v";
    
    return f_name;
}


/** Get Rev type spec of object (dynamic version) */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_workspaceVector<valType>::getTypeSpec( void ) const
{
    return this->getClassTypeSpec();
}


/** Get return type */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_workspaceVector<valType>::getReturnType( void ) const
{
    return WorkspaceVector<valType>::getClassTypeSpec();
}


#endif

