#ifndef Func_BSPInterval_H
#define Func_BSPInterval_H

#include "ModelVector.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * @brief Func_BSPInterval: function creating model vectors
     *
     * This templated function constructs vectors and is used for language
     * constructs such as "v( x1, x2, ..., xn)" and "[ x1, x2, ..., xn ]" when
     * the elements are non-abstract model objects with non-abstract value types.
     */
    template <typename valType>
    class Func_BSPInterval : public TypedFunction< ModelVector< valType > > {
        
    public:
        Func_BSPInterval(void);                                                                 //!< Default constructor
        
        // Basic utility functions
        Func_BSPInterval*                                                                                 clone(void) const;                                          //!< Clone the object
        static const std::string&                                                                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                                                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                                                     getFunctionName(void) const;
        std::vector<std::string>                                                                        getFunctionNameAliases( void ) const;
        const TypeSpec&                                                                                 getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::RbVector<typename valType::valueType > >*            createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                                            getArgumentRules(void) const;                               //!< Get argument rules
        
    protected:

    };
    
}


#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "BSPIntervalFunction.h"
#include "RbUtil.h"
#include "RbVector.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "VectorFunction.h"


/** Default constructor */
template <typename valType>
RevLanguage::Func_BSPInterval<valType>::Func_BSPInterval() : TypedFunction< ModelVector<valType> >()
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename valType>
RevLanguage::Func_BSPInterval<valType>* RevLanguage::Func_BSPInterval<valType>::clone( void ) const
{
    return new Func_BSPInterval( *this );
}


/** Execute function: create deterministic replicate<valType> object */
template <typename valType>
RevBayesCore::TypedFunction< RevBayesCore::RbVector< typename valType::valueType> >* RevLanguage::Func_BSPInterval<valType>::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<typename valType::valueType> >* v  = static_cast<const ModelVector< valType > &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<long> >* n                         = static_cast<const ModelVector< Natural > &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::BSPIntervalFunction<typename valType::valueType>* func = new RevBayesCore::BSPIntervalFunction<typename valType::valueType>( v, n );
    
    return func;
}


/** Get argument rules */
template <typename valType>
const RevLanguage::ArgumentRules& RevLanguage::Func_BSPInterval<valType>::getArgumentRules( void ) const
{
    static ArgumentRules argument_rules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "x", ModelVector< valType >::getClassTypeSpec(), "The values that we replicate.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "n", ModelVector< Natural >::getClassTypeSpec(), "How often we replicate each value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argument_rules;
}


/** Get Rev type of object (static version) */
template <typename valType>
const std::string& RevLanguage::Func_BSPInterval<valType>::getClassType( void )
{
    static std::string rev_type = "Func_BSPInterval<" + valType::getClassType() + ">";
    
    return rev_type;
}


/** Get Rev type spec of object (static version) */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_BSPInterval<valType>::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename valType>
std::string RevLanguage::Func_BSPInterval<valType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "BSPInterval";
    
    return f_name;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
template <typename valType>
std::vector<std::string> RevLanguage::Func_BSPInterval<valType>::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    
    return a_names;
}


/** Get Rev type spec of object (dynamic version) */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_BSPInterval<valType>::getTypeSpec( void ) const
{
    return this->getClassTypeSpec();
}


#endif

