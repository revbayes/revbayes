#ifndef Func_replicate_H
#define Func_replicate_H

#include "ModelVector.h"
#include "RlTypedFunction.h"

#include <cmath>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Func_replicate: function creating model vectors
     *
     * This templated function constructs vectors and is used for language
     * constructs such as "v( x1, x2, ..., xn)" and "[ x1, x2, ..., xn ]" when
     * the elements are non-abstract model objects with non-abstract value types.
     */
    template <typename valType, typename nonNegType>
    class Func_replicate : public TypedFunction< ModelVector< valType > > {
        
    public:
        Func_replicate(void);                                                                 //!< Default constructor
        
        // Basic utility functions
        Func_replicate*                                                                                 clone(void) const;                                          //!< Clone the object
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
#include "ReplicateFunction.h"
#include "RbException.h"
#include "RbUtil.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "VectorFunction.h"


/** Default constructor */
template <typename valType, typename nonNegType>
RevLanguage::Func_replicate<valType, nonNegType>::Func_replicate() : TypedFunction< ModelVector<valType> >()
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename valType, typename nonNegType>
RevLanguage::Func_replicate<valType, nonNegType>* RevLanguage::Func_replicate<valType, nonNegType>::clone( void ) const
{
    return new Func_replicate( *this );
}


/** Execute function: create deterministic replicate<valType> object */
template <typename valType, typename nonNegType>
RevBayesCore::TypedFunction< RevBayesCore::RbVector< typename valType::valueType> >* RevLanguage::Func_replicate<valType, nonNegType>::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode<typename valType::valueType>* v   = static_cast<const valType &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    double n_initial = static_cast<const nonNegType &>( this->args[1].getVariable()->getRevObject() ).getValue();
    int n = static_cast<int>(std::floor(n_initial));
    if (n != n_initial) {
        // Real or Integer types are permitted, but the value must be an integer
        throw RbException("n must have an integer value; try floor(n) or ceil(n)");
    }
    
    // Checking for negatives allows a user to pass a Real without first 
    // coercing into RealPos
    if (n < 0) {
         throw RbException("n may not be negative");
    }
    
    auto* func = new RevBayesCore::ReplicateFunction<typename valType::valueType>( v, n );
    
    return func;
}


/** Get argument rules */
template <typename valType, typename nonNegType>
const RevLanguage::ArgumentRules& RevLanguage::Func_replicate<valType, nonNegType>::getArgumentRules( void ) const
{
    static ArgumentRules argument_rules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "x", valType::getClassTypeSpec(), "The value that we replicate.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "n", nonNegType::getClassTypeSpec(), "How often we replicate the value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argument_rules;
}


/** Get Rev type of object (static version) */
template <typename valType, typename nonNegType>
const std::string& RevLanguage::Func_replicate<valType, nonNegType>::getClassType( void )
{
    static std::string rev_type = "Func_replicate<" + valType::getClassType() + ", " + nonNegType::getClassType() + ">";
    
    return rev_type;
}


/** Get Rev type spec of object (static version) */
template <typename valType, typename nonNegType>
const RevLanguage::TypeSpec& RevLanguage::Func_replicate<valType, nonNegType>::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename valType, typename nonNegType>
std::string RevLanguage::Func_replicate<valType, nonNegType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "rep";
    
    return f_name;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
template <typename valType, typename nonNegType>
std::vector<std::string> RevLanguage::Func_replicate<valType, nonNegType>::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "replicate" );
    
    return a_names;
}


/** Get Rev type spec of object (dynamic version) */
template <typename valType, typename nonNegType>
const RevLanguage::TypeSpec& RevLanguage::Func_replicate<valType, nonNegType>::getTypeSpec( void ) const
{
    return this->getClassTypeSpec();
}


#endif

