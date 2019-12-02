#ifndef Func_replicate_H
#define Func_replicate_H

#include "ModelVector.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * @brief Func_replicate: function creating model vectors
     *
     * This templated function constructs vectors and is used for language
     * constructs such as "v( x1, x2, ..., xn)" and "[ x1, x2, ..., xn ]" when
     * the elements are non-abstract model objects with non-abstract value types.
     */
    template <typename valType>
    class Func_replicate : public TypedFunction< ModelVector< valType> > {
        
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
        
        std::vector<std::string>                                                                        getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                                                                                     getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                                                                                     getHelpDetails(void) const;                                 //!< Get the more detailed description of the function
        std::string                                                                                     getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>                                                      getHelpReferences(void) const;                              //!< Get some references/citations for this function
        std::vector<std::string>                                                                        getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
        std::string                                                                                     getHelpTitle(void) const;                                   //!< Get the title of this help entry
        
    };
    
}


#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "ReplicateFunction.h"
#include "RbUtil.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "VectorFunction.h"


/** Default constructor */
template <typename valType>
RevLanguage::Func_replicate<valType>::Func_replicate() : TypedFunction< ModelVector<valType> >()
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename valType>
RevLanguage::Func_replicate<valType>* RevLanguage::Func_replicate<valType>::clone( void ) const
{
    return new Func_replicate( *this );
}


/** Execute function: create deterministic replicate<valType> object */
template <typename valType>
RevBayesCore::TypedFunction< RevBayesCore::RbVector< typename valType::valueType> >* RevLanguage::Func_replicate<valType>::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode<typename valType::valueType>* v   = static_cast<const valType &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    int n   = (int)static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::ReplicateFunction<typename valType::valueType>* func = new RevBayesCore::ReplicateFunction<typename valType::valueType>( v, n );
    
    return func;
}


/** Get argument rules */
template <typename valType>
const RevLanguage::ArgumentRules& RevLanguage::Func_replicate<valType>::getArgumentRules( void ) const
{
    static ArgumentRules argument_rules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "x", valType::getClassTypeSpec(), "The value that we replicate.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "n", Natural::getClassTypeSpec(), "How often we replicate the value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argument_rules;
}


/** Get Rev type of object (static version) */
template <typename valType>
const std::string& RevLanguage::Func_replicate<valType>::getClassType( void )
{
    static std::string rev_type = "Func_replicate<" + valType::getClassType() + ">";
    
    return rev_type;
}


/** Get Rev type spec of object (static version) */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_replicate<valType>::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename valType>
std::string RevLanguage::Func_replicate<valType>::getFunctionName( void ) const
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
template <typename valType>
std::vector<std::string> RevLanguage::Func_replicate<valType>::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "replicate" );
    
    return a_names;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
template <typename valType>
std::vector<std::string> RevLanguage::Func_replicate<valType>::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
template <typename valType>
std::string RevLanguage::Func_replicate<valType>::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "'rep' creates a vector of 'n' copies of the value 'x'.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
template <typename valType>
std::string RevLanguage::Func_replicate<valType>::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details =  "'rep' creates a vector of 'n' elements, each with value 'x', preserving the type of 'x' in the returned vector." ;
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
template <typename valType>
std::string RevLanguage::Func_replicate<valType>::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "rep(0.1, 3)\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
template <typename valType>
std::vector<RevBayesCore::RbHelpReference> RevLanguage::Func_replicate<valType>::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
template <typename valType>
std::vector<std::string> RevLanguage::Func_replicate<valType>::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "simplex" );
    see_also.push_back( "v" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
template <typename valType>
std::string RevLanguage::Func_replicate<valType>::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Replicate a value";
    
    return title;
}


/** Get Rev type spec of object (dynamic version) */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_replicate<valType>::getTypeSpec( void ) const
{
    return this->getClassTypeSpec();
}


#endif

