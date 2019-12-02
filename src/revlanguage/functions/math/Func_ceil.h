#ifndef Func_ceil_H
#define Func_ceil_H

#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the ceil function.
     *
     * The RevLanguage wrapper of the ceil function connects
     * the variables/parameters of the function and creates the internal CeilFunction object.
     * Please read the CeilFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    template <typename valType, typename retType>
    class Func_ceil : public TypedFunction<retType> {
        
    public:
        Func_ceil( void );
        
        // Basic utility functions
        Func_ceil*                                                          clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< typename retType::valueType>*          createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules

    protected:
        
        std::vector<std::string>                                            getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                                                         getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                                                         getHelpDetails(void) const;                                 //!< Get the more detailed description of the function
        std::string                                                         getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>                          getHelpReferences(void) const;                              //!< Get some references/citations for this function
        std::vector<std::string>                                            getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
        std::string                                                         getHelpTitle(void) const;                                   //!< Get the title of this help entry

    };
    
}


#include "CeilFunction.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

/** default constructor */
template <typename valType, typename retType>
RevLanguage::Func_ceil<valType, retType>::Func_ceil( void ) : TypedFunction<retType>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename valType, typename retType>
RevLanguage::Func_ceil<valType, retType>* RevLanguage::Func_ceil<valType, retType>::clone( void ) const
{
    
    return new Func_ceil<valType, retType>( *this );
}


template <typename valType, typename retType>
RevBayesCore::TypedFunction< typename retType::valueType >* RevLanguage::Func_ceil<valType, retType>::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* arg = static_cast<const valType &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::CeilFunction* f = new RevBayesCore::CeilFunction( arg );
    
    return f;
}


/* Get argument rules */
template <typename valType, typename retType>
const RevLanguage::ArgumentRules& RevLanguage::Func_ceil<valType, retType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set ) 
    {
        
        argumentRules.push_back( new ArgumentRule( "x", valType::getClassTypeSpec(), "The value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename valType, typename retType>
const std::string& RevLanguage::Func_ceil<valType, retType>::getClassType(void)
{
    
    static std::string rev_type = "Func_ceil";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
template <typename valType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func_ceil<valType, retType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( valType::getClassTypeSpec() ), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
template <typename valType, typename retType>
std::string RevLanguage::Func_ceil<valType, retType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "ceil";
    
    return f_name;
}



/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
template <typename valType, typename retType>
std::vector<std::string> RevLanguage::Func_ceil<valType, retType>::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
template <typename valType, typename retType>
std::string RevLanguage::Func_ceil<valType, retType>::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "The 'ceil' function maps the value of a number to the smallest following integer.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
template <typename valType, typename retType>
std::string RevLanguage::Func_ceil<valType, retType>::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details;
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
template <typename valType, typename retType>
std::string RevLanguage::Func_ceil<valType, retType>::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    example += "# compute the ceiling of a real number\n";
    example += "number <- 3.4\n";
    example += "ceiled_number <- ceil(number)\n";
    example += "if (ceiled_number != 4.0) {\n";
    example += "    print(\"Problem when computing a ceiled value.\")\n";
    example += "} else {\n";
    example += "    print(\"Correct computation of a ceiled value.\")\n";
    example += "}\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
template <typename valType, typename retType>
std::vector<RevBayesCore::RbHelpReference> RevLanguage::Func_ceil<valType, retType>::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
template <typename valType, typename retType>
std::vector<std::string> RevLanguage::Func_ceil<valType, retType>::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "abs" );
    see_also.push_back( "floor" );
    see_also.push_back( "round" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
template <typename valType, typename retType>
std::string RevLanguage::Func_ceil<valType, retType>::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Ceiling value of a number";
    
    return title;
}


template <typename valType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func_ceil<valType, retType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif
