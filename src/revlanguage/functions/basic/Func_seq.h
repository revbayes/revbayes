#ifndef Func_seq_H
#define Func_seq_H

#include "Procedure.h"
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Function that creates a numerical sequence.
     *
     * This function is the equivalent to the 'R' function seq().
     * The function creates a numerical sequence from the specified value
     * to the specified value incremented by a specified value:
     * 
     * a <- seq(from=1,to=10,by=2)
     *
     * Which gives:
     * {1,3,5,7,9}
     *
     * This function is very similar to the range function except that
     * it has a user specified increment (or decrement).
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-06-19
     *
     */
    
    template <typename valType>
    class Func_seq : public Procedure {
        
    public:
        Func_seq();
        
        // Basic utility functions
        Func_seq*                                       clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                                 getReturnType(void) const;                                  //!< Get type of return value
        
        
        RevPtr<RevVariable>                             execute(void);                                              //!< Execute function
        
    protected:
        
        std::vector<std::string>                        getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                                     getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                                     getHelpDetails(void) const;                                 //!< Get the more detailed description of the function
        std::string                                     getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>      getHelpReferences(void) const;                              //!< Get some references/citations for this function
        std::vector<std::string>                        getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
        std::string                                     getHelpTitle(void) const;                                   //!< Get the title of this help entry

    };
    
}


#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "RbUtil.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"



template <typename valType>
RevLanguage::Func_seq<valType>::Func_seq() : Procedure() {
    
}

/* Clone object */
template <typename valType>
RevLanguage::Func_seq<valType>* RevLanguage::Func_seq<valType>::clone( void ) const {
    
    return new Func_seq( *this );
}


/** Execute function: We rely on getValue and overloaded push_back to provide functionality */
template <typename valType>
RevLanguage::RevPtr<RevLanguage::RevVariable> RevLanguage::Func_seq<valType>::execute( void )
{
    
    typename valType::valueType from = static_cast<const valType &>( args[0].getVariable()->getRevObject() ).getValue();
    typename valType::valueType to   = static_cast<const valType &>( args[1].getVariable()->getRevObject() ).getValue();
    typename valType::valueType by   = static_cast<const valType &>( args[2].getVariable()->getRevObject() ).getValue();
    
    typename valType::valueType val = from;
    
//    typename valType::valueType eps = valType::EPSILON;
    typename valType::valueType eps = std::numeric_limits<typename valType::valueType>::epsilon();
    
    ModelVector<valType> *seq = new ModelVector<valType>();
//    while ( (val >= from && val <= to) || (val <= from && val >= to) )
    while ( ((val - from) >= -eps && (val - to) <= eps) || ((val - from) <= eps && (val - to) >= -eps) )
    {
        seq->push_back( valType( val ) );
        val += by;
    }
    
    return new RevVariable( seq );
}


/** Get argument rules */
template <typename valType>
const RevLanguage::ArgumentRules& RevLanguage::Func_seq<valType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "from", valType::getClassTypeSpec(), "The first value of the sequence.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "to"  , valType::getClassTypeSpec(), "The last value of the sequence.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "by"  , valType::getClassTypeSpec(), "The step-size between value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Func_seq<valType>::getClassType(void)
{
    
    static std::string rev_type = "Func_seq<" + valType::getClassType() + ">";
    
	return rev_type;
}


/** Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_seq<valType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename valType>
std::string RevLanguage::Func_seq<valType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "seq";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
template <typename valType>
std::vector<std::string> RevLanguage::Func_seq<valType>::getHelpAuthor(void) const
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
std::string RevLanguage::Func_seq<valType>::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Create a sequence of values separate by a given step-size.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
template <typename valType>
std::string RevLanguage::Func_seq<valType>::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details =  "The 'seq' function create a sequence of values, starting with the initial value and then adding the step-size to it until the value reaches the 'to'-value." ;
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
template <typename valType>
std::string RevLanguage::Func_seq<valType>::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    example += "seq(-0.5, 10.5, 2)\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
template <typename valType>
std::vector<RevBayesCore::RbHelpReference> RevLanguage::Func_seq<valType>::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
template <typename valType>
std::vector<std::string> RevLanguage::Func_seq<valType>::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "rep" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
template <typename valType>
std::string RevLanguage::Func_seq<valType>::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Create a sequence values";
    
    return title;
}


/** Get type spec */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_seq<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_seq<valType>::getReturnType( void ) const {
    
    return ModelVector<valType>::getClassTypeSpec();
}


#endif


