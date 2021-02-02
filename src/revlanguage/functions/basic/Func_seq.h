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


