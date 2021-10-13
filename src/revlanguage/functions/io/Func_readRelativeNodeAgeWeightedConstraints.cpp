#include <sstream>
#include <vector>

#include "ArgumentRule.h"
#include "RelativeNodeAgeWeightedConstraints.h"
#include "RelativeNodeAgeWeightedConstraintsReader.h"
#include "Func_readRelativeNodeAgeWeightedConstraints.h"
#include "Real.h"
#include "RlRelativeNodeAgeWeightedConstraints.h"
#include "RlString.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/** Clone object */
Func_readRelativeNodeAgeWeightedConstraints* Func_readRelativeNodeAgeWeightedConstraints::clone( void ) const
{
    
    return new Func_readRelativeNodeAgeWeightedConstraints( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readRelativeNodeAgeWeightedConstraints::execute( void )
{
    
    // get the information from the arguments for reading the file
    std::string fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    std::string sep = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    double th = static_cast<const Real&>( args[2].getVariable()->getRevObject() ).getValue();

    RevBayesCore::RelativeNodeAgeWeightedConstraintsReader* dmr = new RevBayesCore::RelativeNodeAgeWeightedConstraintsReader( fn, sep, 0, th );
    RevBayesCore::RelativeNodeAgeWeightedConstraints* dm = new RevBayesCore::RelativeNodeAgeWeightedConstraints(dmr);
    
    return new RevVariable( new RlRelativeNodeAgeWeightedConstraints(dm) );
}


/** Get argument rules */
const ArgumentRules& Func_readRelativeNodeAgeWeightedConstraints::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "Relative or absolute name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        std::vector<std::string> sep = {"separator","delimiter"};
        argumentRules.push_back( new ArgumentRule( sep, RlString::getClassTypeSpec(), "The separator/delimiter between columns.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "" ) ) );
        argumentRules.push_back( new ArgumentRule( "threshold", Real::getClassTypeSpec(), "weight threshold below which constraints are ignored.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
        
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readRelativeNodeAgeWeightedConstraints::getClassType(void)
{
    
    static std::string rev_type = "Func_readRelativeNodeAgeWeightedConstraints";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readRelativeNodeAgeWeightedConstraints::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readRelativeNodeAgeWeightedConstraints::getFunctionName( void ) const
{
    // create a name variable that is the same for all instances of this class
    std::string f_name = "readRelativeNodeAgeWeightedConstraints";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readRelativeNodeAgeWeightedConstraints::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readRelativeNodeAgeWeightedConstraints::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlRelativeNodeAgeWeightedConstraints::getClassTypeSpec();
    return return_typeSpec;
}

