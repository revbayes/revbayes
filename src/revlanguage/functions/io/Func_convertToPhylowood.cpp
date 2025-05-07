#include "Func_convertToPhylowood.h"

#include <cstddef>
#include <vector>

#include "ArgumentRule.h"
#include "OptionRule.h"
#include "PhylowoodConverter.h"
#include "Probability.h"
#include "RlString.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTree.h"
#include "TypeSpec.h"

using namespace RevLanguage;


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_convertToPhylowood* Func_convertToPhylowood::clone( void ) const
{
    
    return new Func_convertToPhylowood( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_convertToPhylowood::execute( void )
{
    
    
    const std::string& stateFilename = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string& treeFilename  = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    const std::string& geoFilename   = static_cast<const RlString&>( args[2].getVariable()->getRevObject() ).getValue();
    const std::string& phwFilename   = static_cast<const RlString&>( args[3].getVariable()->getRevObject() ).getValue();
    double burn = static_cast<const Probability&>( args[4].getVariable()->getRevObject() ).getValue();
    const std::string& chartype      = static_cast<const RlString&>( args[5].getVariable()->getRevObject() ).getValue();
    const std::string& bgtype      = static_cast<const RlString&>( args[6].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::PhylowoodConverter phw(stateFilename, treeFilename, geoFilename, phwFilename, burn, chartype, bgtype);
    //    RevBayesCore::PhylowoodConverter phw;
    
    return NULL;
}



/** Get argument rules */
const ArgumentRules& Func_convertToPhylowood::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        argumentRules.push_back( new ArgumentRule( "statefile", RlString::getClassTypeSpec()           , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "treefile",  RlString::getClassTypeSpec()           , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "geofile",   RlString::getClassTypeSpec()           , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "outfile",   RlString::getClassTypeSpec()           , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "burnin",    Probability::getClassTypeSpec()        , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );
        
        std::vector<std::string> options;
        options.push_back( "NaturalNumbers" );
        options.push_back( "Standard" );
        argumentRules.push_back( new OptionRule( "chartype", new RlString("NaturalNumbers"), options, "" ) );
        
        std::vector<std::string> options2;
        options2.push_back( "Range" );
        options2.push_back( "Area" );
        argumentRules.push_back( new OptionRule( "bgtype", new RlString("Area"), options2, "" ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_convertToPhylowood::getClassType(void)
{
    
    static std::string rev_type = "Func_convertToPhylowood";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_convertToPhylowood::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_convertToPhylowood::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "convertToPhylowood";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_convertToPhylowood::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_convertToPhylowood::getReturnType( void ) const {
    
    static TypeSpec return_typeSpec = Tree::getClassTypeSpec();
    return return_typeSpec;
}
