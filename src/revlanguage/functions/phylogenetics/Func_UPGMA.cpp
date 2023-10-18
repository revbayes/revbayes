#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Func_UPGMA.h"
#include "Procedure.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "UPGMA.h"

using namespace RevLanguage;

/** Default constructor */
Func_UPGMA::Func_UPGMA( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_UPGMA* Func_UPGMA::clone( void ) const
{
    
    return new Func_UPGMA( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_UPGMA::execute( void )
{
    const AbstractHomologousDiscreteCharacterData& char_data = static_cast<const AbstractHomologousDiscreteCharacterData &>( args[0].getVariable()->getRevObject() );
    
    bool exclude_ambiguous = false;
    RevBayesCore::DistanceMatrix dist_matrix = char_data.getValue().getPairwiseSequenceDifference( exclude_ambiguous );
    
    RevBayesCore::UPGMA upgma;
    RevBayesCore::Tree* upgma_tree = upgma.constructTree( dist_matrix );
    
    return new RevVariable( new TimeTree( upgma_tree ) );
}


/** Get argument rules */
const ArgumentRules& Func_UPGMA::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_UPGMA::getClassType(void)
{
    
    static std::string rev_type = "Func_UPGMA";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_UPGMA::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_UPGMA::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "UPGMA";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_UPGMA::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_UPGMA::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = TimeTree::getClassTypeSpec();
    
    return return_typeSpec;
}

