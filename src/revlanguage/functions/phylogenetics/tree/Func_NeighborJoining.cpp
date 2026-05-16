#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Func_NeighborJoining.h"
#include "Procedure.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBranchLengthTree.h"
#include "TypeSpec.h"
#include "NeighborJoining.h"

using namespace RevLanguage;

/** Default constructor */
Func_NeighborJoining::Func_NeighborJoining( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_NeighborJoining* Func_NeighborJoining::clone( void ) const
{
    
    return new Func_NeighborJoining( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_NeighborJoining::execute( void )
{
    const AbstractHomologousDiscreteCharacterData& char_data = static_cast<const AbstractHomologousDiscreteCharacterData &>( args[0].getVariable()->getRevObject() );
    
    bool exclude_ambiguous  = false;
    bool relative           = true; // normalize by number of sites
    RevBayesCore::DistanceMatrix dist_matrix = char_data.getValue().getPairwiseSequenceDifference( exclude_ambiguous, relative );
        
    RevBayesCore::NeighborJoining nj;
    RevBayesCore::Tree* nj_tree = nj.constructTree( dist_matrix );
    
    return new RevVariable( new BranchLengthTree( nj_tree ) );
}


/** Get argument rules */
const ArgumentRules& Func_NeighborJoining::getArgumentRules( void ) const
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
const std::string& Func_NeighborJoining::getClassType(void)
{
    
    static std::string rev_type = "Func_NeighborJoining";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_NeighborJoining::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_NeighborJoining::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "NeighborJoining";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_NeighborJoining::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_NeighborJoining::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = BranchLengthTree::getClassTypeSpec();
    
    return return_typeSpec;
}

