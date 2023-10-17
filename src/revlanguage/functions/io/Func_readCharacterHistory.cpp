#include <math.h>
#include <stddef.h>
#include <sstream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "CharacterHistoryDiscrete.h"
#include "Func_readCharacterHistory.h"
#include "ModelVector.h"
#include "NewickConverter.h"
#include "RlCharacterHistory.h"
#include "RLString.h"
#include "RlTree.h"
#include "WorkspaceToCoreWrapperObject.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readCharacterHistory* Func_readCharacterHistory::clone( void ) const
{
    
    return new Func_readCharacterHistory( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readCharacterHistory::execute( void )
{
    
    size_t arg_index = 0;
    
    // get the filename for the tree with MAP character history
    RevBayesCore::path simmap_filename = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::NewickConverter con;
    
    /* Open file */
    std::ifstream in_file( simmap_filename.string() );
    
    if ( !in_file )
        throw RbException()<<"Could not open file "<< simmap_filename;
    
    /* Initialize */
    RevBayesCore::RbVector<RevBayesCore::CharacterHistoryDiscrete> histories;
    
    /* line-processing loop */
    while ( in_file.good() )
    {
        
        // Read a line
        std::string line;
        RevBayesCore::safeGetline( in_file, line );
        
        // skip empty lines
        if (line.length() == 0)
        {
            continue;
        }
        
        RevBayesCore::CharacterHistoryDiscrete *simmap_tree = con.convertSimmapFromNewick( line );
        histories.push_back( *simmap_tree );
    }
    
    // return the tree with annotations
    return new RevVariable( new ModelVector<CharacterHistory>(  ) );
}



/** Get argument rules */
const ArgumentRules& Func_readCharacterHistory::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "The name of the file with the character history.",        ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_readCharacterHistory::getClassType(void)
{
    
    static std::string rev_type = "Func_readCharacterHistory";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readCharacterHistory::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readCharacterHistory::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readCharacterHistory";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readCharacterHistory::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readCharacterHistory::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = ModelVector<CharacterHistory>::getClassTypeSpec();
    return return_typeSpec;
}

