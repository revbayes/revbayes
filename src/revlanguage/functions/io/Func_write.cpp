#include <cstddef>
#include <string>
#include <iostream>
#include <vector>

#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "Func_write.h"
#include "RbFileManager.h"
#include "RevNullObject.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RbBoolean.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

#if defined (RB_MPI)
#include <mpi.h>
#endif

using namespace RevLanguage;


Func_write::Func_write( void ) :
process_ID( 0 )
{
#if defined (RB_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &process_ID);
#endif
}

using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_write* Func_write::clone( void ) const
{
    
    return new Func_write( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_write::execute( void )
{
    
    // get the information from the arguments for reading the file
    RevBayesCore::path fn = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    bool  append = static_cast<const RlBoolean&>( args[2].getVariable()->getRevObject() ).getValue();
    const std::string& separator = static_cast<const RlString&>( args[3].getVariable()->getRevObject() ).getValue();
    
    if ( process_ID == 0 )
    {
        
        if ( fn != "" )
        {
            RevBayesCore::createDirectoryForFile(fn);
            
            std::ofstream out_stream;
            
            if ( append == true )
            {
                
                // open the stream to the file
                out_stream.open( fn.string(), std::fstream::out | std::fstream::app);
            }
            else
            {
                
                // open the stream to the file
                out_stream.open( fn.string(), std::fstream::out);
            }
            
            // print the arguments
            args[0].getVariable()->getRevObject().printValue(out_stream, !false);
            for (size_t i = 4; i < args.size(); i++)
            {
                out_stream << separator;
                args[i].getVariable()->getRevObject().printValue( out_stream , !false );
            }
            
            out_stream.flush();
            out_stream.close();
        }
        else
        {
            
            std::ostream& o = std::cout;
            
            // print the arguments
            args[0].getVariable()->getRevObject().printValue( o, false );
            for (size_t i = 4; i < args.size(); i++)
            {
                o << separator;
                args[i].getVariable()->getRevObject().printValue( o, false );
            }
            o << std::endl;
        }
    }
    
#ifdef RB_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_write::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "", RevObject::getClassTypeSpec(), "A variable to write.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new Ellipsis( "Additional variables to write.", RevObject::getClassTypeSpec() ) );
        argumentRules.push_back( new ArgumentRule( "filename" , RlString::getClassTypeSpec() , "Output will be written to this file, or to the screen if no value is provided.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("") ) );
        argumentRules.push_back( new ArgumentRule( "append"   , RlBoolean::getClassTypeSpec(), "Append or overwrite existing file?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        argumentRules.push_back( new ArgumentRule( "separator", RlString::getClassTypeSpec() , "Character(s) with which to separate written variables.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("\\t (tab)") ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_write::getClassType(void)
{
    
    static std::string rev_type = "Func_write";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_write::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_write::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "print" );
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_write::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "write";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_write::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_write::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RevNullObject::getClassTypeSpec();
    return return_typeSpec;
}




