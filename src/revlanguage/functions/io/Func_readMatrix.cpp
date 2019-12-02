#include "AbstractHomologousDiscreteCharacterData.h"
#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "MatrixReal.h"
#include "MatrixReader.h"
#include "Ellipsis.h"
#include "Func_readMatrix.h"
#include "ModelVector.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RevNullObject.h"
#include "RlBoolean.h"
#include "RlMatrixReal.h"
#include "RlString.h"
#include "RlUtils.h"
#include "StringUtilities.h"
#include "RlUserInterface.h"

#include <map>
#include <set>
#include <sstream>


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readMatrix* Func_readMatrix::clone( void ) const
{
	
	return new Func_readMatrix( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readMatrix::execute( void )
{
	
	// get the information from the arguments for reading the file
	const RlString&    fn  = static_cast<const RlString&>( args[0].getVariable()->getRevObject() );
    const std::string& del = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();

    RevBayesCore::MatrixReader* mr = new RevBayesCore::MatrixReader( fn.getValue(), del[0] );
    RevBayesCore::MatrixReal    m  = mr->getMatrix();
    MatrixReal*                 rm = new MatrixReal( m );

    return new RevVariable( rm );
    
}


/** Get argument rules */
const ArgumentRules& Func_readMatrix::getArgumentRules( void ) const
{
	
	static ArgumentRules argumentRules = ArgumentRules();
	static bool rules_set = false;
	
	if (!rules_set)
	{
		
		argumentRules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "Relative or absolute name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "delimiter", RlString::getClassTypeSpec(), "The delimiter between columns.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "\t" ) ) );
		rules_set = true;
		
	}
	
	return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readMatrix::getClassType(void)
{
	
	static std::string rev_type = "Func_readMatrix";
	
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readMatrix::getClassTypeSpec(void)
{
	
	static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
	
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readMatrix";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readMatrix::getTypeSpec( void ) const
{
	
	static TypeSpec type_spec = getClassTypeSpec();
	
	return type_spec;
}


/** Get return type */
const TypeSpec& Func_readMatrix::getReturnType( void ) const
{
	
	static TypeSpec return_typeSpec = MatrixReal::getClassTypeSpec();
	return return_typeSpec;
}

