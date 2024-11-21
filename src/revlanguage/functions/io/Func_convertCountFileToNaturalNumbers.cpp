#include <sstream>
#include <cstdint>
#include <vector>

#include "ArgumentRule.h"
#include "CountFileToNaturalNumbersConverter.h"
#include "Func_convertCountFileToNaturalNumbers.h"
#include "RbException.h"
#include "RlString.h"
#include "StringUtilities.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "Integer.h"
#include "Natural.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_convertCountFileToNaturalNumbers* Func_convertCountFileToNaturalNumbers::clone( void ) const
{

	return new Func_convertCountFileToNaturalNumbers( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_convertCountFileToNaturalNumbers::execute( void )
{

	// get the information from the arguments for reading the file
	const RlString& fi = static_cast<const RlString&>( args[0].getVariable()->getRevObject() );
    RevBayesCore::TypedDagNode<std::int64_t>* n_individuals = static_cast<const Integer &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
	const RlString& fo = static_cast<const RlString&>( args[2].getVariable()->getRevObject() );
   
    RevBayesCore::CountFileToNaturalNumbersConverter nn;
    nn.cfconverter( fi.getValue(), n_individuals->getValue(), fo.getValue() );

	return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_convertCountFileToNaturalNumbers::getArgumentRules( void ) const
{

	static ArgumentRules argument_rules = ArgumentRules();
	static bool rules_set = false;

	if ( rules_set == false )
	{
		argument_rules.push_back( new ArgumentRule( "count file", RlString::getClassTypeSpec(), "A count file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "number of individuals", Natural::getClassTypeSpec(), "The number of (virtual or effective) individuals.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
		argument_rules.push_back( new ArgumentRule( "output file", RlString::getClassTypeSpec(), "A NaturalNumbers-type file, where each number corresponds to a PoMo state.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

		rules_set = true;

	}

	return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_convertCountFileToNaturalNumbers::getClassType(void)
{

	static std::string rev_type = "Func_convertCountFileToNaturalNumbers";

	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_convertCountFileToNaturalNumbers::getClassTypeSpec(void)
{

	static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_convertCountFileToNaturalNumbers::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "convertCountFileToNaturalNumbers";

    return f_name;
}


/** Get type spec */
const TypeSpec& Func_convertCountFileToNaturalNumbers::getTypeSpec( void ) const
{

	static TypeSpec type_spec = getClassTypeSpec();

	return type_spec;
}


/** Get return type */
const TypeSpec& Func_convertCountFileToNaturalNumbers::getReturnType( void ) const
{

	static TypeSpec return_typeSpec = RevNullObject::getClassTypeSpec();
	return return_typeSpec;
}
