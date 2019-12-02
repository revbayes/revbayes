#include "HomologousDiscreteCharacterData.h"
#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "PoMoCountFileReader.h"
#include "Ellipsis.h"
#include "Func_readPoMoCountFile.h"
#include "ModelVector.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RevNullObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBoolean.h"
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
Func_readPoMoCountFile* Func_readPoMoCountFile::clone( void ) const
{

	return new Func_readPoMoCountFile( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readPoMoCountFile::execute( void )
{

	// get the information from the arguments for reading the file
	const RlString& fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() );
	RevBayesCore::TypedDagNode<long>* virtualPopulationSize = static_cast<const Integer &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<long>* n_states = static_cast<const Integer &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

	RevBayesCore::PoMoCountFileReader* pcfr = new RevBayesCore::PoMoCountFileReader( fn.getValue(), virtualPopulationSize->getValue(), ' ' );

	RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::PoMoState> *pomoAln = new RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::PoMoState>( *(pcfr->getMatrix() ) );

	AbstractHomologousDiscreteCharacterData *rlPoMoAln = new AbstractHomologousDiscreteCharacterData( pomoAln );

	return new RevVariable( rlPoMoAln );
}


/** Get argument rules */
const ArgumentRules& Func_readPoMoCountFile::getArgumentRules( void ) const
{

	static ArgumentRules argument_rules = ArgumentRules();
	static bool rules_set = false;

	if ( rules_set == false )
	{
		argument_rules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "Relative or absolute name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
		argument_rules.push_back( new ArgumentRule( "virtualPopulationSize", Natural::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "numStates", Natural::getClassTypeSpec(), "The number of states (e.g. 4 for A,C,G and T).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(4) ) );

		rules_set = true;

	}

	return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_readPoMoCountFile::getClassType(void)
{

	static std::string rev_type = "Func_readPoMoCountFile";

	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readPoMoCountFile::getClassTypeSpec(void)
{

	static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readPoMoCountFile::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readPoMoCountFile";

    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readPoMoCountFile::getTypeSpec( void ) const
{

	static TypeSpec type_spec = getClassTypeSpec();

	return type_spec;
}


/** Get return type */
const TypeSpec& Func_readPoMoCountFile::getReturnType( void ) const
{

	static TypeSpec return_typeSpec = AbstractHomologousDiscreteCharacterData::getClassTypeSpec();
	return return_typeSpec;
}
