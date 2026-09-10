#include <sstream>
#include <vector>
#include <string>

#include "ArgumentRule.h"
#include "FastaFileToNaturalNumbersConverter.h"
#include "Func_convertFastaFileToNaturalNumbers.h"
#include "RbException.h"
#include "RlString.h"
#include "StringUtilities.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
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
Func_convertFastaFileToNaturalNumbers* Func_convertFastaFileToNaturalNumbers::clone( void ) const
{

	return new Func_convertFastaFileToNaturalNumbers( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_convertFastaFileToNaturalNumbers::execute( void )
{

	// get the information from the arguments for reading the file
	const RlString& fi                                                         = static_cast<const RlString&>(                     args[0].getVariable()->getRevObject() );
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<std::string> >* taxa     = static_cast<const ModelVector<RlString> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<std::string> >* alleles  = static_cast<const ModelVector<RlString> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    //std::vector<std::string> taxa     = static_cast<const ModelVector<RlString> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    //std::vector<std::string> alleles  = static_cast<const ModelVector<RlString> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode<std::int64_t>*                           n_individuals  = static_cast<const Integer &>(               this->args[3].getVariable()->getRevObject() ).getDagNode();
	const RlString& fo                                                         = static_cast<const RlString&>(                     args[4].getVariable()->getRevObject() );
   
    RevBayesCore::FastaFileToNaturalNumbersConverter nn;
    nn.faconverter( fi.getValue(), taxa->getValue(), alleles->getValue(),  n_individuals->getValue(), fo.getValue() );

	return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_convertFastaFileToNaturalNumbers::getArgumentRules( void ) const
{

	static ArgumentRules argument_rules = ArgumentRules();
	static bool rules_set = false;

	if ( rules_set == false )
	{
		argument_rules.push_back( new ArgumentRule( "fasta file", RlString::getClassTypeSpec(), "A fasta file: sequence names must all start with one of the taxa names; this how they are assigned to each taxa. E.g., >taxa1plusotherinfo, >taxa1_plusotherinfo or just >taxa1.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "taxa"      , ModelVector<RlString>::getClassTypeSpec() , "A string vector listing the the taxa names.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "alleles"   , ModelVector<RlString>::getClassTypeSpec(), "A string vector listing the alleles.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "number of individuals", Natural::getClassTypeSpec(), "The number of (virtual or effective) individuals.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
		argument_rules.push_back( new ArgumentRule( "output file", RlString::getClassTypeSpec(), "A NaturalNumbers-type file, where each number corresponds to a PoMo state.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

		rules_set = true;

	}

	return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_convertFastaFileToNaturalNumbers::getClassType(void)
{

	static std::string rev_type = "Func_convertFastaFileToNaturalNumbers";

	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_convertFastaFileToNaturalNumbers::getClassTypeSpec(void)
{

	static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_convertFastaFileToNaturalNumbers::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "convertFastaFileToNaturalNumbers";

    return f_name;
}


/** Get type spec */
const TypeSpec& Func_convertFastaFileToNaturalNumbers::getTypeSpec( void ) const
{

	static TypeSpec type_spec = getClassTypeSpec();

	return type_spec;
}


/** Get return type */
const TypeSpec& Func_convertFastaFileToNaturalNumbers::getReturnType( void ) const
{

	static TypeSpec return_typeSpec = RevNullObject::getClassTypeSpec();
	return return_typeSpec;
}
