#include <sstream>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DiscreteTaxonData.h"
#include "Func_readPoMoCountFile.h"
#include "HomologousDiscreteCharacterData.h"
#include "Integer.h"
#include "Natural.h"
#include "NaturalNumbersState.h"
#include "OptionRule.h"
#include "PoMoCountFileReader.h"
#include "PoMoState.h"
#include "RbException.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"


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
    const std::string& file_in      = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    long  virtual_population_size   = static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getValue();

    const std::string& format       = static_cast<const RlString&>( args[2].getVariable()->getRevObject() ).getValue();
    RevBayesCore::PoMoCountFileReader::FORMAT reader_format = RevBayesCore::PoMoCountFileReader::PoMo;
    if ( format == "NaturalNumbers" )
    {
        reader_format = RevBayesCore::PoMoCountFileReader::NaturalNumbers;
    }
    
    const std::string& weighting_method = static_cast<const RlString&>(       args[3].getVariable()->getRevObject() ).getValue();
    long  effective_population_size     = static_cast<const Natural &>( this->args[4].getVariable()->getRevObject() ).getValue();

	RevBayesCore::PoMoCountFileReader* pcfr = new RevBayesCore::PoMoCountFileReader( file_in, virtual_population_size, reader_format, weighting_method, effective_population_size );

	AbstractHomologousDiscreteCharacterData *rlPoMoAln = new AbstractHomologousDiscreteCharacterData( pcfr->getMatrix() );

	return new RevVariable( rlPoMoAln );
}


/** Get argument rules */
const ArgumentRules& Func_readPoMoCountFile::getArgumentRules( void ) const
{

	static ArgumentRules argument_rules = ArgumentRules();
	static bool rules_set = false;

    if ( rules_set == false )
    {
        argument_rules.push_back( new ArgumentRule( "countFile", RlString::getClassTypeSpec(), "A count file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "virtualPopulationSize", Natural::getClassTypeSpec(), "The number of virtual individuals in the population.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
            
        std::vector<std::string> format_options;
        format_options.push_back( "PoMo" );
        format_options.push_back( "NaturalNumbers" );
        argument_rules.push_back( new OptionRule( "format", new RlString("PoMo"), format_options, "The output data type format." ) );
        
        std::vector<std::string> weighting_options;
        weighting_options.push_back( "Fixed" );
        weighting_options.push_back( "Binomial" );
        weighting_options.push_back( "Sampled" );
        weighting_options.push_back( "Hypergeometric" );
        weighting_options.push_back( "None" );
        argument_rules.push_back( new OptionRule( "samplingCorrection", new RlString("Fixed"), weighting_options, "The sampling correction to map the observed counts to PoMo states." ) );
        
        argument_rules.push_back( new ArgumentRule( "effectivePopulationSize", Natural::getClassTypeSpec(), "A tentative population size used by the hypergeometric method to correct the counts.", ArgumentRule::BY_VALUE, ArgumentRule::ANY,  new Natural(10000) ) );

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
