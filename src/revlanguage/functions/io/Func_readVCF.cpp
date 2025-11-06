#include "HomologousDiscreteCharacterData.h"
#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "VCFReader.h"
#include "Ellipsis.h"
#include "Func_readVCF.h"
#include "OptionRule.h"
#include "RbException.h"
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
Func_readVCF* Func_readVCF::clone( void ) const
{
    
    return new Func_readVCF( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readVCF::execute( void )
{
    size_t arg_index = 0;
    
    // get the information from the arguments for reading the file
    const RlString& fn = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() );
//    RevBayesCore::TypedDagNode<std::int64_t>* virtualPopulationSize = static_cast<const Integer &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
//    RevBayesCore::TypedDagNode<std::int64_t>* n_states = static_cast<const Integer &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    
    const std::string type = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();

    RevBayesCore::VCFReader vcf_reader = RevBayesCore::VCFReader( fn.getValue() );
    
    AbstractHomologousDiscreteCharacterData *rl_aln = NULL;
    
    if ( type == "DNA" )
    {
        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::DnaState> *aln = vcf_reader.readDNAMatrix();
        rl_aln = new AbstractHomologousDiscreteCharacterData( aln );
    }
    else if ( type == "binary")
    {
        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::BinaryState> *aln = vcf_reader.readBinaryMatrix();
        rl_aln = new AbstractHomologousDiscreteCharacterData( aln );
    }
    return new RevVariable( rl_aln );
}


/** Get argument rules */
const ArgumentRules& Func_readVCF::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        argument_rules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "Relative or absolute name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        std::vector<std::string> character_options;
        character_options.push_back( "DNA" );
        character_options.push_back( "binary" );
        argument_rules.push_back( new OptionRule( "type", new RlString("binary"), character_options, "The type of data to be constructed." ) );
        
//        argument_rules.push_back( new ArgumentRule( "virtualPopulationSize", Natural::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
//        argument_rules.push_back( new ArgumentRule( "numStates", Natural::getClassTypeSpec(), "The number of states (e.g. 4 for A,C,G and T).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(4) ) );
        
        rules_set = true;
        
    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_readVCF::getClassType(void)
{
    
    static std::string rev_type = "Func_readVCF";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readVCF::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readVCF::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readVCF";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readVCF::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readVCF::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = AbstractHomologousDiscreteCharacterData::getClassTypeSpec();
    return return_typeSpec;
}
