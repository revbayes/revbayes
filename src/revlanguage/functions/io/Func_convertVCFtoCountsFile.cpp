#include "HomologousDiscreteCharacterData.h"
#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "VCFReader.h"
#include "Ellipsis.h"
#include "Func_convertVCFtoCountsFile.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RevNullObject.h"
#include "RlAbstractDiscreteTaxonData.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlTaxon.h"
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
Func_convertVCFtoCountsFile* Func_convertVCFtoCountsFile::clone( void ) const
{
    
    return new Func_convertVCFtoCountsFile( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_convertVCFtoCountsFile::execute( void )
{
    size_t arg_index = 0;
    
    // get the information from the arguments for reading the file
    const std::string& fn_vcf       = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    const std::string& fn_counts    = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();

    const std::string type = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();

    const std::string& ploidy_str       = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    RevBayesCore::VCFReader::PLOIDY ploidy = RevBayesCore::VCFReader::DIPLOID;
    if ( ploidy_str == "haploid" )
    {
        ploidy = RevBayesCore::VCFReader::HAPLOID;
    }
    
    const std::string& unkown_str       = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    RevBayesCore::VCFReader::UNKOWN_TREATMENT unkown = RevBayesCore::VCFReader::MISSING;
    if ( unkown_str == "alternative" )
    {
        unkown = RevBayesCore::VCFReader::ALTERNATIVE;
    }
    else if ( unkown_str == "reference" )
    {
        unkown = RevBayesCore::VCFReader::REFERENCE;
    }
    
    RevBayesCore::VCFReader vcf_reader = RevBayesCore::VCFReader( fn_vcf, ploidy, unkown, false );

    const RevBayesCore::RbVector<RevBayesCore::Taxon>& taxa  = static_cast< const ModelVector<Taxon> &>( args[arg_index++].getVariable()->getRevObject() ).getValue();

    const std::string& chrom = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::AbstractDiscreteTaxonData* ref_genome = NULL;
    RevObject& tmp_ref = args[arg_index++].getVariable()->getRevObject();
    if ( tmp_ref != RevNullObject::getInstance() )
    {
        ref_genome = static_cast<const AbstractDiscreteTaxonData&>( tmp_ref ).getValue().clone();
    }

    
    long thinning    = static_cast< const Natural&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    long skip_first  = static_cast< const Natural&>( args[arg_index++].getVariable()->getRevObject() ).getValue();

    vcf_reader.convertToCountsFile( fn_counts, taxa, type, chrom, ref_genome, thinning, skip_first );
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_convertVCFtoCountsFile::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        argument_rules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "Relative or absolute name of the vcf file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "countsfile", RlString::getClassTypeSpec(), "Relative or absolute name of the counts file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        std::vector<std::string> character_options;
        character_options.push_back( "DNA" );
        character_options.push_back( "binary" );
        argument_rules.push_back( new OptionRule( "type", new RlString("binary"), character_options, "The type of data to be constructed." ) );
        
        std::vector<std::string> ploidy_options;
        ploidy_options.push_back( "diploid" );
        ploidy_options.push_back( "haploid" );
        argument_rules.push_back( new OptionRule( "ploidy", new RlString("diploid"), ploidy_options, "The ploidy type." ) );
        
        std::vector<std::string> unknown_options;
        unknown_options.push_back( "missing" );
        unknown_options.push_back( "reference" );
        unknown_options.push_back( "alternative" );
        argument_rules.push_back( new OptionRule( "unkownTreatment", new RlString("missing"), unknown_options, "How to treat the '.' character, i.e., if the state was unkown." ) );
        
        argument_rules.push_back( new ArgumentRule( "taxa"     , ModelVector<Taxon>::getClassTypeSpec()                     , "The taxa to match the individuals to species/populations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        argument_rules.push_back( new ArgumentRule( "chrom",       RlString::getClassTypeSpec(), "Name of the chromosome we want to extract. If empty, then all chromosomes are used.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "" ) ) );

        argument_rules.push_back( new ArgumentRule( "reference",   AbstractDiscreteTaxonData::getClassTypeSpec(), "The reference genome if we want to fill in the characters that are not coded in the VCF file, i.e., if the VCF file only contains SNPs but we need the full sequences.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );

        
        argument_rules.push_back( new ArgumentRule( "thinning", Natural::getClassTypeSpec(), "If thinning is larger than 1, then we only take the i-th entry of the VCF.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );
        argument_rules.push_back( new ArgumentRule( "skipFirst", Natural::getClassTypeSpec(), "Skip the first n entries.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );

        rules_set = true;
        
    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_convertVCFtoCountsFile::getClassType(void)
{
    
    static std::string rev_type = "Func_convertVCFtoCountsFile";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_convertVCFtoCountsFile::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_convertVCFtoCountsFile::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "convertVCFtoCountsFile";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_convertVCFtoCountsFile::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_convertVCFtoCountsFile::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = AbstractHomologousDiscreteCharacterData::getClassTypeSpec();
    return return_typeSpec;
}
