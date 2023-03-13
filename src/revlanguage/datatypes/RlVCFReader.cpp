#include "VCFReader.h"

#include <ostream>
#include <string>
#include <vector>

#include "RlAbstractDiscreteTaxonData.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Delimiter.h"
#include "OptionRule.h"
#include "RevObject.h"
#include "Real.h"
#include "RlVCFReader.h"
#include "RlDemographicFunction.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlMatrixRealPos.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevLanguage { class Argument; }


using namespace RevLanguage;

VCFReader::VCFReader() : WorkspaceToCoreWrapperObject<RevBayesCore::VCFReader>()
{
    
    // computing stats
    
    ArgumentRules* stats_arg_rules = new ArgumentRules();
    
    stats_arg_rules->push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "Relative or absolute base for the name of the statistics files.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    stats_arg_rules->push_back( new ArgumentRule( "taxa"     , ModelVector<Taxon>::getClassTypeSpec()                     , "The taxa to match the individuals to species/populations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

    methods.addFunction(new MemberProcedure( "computeStatistics", RevNullObject::getClassTypeSpec(), stats_arg_rules) );
    
    
    
    
    ArgumentRules* convert_to_nexus_arg_rules = new ArgumentRules();
    
    convert_to_nexus_arg_rules->push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "Relative or absolute base for the name of the statistics files.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    convert_to_nexus_arg_rules->push_back( new ArgumentRule( "taxa"     , ModelVector<Taxon>::getClassTypeSpec()                     , "The taxa to match the individuals to species/populations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    convert_to_nexus_arg_rules->push_back( new ArgumentRule( "chrom",       RlString::getClassTypeSpec(), "Name of the chromosome we want to extract. If empty, then all chromosomes are used.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "" ) ) );
    convert_to_nexus_arg_rules->push_back( new ArgumentRule( "reference",   AbstractDiscreteTaxonData::getClassTypeSpec(), "The reference genome if we want to fill in the characters that are not coded in the VCF file, i.e., if the VCF file only contains SNPs but we need the full sequences.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );

    convert_to_nexus_arg_rules->push_back( new ArgumentRule( "thinning", Natural::getClassTypeSpec(), "If thinning is larger than 1, then we only take the i-th entry of the VCF.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );
    convert_to_nexus_arg_rules->push_back( new ArgumentRule( "skipFirst", Natural::getClassTypeSpec(), "Skip the first n entries.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );

    methods.addFunction(new MemberProcedure( "convertToNexusFile", RevNullObject::getClassTypeSpec(), convert_to_nexus_arg_rules) );

    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
VCFReader* VCFReader::clone(void) const
{
    
    return new VCFReader(*this);
}


void VCFReader::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // get the parameter values
    const std::string& fn_vcf       = static_cast<const RlString&>( filename->getRevObject() ).getValue();
    
    
    const std::string& ploidy_str       = static_cast<const RlString&>( ploidy->getRevObject() ).getValue();
    RevBayesCore::VCFReader::PLOIDY ploidy_type = RevBayesCore::VCFReader::DIPLOID;
    if ( ploidy_str == "haploid" )
    {
        ploidy_type = RevBayesCore::VCFReader::HAPLOID;
    }
    
    const std::string& unkown_str       = static_cast<const RlString&>( unknown->getRevObject() ).getValue();
    RevBayesCore::VCFReader::UNKOWN_TREATMENT unkown = RevBayesCore::VCFReader::MISSING;
    if ( unkown_str == "alternative" )
    {
        unkown = RevBayesCore::VCFReader::ALTERNATIVE;
    }
    else if ( unkown_str == "reference" )
    {
        unkown = RevBayesCore::VCFReader::REFERENCE;
    }
    
    value = new RevBayesCore::VCFReader( fn_vcf, ploidy_type, unkown, false );
    
}


/* Map calls to member methods */
RevPtr<RevVariable> VCFReader::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "computeStatistics")
    {
        found = true;
        size_t arg_index = 0;
        
        const std::string& fn    = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
        const RevBayesCore::RbVector<RevBayesCore::Taxon>& taxa  = static_cast< const ModelVector<Taxon> &>( args[arg_index++].getVariable()->getRevObject() ).getValue();

        value->computeMonomorphicVariableStatistics(fn, taxa);
        
        return NULL;
    }
    else if (name == "convertToNexusFile")
    {
        found = true;
        size_t arg_index = 0;
        
        const std::string& fn    = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
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

        value->convertToNexusFile( fn, "DNA", chrom, ref_genome, taxa, thinning, skip_first );
        
        return NULL;
    }
    
    return WorkspaceToCoreWrapperObject<RevBayesCore::VCFReader>::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& VCFReader::getClassType(void)
{
    
    static std::string rev_type = "VCFReader";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& VCFReader::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::VCFReader>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string VCFReader::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "VCFReader";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& VCFReader::getParameterRules(void) const
{
    
    static MemberRules argument_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "Relative or absolute name of the vcf file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        std::vector<std::string> ploidy_options;
        ploidy_options.push_back( "diploid" );
        ploidy_options.push_back( "haploid" );
        argument_rules.push_back( new OptionRule( "ploidy", new RlString("diploid"), ploidy_options, "The ploidy type." ) );
        
        std::vector<std::string> unknown_options;
        unknown_options.push_back( "missing" );
        unknown_options.push_back( "reference" );
        unknown_options.push_back( "alternative" );
        argument_rules.push_back( new OptionRule( "unkownTreatment", new RlString("missing"), unknown_options, "How to treat the '.' character, i.e., if the state was unkown." ) );
        
        rules_set = true;
    }
    
    return argument_rules;
}


/** Get type spec */
const TypeSpec& VCFReader::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void VCFReader::printValue(std::ostream &o) const
{
    
    o << "VCFReader";
}


/** Set a member variable */
void VCFReader::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "file")
    {
        filename = var;
    }
    else if ( name == "unkownTreatment")
    {
        unknown = var;
    }
    else if ( name == "ploidy")
    {
        ploidy = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}
