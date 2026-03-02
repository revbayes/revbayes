#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_computeWithinSpeciesVarianceFromCharacterData.h"
#include "ModelVector.h"
#include "ComputeWithinSpeciesVarianceFromCharacterDataFunction.h"
#include "RlContinuousCharacterData.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "Real.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "RlString.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class ContinuousCharacterData; }
namespace RevBayesCore { class Taxon; }
namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/** Default constructor */
Func_computeWithinSpeciesVarianceFromCharacterData::Func_computeWithinSpeciesVarianceFromCharacterData( void ) : TypedFunction<ModelVector<Real> >()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_computeWithinSpeciesVarianceFromCharacterData* Func_computeWithinSpeciesVarianceFromCharacterData::clone( void ) const
{

    return new Func_computeWithinSpeciesVarianceFromCharacterData( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >* Func_computeWithinSpeciesVarianceFromCharacterData::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode<RevBayesCore::ContinuousCharacterData>* data = static_cast<const ContinuousCharacterData &>( args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode<std::int64_t>* v_site = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode<std::int64_t>* n_site = static_cast<const Natural &>( args[2].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::RbVector<RevBayesCore::Taxon>& taxa = static_cast<const ModelVector<Taxon> &>( args[3].getVariable()->getRevObject() ).getValue();

    const std::string& tr = static_cast<const RlString &>( args[4].getVariable()->getRevObject() ).getValue();

    RevBayesCore::ComputeWithinSpeciesVarianceFromCharacterDataFunction::MISSING_TREATMENT mtr;
    if (tr == "mean")
    {
        mtr = RevBayesCore::ComputeWithinSpeciesVarianceFromCharacterDataFunction::MISSING_TREATMENT::MEAN;
    }
    else if (tr == "median")
    {
        mtr = RevBayesCore::ComputeWithinSpeciesVarianceFromCharacterDataFunction::MISSING_TREATMENT::MEDIAN;
    }
    else if (tr == "none")
    {
        mtr = RevBayesCore::ComputeWithinSpeciesVarianceFromCharacterDataFunction::MISSING_TREATMENT::NONE;
    }
    else
    {
        throw RbException("argument missingVarianceTreatment must be one of \"mean\", \"median\" or \"none\"");
    }


    RevBayesCore::ComputeWithinSpeciesVarianceFromCharacterDataFunction* f = new RevBayesCore::ComputeWithinSpeciesVarianceFromCharacterDataFunction( data, v_site, n_site, taxa, mtr );

    return f;
}


/** Get argument rules */
const ArgumentRules& Func_computeWithinSpeciesVarianceFromCharacterData::getArgumentRules( void ) const
{

    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;

    if ( rules_set == false )
    {

        argument_rules.push_back( new ArgumentRule( "data", ContinuousCharacterData::getClassTypeSpec(), "The character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "variance_site", Natural::getClassTypeSpec(), "The site that stores the trait variance.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "num_sample_site", Natural::getClassTypeSpec(), "The site that stores the number of samples per species.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "taxa", ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa which have species and individual names.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        std::vector<std::string> missingTreatmentTypes;
        missingTreatmentTypes.push_back( "mean" );
        missingTreatmentTypes.push_back( "median" );
        missingTreatmentTypes.push_back( "none" );
        argument_rules.push_back( new OptionRule ("missingVarianceTreatment", new RlString("none"), missingTreatmentTypes, "The within-species variance for species with only one sample. Options \"mean\" and \"median\" return the mean/median of within-species variance of all species with multiple samples. Option \"none\" returns -1 and requires the user to manually specify the within-species variance afterwards.") );

        rules_set = true;
    }

    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_computeWithinSpeciesVarianceFromCharacterData::getClassType(void)
{

    static std::string rev_type = "Func_computeWithinSpeciesVarianceFromCharacterData";

    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_computeWithinSpeciesVarianceFromCharacterData::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_computeWithinSpeciesVarianceFromCharacterData::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnComputeWithinSpeciesVarianceFromCharacterData";

    return f_name;
}

std::vector<std::string> Func_computeWithinSpeciesVarianceFromCharacterData::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "fnWithinSpeciesVarianceFromData" );

    return a_names;
}

/** Get type spec */
const TypeSpec& Func_computeWithinSpeciesVarianceFromCharacterData::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
