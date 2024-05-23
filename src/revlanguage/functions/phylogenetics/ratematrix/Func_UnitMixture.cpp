#include "Func_UnitMixture.h"

#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlSiteMixtureModel.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "GenericFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RbException.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlBoolean.h"
#include "Integer.h"
#include "Simplex.h"
#include "TypeSpec.h"
#include "ModelVector.h"
#include "RbVector.h"
#include "RbVectorImpl.h"


#include "DistributionChisq.h"
#include "RbMathFunctions.h"
#include "GeneratorToSiteModel.h"

using std::vector;
using std::shared_ptr;

namespace Core = RevBayesCore;


// FIXME -- Make an argument rule so that we can check if the model is a RateMatrix when we create the function.
Core::SiteMixtureModel* UnitMixtureFunc(const Core::RateGenerator& model,
                                         const Core::Simplex& root_freqs,
                                         double rate)
{
    if (root_freqs.size() == model.getNumberOfStates())
    {
	auto site_model = std::shared_ptr<const Core::SiteModel>(new Core::GeneratorToSiteModel(model, root_freqs, rate));
	return new Core::SiteMixtureModel({site_model}, {1.0});
    }
    else
        throw RbException()<<"fnUnitMixture: root frequencies have "<<root_freqs.size()<<" states, but model has "<<model.getNumberOfStates()<<" states.";
}

Core::SiteMixtureModel* UnitMixtureFuncNoFreqs(const Core::RateGenerator& model, double rate)
{
    if (auto matrix = dynamic_cast<const Core::RateMatrix*>(&model))
    {
	auto site_model = std::shared_ptr<const Core::SiteModel>(new Core::GeneratorToSiteModel(*matrix, rate));
	return new Core::SiteMixtureModel({site_model}, {1.0});
    }
    else
        throw RbException()<<"fnUnitMixture: root frequencies not supplied, but model is not a RateMatrix";
}

using namespace RevLanguage;


/** default constructor */
Func_UnitMixture::Func_UnitMixture( void ) : TypedFunction<SiteMixtureModel>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_UnitMixture* Func_UnitMixture::clone( void ) const
{
    return new Func_UnitMixture( *this );
}


Core::TypedFunction< Core::SiteMixtureModel >* Func_UnitMixture::createFunction( void ) const
{
    Core::TypedDagNode< Core::RateGenerator >* model = dynamic_cast<const RateGenerator &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    Core::TypedDagNode< double >* rate = dynamic_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    if (auto& root_freqs_var = *args[1].getVariable(); root_freqs_var.getRevObject() != RevNullObject::getInstance())
    {
        Core::TypedDagNode< Core::Simplex >* root_freqs = dynamic_cast<const Simplex &>( root_freqs_var.getRevObject() ).getDagNode();
        return Core::generic_function_ptr< Core::SiteMixtureModel >( UnitMixtureFunc, model, root_freqs, rate );
    }
    else
        return Core::generic_function_ptr< Core::SiteMixtureModel >( UnitMixtureFuncNoFreqs, model, rate );
}


/* Get argument rules */
const ArgumentRules& Func_UnitMixture::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "model"           , RateGenerator::getClassTypeSpec(), "The mixture models to mix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rootFrequencies" , Simplex::getClassTypeSpec(), "State frequencies at the root.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, nullptr ) );
        argumentRules.push_back( new ArgumentRule( "rate"            , RealPos::getClassTypeSpec(), "Scale by model by this rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0)) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_UnitMixture::getClassType(void)
{
    static std::string rev_type = "Func_UnitMixture";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_UnitMixture::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_UnitMixture::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnUnitMixture";

    return f_name;
}


const TypeSpec& Func_UnitMixture::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
