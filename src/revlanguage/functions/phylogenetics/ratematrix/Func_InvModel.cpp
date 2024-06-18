#include "Func_InvModel.h"

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

#include "DistributionChisq.h"
#include "RbMathFunctions.h"
#include "ConcreteTimeReversibleRateMatrix.h"
#include "SiteMixtureModel.h"
#include "GeneratorToSiteModel.h"

using std::vector;
using std::shared_ptr;

namespace Core = RevBayesCore;

Core::SiteMixtureModel* InvModelFunc(const Core::SiteMixtureModel& sub_model, double p)
{
    if (p < 0 or p > 1)
        throw RbException()<<"fnInv: pInv should be in [0,1], but pInv = "<<p;

    int n = sub_model.getNumberOfStates();

    // 1. Compute average root frequencies of the submodel
    auto f = sub_model.componentProbs();
    vector<double> pi(n,0);
    for(int m=0; m<f.size(); m++)
    {
        auto sub_pi = sub_model.getComponent(m).getRootFrequencies();
        for(int l=0;l<pi.size();l++)
            pi[l] += sub_pi[l] * f[m];
    }

    // 2. Construct the INV model component
    vector<double> er(n*(n-1)/2,0);
    auto INV1 = Core::ConcreteTimeReversibleRateMatrix(er,pi,{});
    auto INV2 = std::shared_ptr<const Core::SiteModel>(Core::GeneratorToSiteModel(INV1).clone());
    auto INV = std::shared_ptr<const Core::SiteMixtureModel>(new Core::SiteMixtureModel({INV2},{1.0}));

    auto VAR = shared_ptr<const RevBayesCore::SiteMixtureModel>(sub_model.clone());
    
    auto inv_model = Core::mix_mixture( {VAR, INV}, {1-p, p} )->clone();

    // Here we scaling the model so that the invariant site category does not change the rate.
    // This fits what RevBayes is currently doing.

    // An alternative approach would be to scale the model to the desired rate after it is fully constructed,
    // in phyloCTMC. That which would make it unnecessary to scale the rate here.  It would allow models with
    // more invariant sites to have a smaller rate, which seems natural and is desirable in some contexts.

    if (auto rate = sub_model.rate())
	inv_model->setRate( *rate );

    return inv_model;
}

using namespace RevLanguage;


/** default constructor */
Func_InvModel::Func_InvModel( void ) : TypedFunction<SiteMixtureModel>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_InvModel* Func_InvModel::clone( void ) const
{
    return new Func_InvModel( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::SiteMixtureModel >* Func_InvModel::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< RevBayesCore::SiteMixtureModel >* submodel = static_cast<const SiteMixtureModel &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* pInv = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr< RevBayesCore::SiteMixtureModel >( InvModelFunc, submodel, pInv );
}


/* Get argument rules */
const ArgumentRules& Func_InvModel::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "submodel"       , SiteMixtureModel::getClassTypeSpec(), "Sub-rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "pInv"           , RealPos::getClassTypeSpec(), "The fraction of invariable sites.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_InvModel::getClassType(void)
{

    static std::string rev_type = "Func_InvModel";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_InvModel::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_InvModel::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnInvASRV";

    return f_name;
}


const TypeSpec& Func_InvModel::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
