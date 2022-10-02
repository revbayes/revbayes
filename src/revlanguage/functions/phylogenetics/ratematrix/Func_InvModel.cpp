#include "Func_GammaRateModel.h"

#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlMixtureModel.h"
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
#include "ConcreteMixtureModel.h"

using std::vector;
using std::unique_ptr;

// wait, so, how do we encode the operations like scale(model), rate(model), etc?
// well, ratemodel should perhaps be its OWN interface!

// how about things like the m3 model?
// also, could we create this using a FUNCTION?

RevBayesCore::ConcreteMixtureModel* InvModelFunc(const RevBayesCore::MixtureModel& sub_model, double p)
{
    using namespace RevBayesCore;

    if (p < 0 or p > 1)
        throw RbException()<<"fnInv: pInv should be in [0,1], but pInv = "<<p;

    vector<unique_ptr<RevBayesCore::MixtureModel>> sub_models;
    vector<double> fractions = {1-p, p};

    // The non-invariant part.
    sub_models.push_back( unique_ptr<RevBayesCore::MixtureModel>(sub_model.clone()) );

    // sub_models.push_back( invariant );
    
    return new RevBayesCore::ConcreteMixtureModel(sub_models, fractions);
}

using namespace RevLanguage;


/** default constructor */
Func_GammaRateModel::Func_GammaRateModel( void ) : TypedFunction<MixtureModel>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_GammaRateModel* Func_GammaRateModel::clone( void ) const
{
    return new Func_GammaRateModel( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::MixtureModel >* Func_GammaRateModel::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< RevBayesCore::MixtureModel >* submodel = static_cast<const MixtureModel &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* pInv = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr< RevBayesCore::MixtureModel >( InvModelFunc, submodel, pInv );
}


/* Get argument rules */
const ArgumentRules& Func_GammaRateModel::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "submodel"       , MixtureModel::getClassTypeSpec(), "Sub-rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "pInv"           , RealPos::getClassTypeSpec(), "The fraction of invariable sites.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_GammaRateModel::getClassType(void)
{

    static std::string rev_type = "Func_GammaRateModel";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_GammaRateModel::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_GammaRateModel::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnGamma";

    return f_name;
}


const TypeSpec& Func_GammaRateModel::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
