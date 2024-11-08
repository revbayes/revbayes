#include "Func_GammaRateModel.h"

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
#include "SiteMixtureModel.h"

using std::vector;
using std::unique_ptr;

namespace Core = RevBayesCore;

vector<double> gamma_rates(double a, int nCats, bool median)
{
    using namespace RevBayesCore;

    double b = a;

    double factor = a / b * nCats;
    std::vector<double> rates(nCats, 0);

    if (median) {
        /* the median value for each category is used to represent all of the values
        in that category */
        double interval = 1.0 / (2.0 * nCats);
        for (int i=0; i<nCats; i++)
            rates[i] = RbStatistics::ChiSquare::quantile((i * 2.0 + 1.0) * interval, 2.0 * a) / (2.0 * b);
        double t = 0.0;
        for (int i=0; i<nCats; i++)
            t += rates[i];
        for (int i=0; i<nCats; i++)
            rates[i] *= factor / t;
    }
    else
    {
        /* the mean value for each category is used to represent all of the values
        in that category */
        /* calculate the points in the gamma distribution */
        for (int i=0; i<nCats-1; i++)
            rates[i] = RbStatistics::ChiSquare::quantile((i + 1.0) / nCats, 2.0 * a) / (2.0 * b);
        /* calculate the cumulative values */
        for (int i=0; i<nCats-1; i++)
            rates[i] = RbMath::incompleteGamma(rates[i] * b, a + 1.0);
        rates[nCats-1] = 1.0;
        /* calculate the relative values and rescale */
        for (int i=nCats-1; i>0; i--){
            rates[i] -= rates[i-1];
            rates[i] *= factor;
        }
        rates[0] *= factor;
    }

    return rates;
}

Core::SiteMixtureModel* GammaRateModelFunc(const Core::SiteMixtureModel& submodel, double a, int nCats, Core::Boolean median)
{
    using namespace RevBayesCore;

    if (nCats < 1)
        throw RbException()<<"GammaRateModel: number of components must be at least 1, but got "<<nCats;

    if (a <= 0)
        throw RbException()<<"GammaRateModel: alpha should be positive, but alpha = "<<a;

    vector<double> fraction(nCats, 1.0/nCats);

    vector<double> rates = gamma_rates(a, nCats, median);

    // It would be nice to be able to pass a smart pointer here.
    return scaled_mixture( submodel, fraction, rates )->clone();
}

using namespace RevLanguage;


/** default constructor */
Func_GammaRateModel::Func_GammaRateModel( void ) : TypedFunction<SiteMixtureModel>( )
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


Core::TypedFunction< Core::SiteMixtureModel >* Func_GammaRateModel::createFunction( void ) const
{
    Core::TypedDagNode< Core::SiteMixtureModel >* submodel = dynamic_cast<const SiteMixtureModel &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* a = dynamic_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< long >* nCats = dynamic_cast<const Integer &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< Core::Boolean >* median = dynamic_cast<const RlBoolean &>( this->args[3].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::SiteMixtureModel >( GammaRateModelFunc, submodel, a, nCats, median);
}


/* Get argument rules */
const ArgumentRules& Func_GammaRateModel::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "submodel"       , SiteMixtureModel::getClassTypeSpec(), "Sub-model.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "alpha"          , RealPos::getClassTypeSpec(), "The alpha parameter of the gamma distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "n"              , Integer::getClassTypeSpec(), "The number of bins to approximate the gamma distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Integer(4) ) );
        argumentRules.push_back( new ArgumentRule( "median"         , RlBoolean::getClassTypeSpec(), "Should we use the median (or the mean)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(false) ) );

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
    std::string f_name = "fnGammaASRV";

    return f_name;
}


const TypeSpec& Func_GammaRateModel::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
