#include "Func_MixtureModel.h"

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
#include "SiteMixtureModel.h"

using std::vector;
using std::unique_ptr;

namespace Core = RevBayesCore;

Core::SiteMixtureModel* MixtureModelFunc(const Core::RbVector<Core::SiteMixtureModel>& models,
					 const Core::Simplex& fractions)
{
    if (models.size() == 0)
        throw RbException()<<"Cannot have a mixture of 0 models!";

    if (fractions.size() != models.size())
        throw RbException()<<"Got "<<models.size()<<" models but "<<fractions.size()<<" probabilities.";

    // convert vector<T> to vector<shared_ptr<const T>>
    std::vector<std::shared_ptr<const Core::SiteMixtureModel>> models2;
    for(auto& model: models)
	models2.push_back(std::shared_ptr<const Core::SiteMixtureModel>(model.clone()));

    // It would be nice to be able to pass a smart pointer here.
    return Core::mix_mixture(models2, fractions)->clone();
}

Core::SiteMixtureModel* MixtureModelFuncEqualWeights(const Core::RbVector<Core::SiteMixtureModel>& models)
{
    int n = models.size();
    std::vector<double> fractions(n, 1.0/n);

    return MixtureModelFunc(models, fractions);
}

using namespace RevLanguage;


/** default constructor */
Func_MixtureModel::Func_MixtureModel( void ) : TypedFunction<SiteMixtureModel>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_MixtureModel* Func_MixtureModel::clone( void ) const
{
    return new Func_MixtureModel( *this );
}


Core::TypedFunction< Core::SiteMixtureModel >* Func_MixtureModel::createFunction( void ) const
{
    Core::TypedDagNode< Core::RbVector<Core::SiteMixtureModel> >* models = static_cast<const ModelVector<SiteMixtureModel> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    if (auto& fractions_var = *args[1].getVariable(); fractions_var.getRevObject() != RevNullObject::getInstance())
    {
	Core::TypedDagNode< Core::Simplex >* fractions = static_cast<const Simplex &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
        return Core::generic_function_ptr< Core::SiteMixtureModel >( MixtureModelFunc, models, fractions );
    }
    else
        return Core::generic_function_ptr< Core::SiteMixtureModel >( MixtureModelFuncEqualWeights, models );
}


/* Get argument rules */
const ArgumentRules& Func_MixtureModel::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "models"    , ModelVector<SiteMixtureModel>::getClassTypeSpec(), "The mixture models to mix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "fractions" , Simplex::getClassTypeSpec(), "The probability of each model.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, nullptr) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_MixtureModel::getClassType(void)
{
    static std::string rev_type = "Func_MixtureModel";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_MixtureModel::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_MixtureModel::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnMixtureASRV";

    return f_name;
}


const TypeSpec& Func_MixtureModel::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
