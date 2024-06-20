#include "Func_ScaleSiteMixtureModel.h"

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

Core::SiteMixtureModel* ScaleSiteMixtureModel(const Core::SiteMixtureModel& model,
					      double factor)
{
    // convert vector<T> to vector<shared_ptr<const T>>
    auto model2 = model.clone();
    model2->scale(factor);
    return model2;
}

using namespace RevLanguage;


/** default constructor */
Func_ScaleSiteMixtureModel::Func_ScaleSiteMixtureModel( void ) : TypedFunction<SiteMixtureModel>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_ScaleSiteMixtureModel* Func_ScaleSiteMixtureModel::clone( void ) const
{
    return new Func_ScaleSiteMixtureModel( *this );
}


Core::TypedFunction< Core::SiteMixtureModel >* Func_ScaleSiteMixtureModel::createFunction( void ) const
{
    Core::TypedDagNode< Core::SiteMixtureModel >* model = static_cast<const SiteMixtureModel &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    Core::TypedDagNode< double >* factor = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();


    return Core::generic_function_ptr< Core::SiteMixtureModel >( ScaleSiteMixtureModel, model, factor );
}


/* Get argument rules */
const ArgumentRules& Func_ScaleSiteMixtureModel::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "model"    , SiteMixtureModel::getClassTypeSpec(), "The mixture model to scale.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rate"     , RealPos::getClassTypeSpec(), "The factor by which to scale the speed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_ScaleSiteMixtureModel::getClassType(void)
{
    static std::string rev_type = "Func_ScaleSiteMixtureModel";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_ScaleSiteMixtureModel::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_ScaleSiteMixtureModel::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnScale";

    return f_name;
}


const TypeSpec& Func_ScaleSiteMixtureModel::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
