#include "Func_ScaleVectorSiteMixtureModel.h"

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

Core::RbVector<Core::SiteMixtureModel>* ScaleVectorSiteMixtureModels(const Core::RbVector<Core::SiteMixtureModel>& models,
								    const std::vector<double>& rates)
{
    if (models.size() != rates.size())
	throw RbException()<<"fnScale: got "<<models.size()<<" models, but "<<rates.size()<<" rates";

    auto models2 = new Core::RbVector<Core::SiteMixtureModel>;

    for(int i=0; i<rates.size(); i++)
    {
	models2->push_back(models[i]);
	models2->back().scale(rates[i]);
    }

    return models2;
}

Core::RbVector<Core::SiteMixtureModel>* ScaleVectorSiteMixtureModel(const Core::SiteMixtureModel& model,
								    const std::vector<double>& rates)
{
    auto models2 = new Core::RbVector<Core::SiteMixtureModel>;

    for(int i=0; i<rates.size(); i++)
    {
	models2->push_back(model);
	models2->back().scale(rates[i]);
    }

    return models2;
}

using namespace RevLanguage;


/** default constructor */
Func_ScaleVectorSiteMixtureModel::Func_ScaleVectorSiteMixtureModel( void ) : TypedFunction<ModelVector<SiteMixtureModel>>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_ScaleVectorSiteMixtureModel* Func_ScaleVectorSiteMixtureModel::clone( void ) const
{
    return new Func_ScaleVectorSiteMixtureModel( *this );
}


Core::TypedFunction< Core::RbVector<Core::SiteMixtureModel> >* Func_ScaleVectorSiteMixtureModel::createFunction( void ) const
{
    Core::TypedDagNode< Core::RbVector<double> >* rates = dynamic_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    if ( this->args[0].getVariable()->getRevObject().isType( ModelVector<SiteMixtureModel>::getClassTypeSpec() ) )
    {
	Core::TypedDagNode< Core::RbVector<Core::SiteMixtureModel> >* models = dynamic_cast<const ModelVector<SiteMixtureModel> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
	return Core::generic_function_ptr< Core::RbVector<Core::SiteMixtureModel> >( ScaleVectorSiteMixtureModels, models, rates );
    }
    else if (this->args[0].getVariable()->getRevObject().isType( SiteMixtureModel::getClassTypeSpec() ) )
    {
	Core::TypedDagNode< Core::SiteMixtureModel >* model = dynamic_cast<const SiteMixtureModel &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
	return Core::generic_function_ptr< Core::RbVector<Core::SiteMixtureModel> >( ScaleVectorSiteMixtureModel, model, rates );
    }
    else
	std::abort();
}


/* Get argument rules */
const ArgumentRules& Func_ScaleVectorSiteMixtureModel::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
	std::vector<TypeSpec> modelsTypes;
	modelsTypes.push_back( ModelVector<SiteMixtureModel>::getClassTypeSpec()  );
	modelsTypes.push_back( SiteMixtureModel::getClassTypeSpec()  );
        argumentRules.push_back( new ArgumentRule( "models"    , modelsTypes, "The mixture models to scale.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        argumentRules.push_back( new ArgumentRule( "rates"     , ModelVector<RealPos>::getClassTypeSpec(), "The factors by which to scale the speed of each model.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_ScaleVectorSiteMixtureModel::getClassType(void)
{
    static std::string rev_type = "Func_ScaleVectorSiteMixtureModel";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_ScaleVectorSiteMixtureModel::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_ScaleVectorSiteMixtureModel::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnScale";

    return f_name;
}


const TypeSpec& Func_ScaleVectorSiteMixtureModel::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
