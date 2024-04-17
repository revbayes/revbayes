#include "Func_MixtureModel.h"

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
#include "ModelVector.h"
#include "RbVector.h"
#include "RbVectorImpl.h"


#include "DistributionChisq.h"
#include "RbMathFunctions.h"
#include "ConcreteMixtureModel.h"
#include "UnitMixtureModel.h"

using std::vector;
using std::unique_ptr;

namespace Core = RevBayesCore;

Core::ConcreteMixtureModel* MixtureModelFunc(const Core::RbVector<Core::SiteMixtureModel>& models,
                                             const Core::Simplex& fractions,
                                             const Core::RbVector<double>& scales)
{
    if (models.size() == 0)
        throw RbException()<<"Cannot have a mixture of 0 models!";

    if (fractions.size() != models.size())
        throw RbException()<<"Got "<<models.size()<<" models but "<<fractions.size()<<" probabilities.";

    if (scales.size() != models.size())
        throw RbException()<<"Got "<<models.size()<<" models but "<<scales.size()<<" scaling factors.";

    vector<unique_ptr<Core::SiteMixtureModel>> model_ptrs;
    for(int i=0;i<models.size();i++)
    {
        auto model_ptr = std::unique_ptr<Core::SiteMixtureModel>(models[i].clone());
        model_ptr->scale(scales[i]);
        model_ptrs.push_back( std::move(model_ptr) );
    }
    
    return new Core::ConcreteMixtureModel(model_ptrs, fractions);
}

RevBayesCore::ConcreteMixtureModel* MixtureModelFunc2(const RevBayesCore::RbVector<RevBayesCore::SiteMixtureModel>& models,
                                                      const RevBayesCore::Simplex& fractions)
{
    vector<double> scales(models.size(), 1.0);
    return MixtureModelFunc(models, fractions, scales);
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

    Core::TypedDagNode< Core::Simplex >* fractions = static_cast<const Simplex &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    if (auto& scales_var = *args[2].getVariable(); scales_var.getRevObject() != RevNullObject::getInstance())
    {
        Core::TypedDagNode< Core::RbVector<double> >* scales = dynamic_cast<const ModelVector<RealPos> &>( scales_var.getRevObject() ).getDagNode();
        return Core::generic_function_ptr< Core::SiteMixtureModel >( MixtureModelFunc, models, fractions, scales );
    }
    else
        return Core::generic_function_ptr< Core::SiteMixtureModel >( MixtureModelFunc2, models, fractions );
}


/* Get argument rules */
const ArgumentRules& Func_MixtureModel::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "models"    , ModelVector<SiteMixtureModel>::getClassTypeSpec(), "The mixture models to mix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "fractions" , Simplex::getClassTypeSpec(), "The probability of each model.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rates"     , ModelVector<RealPos>::getClassTypeSpec(), "The rate to scale each model.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, nullptr) );

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
    std::string f_name = "fnMixture";

    return f_name;
}


const TypeSpec& Func_MixtureModel::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
