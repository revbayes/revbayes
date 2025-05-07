#include "Func_ConvertVectorRateMatrix.h"

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
using std::unique_ptr;

namespace Core = RevBayesCore;

Core::RbVector<Core::SiteMixtureModel>* ConvertVectorRateMatrixFunc(const Core::RbVector<Core::RateGenerator>& rate_generators)
{
    auto models = new Core::RbVector<Core::SiteMixtureModel>();
    
    int i=0;
    for(auto& generator: rate_generators)
    {
	i++;
	if (auto matrix = dynamic_cast<const Core::RateMatrix*>(&generator))
	{
	    auto site_model = std::make_shared<const Core::GeneratorToSiteModel>(*matrix);
	    models->push_back( Core::SiteMixtureModel({site_model},{1}) );
	}
	else
	    throw RbException()<<"Converting RateGenerator[] to SiteMixtureModel[]: entry "<<i<<" is not a RateMatrix!";
    }

    return models;
}

using namespace RevLanguage;


/** default constructor */
Func_ConvertVectorRateMatrix::Func_ConvertVectorRateMatrix( void ) : TypedFunction< ModelVector<SiteMixtureModel> >( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_ConvertVectorRateMatrix* Func_ConvertVectorRateMatrix::clone( void ) const
{
    return new Func_ConvertVectorRateMatrix( *this );
}


Core::TypedFunction< Core::RbVector<Core::SiteMixtureModel> >* Func_ConvertVectorRateMatrix::createFunction( void ) const
{
    Core::TypedDagNode< Core::RbVector<Core::RateGenerator> >* matrices = dynamic_cast<const ModelVector<RateGenerator> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::RbVector<Core::SiteMixtureModel> >( ConvertVectorRateMatrixFunc, matrices);
}


/* Get argument rules */
const ArgumentRules& Func_ConvertVectorRateMatrix::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "matrices"           , ModelVector<RateGenerator>::getClassTypeSpec(), "The matrix to convert.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_ConvertVectorRateMatrix::getClassType(void)
{
    static std::string rev_type = "Func_ConvertVectorRateMatrix";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_ConvertVectorRateMatrix::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_ConvertVectorRateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "_RateGenerator[]2SiteMixtureModel[]";

    return f_name;
}


const TypeSpec& Func_ConvertVectorRateMatrix::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
