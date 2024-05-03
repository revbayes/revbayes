#include "Func_ConvertRateMatrix.h"

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

Core::SiteMixtureModel* ConvertRateMatrixFunc(const Core::RateGenerator& model)
{
    if (auto matrix = dynamic_cast<const Core::RateMatrix*>(&model))
    {
	auto site_model = std::make_shared<const Core::GeneratorToSiteModel>(*matrix);
        return new Core::SiteMixtureModel({site_model},{1.0});
    }
    else
        throw RbException()<<"_RateMatrix2MixtureModel: model is not a RateMatrix!";
}

using namespace RevLanguage;


/** default constructor */
Func_ConvertRateMatrix::Func_ConvertRateMatrix( void ) : TypedFunction<SiteMixtureModel>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_ConvertRateMatrix* Func_ConvertRateMatrix::clone( void ) const
{
    return new Func_ConvertRateMatrix( *this );
}


Core::TypedFunction< Core::SiteMixtureModel >* Func_ConvertRateMatrix::createFunction( void ) const
{
    Core::TypedDagNode< Core::RateGenerator >* model = dynamic_cast<const RateMatrix &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::SiteMixtureModel >( ConvertRateMatrixFunc, model);
}


/* Get argument rules */
const ArgumentRules& Func_ConvertRateMatrix::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "matrix"           , RateMatrix::getClassTypeSpec(), "The matrix to convert.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_ConvertRateMatrix::getClassType(void)
{
    static std::string rev_type = "Func_ConvertRateMatrix";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_ConvertRateMatrix::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_ConvertRateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "_RateMatrix2SiteMixtureModel";

    return f_name;
}


const TypeSpec& Func_ConvertRateMatrix::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
