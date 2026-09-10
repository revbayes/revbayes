#include <iosfwd>
#include <string>
#include <vector>

#include "PoMoRateMatrixFunction.h"
#include "Func_pomo.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Real.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "RateGenerator.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlRateGenerator.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"
#include "RealPos.h" // IWYU pragma: keep

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

/** default constructor */
Func_pomo::Func_pomo( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_pomo* Func_pomo::clone( void ) const
{
    
    return new Func_pomo( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_pomo::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator >* q = static_cast<const RateMatrix &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    //RevBayesCore::TypedDagNode< double >* root_pol = static_cast<const double &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* fit = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< std::int64_t >* n = static_cast<const Natural &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::PoMoRateMatrixFunction* f = new RevBayesCore::PoMoRateMatrixFunction( n, q, fit );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_pomo::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "mutationRates", RateGenerator::getClassTypeSpec()    , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "fitness"      , ModelVector<Real>::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "virtualNe"    , Natural::getClassTypeSpec()          , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_pomo::getClassType(void)
{
    
    static std::string rev_type = "Func_pomo";
    
	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_pomo::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_pomo::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPoMo";
    
    return f_name;
}


const TypeSpec& Func_pomo::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
