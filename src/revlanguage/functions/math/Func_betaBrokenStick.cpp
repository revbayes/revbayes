#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "Func_betaBrokenStick.h"
#include "Integer.h"
#include "NormalizeVectorFunction.h"
#include "ModelVector.h"
#include "Probability.h"
#include "RbUtil.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"

#include "BetaBrokenStickFunction.h"


using namespace RevLanguage;

Func_betaBrokenStick::Func_betaBrokenStick() : TypedFunction<Simplex>() {

}

/* Clone object */
Func_betaBrokenStick* Func_betaBrokenStick::clone( void ) const
{

    return new Func_betaBrokenStick( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Simplex>* Func_betaBrokenStick::createFunction( void ) const
{

    RevBayesCore::TypedDagNode<double>* alpha = static_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* beta = static_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* numCats = static_cast<const Integer &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::BetaBrokenStickFunction *func = new RevBayesCore::BetaBrokenStickFunction( alpha, beta, numCats );

    return func;
}


/** Get argument rules */
const ArgumentRules& Func_betaBrokenStick::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {

        argumentRules.push_back( new ArgumentRule( "alpha"   , RealPos::getClassTypeSpec()  , "The alpha parameter of the beta distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "beta"    , RealPos::getClassTypeSpec()  , "The beta parameter of the beta distribution", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numCats", Natural::getClassTypeSpec()   , "The number of categories.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }

    return argumentRules;
}


/** Get class name of object */
const std::string& Func_betaBrokenStick::getClassName(void)
{

    static std::string rbClassName = "Func_betaBrokenStick";

    return rbClassName;
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_betaBrokenStick::getClassTypeSpec(void)
{

    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rbClass;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_betaBrokenStick::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnBetaBrokenStick";

    return f_name;
}


/** Get type spec */
const TypeSpec& Func_betaBrokenStick::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
