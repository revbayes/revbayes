#include "Func_PseudoDataSub2.h"

#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlPseudoData.h"
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

using std::vector;
using std::unique_ptr;

namespace Core = RevBayesCore;

Core::PseudoData<double>* PseudoDataSub2Func(double first, Core::PseudoData<double> second)
{
    Core::PseudoData<double>::func_t sub2_func = [=](const double& x)
        {
            return second(first - x);
        };
    return new Core::PseudoData<double>(sub2_func);
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoDataSub2::Func_PseudoDataSub2( void ) : TypedFunction<PseudoData<Real>>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoDataSub2* Func_PseudoDataSub2::clone( void ) const
{
    return new Func_PseudoDataSub2( *this );
}


Core::TypedFunction< Core::PseudoData<double> >* Func_PseudoDataSub2::createFunction( void ) const
{
    auto first = dynamic_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    auto second = dynamic_cast<const PseudoData<Real> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoData<double> >( PseudoDataSub2Func, first, second );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataSub2::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "first", RealPos::getClassTypeSpec(),          "The amount to subtrace the base from.",          ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "secondPseudoData",  PseudoData<Real>::getClassTypeSpec(), "The base likelihood function.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataSub2::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataSub2";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataSub2::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataSub2::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "_sub";

    return f_name;
}


const TypeSpec& Func_PseudoDataSub2::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
