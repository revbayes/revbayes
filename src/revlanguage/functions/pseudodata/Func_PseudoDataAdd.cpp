#include "Func_PseudoDataAdd.h"

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

Core::PseudoData<double>* PseudoDataAddFunc(Core::PseudoData<double> base, double delta)
{
    Core::PseudoData<double>::func_t add_func = [=](const double& x)
        {
            return base(x - delta);
        };
    return new Core::PseudoData<double>(add_func);
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoDataAdd::Func_PseudoDataAdd( void ) : TypedFunction<PseudoData<Real>>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoDataAdd* Func_PseudoDataAdd::clone( void ) const
{
    return new Func_PseudoDataAdd( *this );
}


Core::TypedFunction< Core::PseudoData<double> >* Func_PseudoDataAdd::createFunction( void ) const
{
    auto base = dynamic_cast<const PseudoData<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    auto delta = dynamic_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoData<double> >( PseudoDataAddFunc, base, delta );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataAdd::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "base",  PseudoData<Real>::getClassTypeSpec(), "The base likelihood function.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "delta", RealPos::getClassTypeSpec(),          "The amount to shift.",          ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataAdd::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataAdd";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataAdd::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataAdd::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "_add";

    return f_name;
}


const TypeSpec& Func_PseudoDataAdd::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
