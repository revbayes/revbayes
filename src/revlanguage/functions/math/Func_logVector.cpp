#include <iosfwd>
#include <string>
#include <vector>

#include "LogVectorFunction.h"
#include "Func_logVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** default constructor */
Func_logVector::Func_logVector( void ) : TypedFunction< ModelVector<Real> >( )
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_logVector* Func_logVector::clone( void ) const
{

    return new Func_logVector( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >* Func_logVector::createFunction( void ) const
{

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* a = static_cast<const ModelVector<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* b = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::LogVectorFunction* f = new RevBayesCore::LogVectorFunction( a, b );

    return f;
}


/* Get argument rules */
const ArgumentRules& Func_logVector::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {

        argumentRules.push_back( new ArgumentRule( "x", ModelVector<Real>::getClassTypeSpec(), "A vector of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "base", RealPos::getClassTypeSpec(), "The base of the logarithm.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_logVector::getClassType(void)
{

    static std::string rev_type = "Func_logVector";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_logVector::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<ModelVector<Real> >::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_logVector::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "log";

    return f_name;
}


const TypeSpec& Func_logVector::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
