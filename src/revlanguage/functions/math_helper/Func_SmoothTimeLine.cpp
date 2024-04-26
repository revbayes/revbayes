#include <cstddef>
#include <iosfwd>
#include <vector>

#include "SmoothTimeLineFunction.h"
#include "Func_SmoothTimeLine.h"
#include "Natural.h"
#include "RlBoolean.h"
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
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** default constructor */
Func_SmoothTimeLine::Func_SmoothTimeLine( void ) : TypedFunction< ModelVector<RealPos> >( )
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_SmoothTimeLine* Func_SmoothTimeLine::clone( void ) const
{

    return new Func_SmoothTimeLine( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >* Func_SmoothTimeLine::createFunction( void ) const
{

    RevBayesCore::TypedDagNode<double>* max_t                              = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* times = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* values = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::SmoothTimeLineFunction* f = new RevBayesCore::SmoothTimeLineFunction( max_t, times, values );
      return f;

}


/* Get argument rules */
const ArgumentRules& Func_SmoothTimeLine::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( rules_set == false )
    {


        argumentRules.push_back( new ArgumentRule( "maxTime"           , Real::getClassTypeSpec()              , "The maximum time after which we smoothen the timelime out.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "times"             , ModelVector<RealPos>::getClassTypeSpec() , "The times at which the values change.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "values"            , ModelVector<RealPos>::getClassTypeSpec()              , "The values for each time bin.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_SmoothTimeLine::getClassType(void)
{

    static std::string rev_type = "Func_SmoothTimeLine";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_SmoothTimeLine::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_SmoothTimeLine::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnSmoothTimeLine";

    return f_name;
}


const TypeSpec& Func_SmoothTimeLine::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
