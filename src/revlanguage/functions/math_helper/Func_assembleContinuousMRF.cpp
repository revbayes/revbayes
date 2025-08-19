#include <cstddef>
#include <iosfwd>
#include <vector>

#include "AssembleOrder1ContinuousMRFFunction.h"
#include "AssembleOrder2ContinuousMRFFunction.h"
#include "Func_assembleContinuousMRF.h"
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
Func_assembleContinuousMRF::Func_assembleContinuousMRF( void ) : TypedFunction< ModelVector<RealPos> >( )
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_assembleContinuousMRF* Func_assembleContinuousMRF::clone( void ) const
{

    return new Func_assembleContinuousMRF( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >* Func_assembleContinuousMRF::createFunction( void ) const
{

    RevBayesCore::TypedDagNode<double>* theta1                              = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* increments = static_cast<const ModelVector<Real> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    bool theta1_is_log_scale                                                = static_cast<const RlBoolean &>( args[4].getVariable()->getRevObject() ).getValue();
    size_t order                                                            = static_cast<const Natural &>( args[5].getVariable()->getRevObject() ).getValue();

    
    RevBayesCore::TypedDagNode<double>* beta                                = static_cast<const Real &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* predictors = NULL;
    if (this->args[3].getVariable()->getRevObject() != RevNullObject::getInstance())
    {
        predictors = static_cast<const ModelVector<Real> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    }
    
    if (order == 1)
    {
      RevBayesCore::AssembleOrder1ContinuousMRFFunction* f = new RevBayesCore::AssembleOrder1ContinuousMRFFunction( theta1, increments, predictors, beta, theta1_is_log_scale );
      return f;

    }
    else if (order == 2)
    {
      RevBayesCore::AssembleOrder2ContinuousMRFFunction* f = new RevBayesCore::AssembleOrder2ContinuousMRFFunction( theta1, increments, theta1_is_log_scale );
      return f;

    }
    else
    {
      throw(RbException("Unsupported order selected in assembleContinuousMRF"));
    }
}


/* Get argument rules */
const ArgumentRules& Func_assembleContinuousMRF::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( rules_set == false )
    {


        argumentRules.push_back( new ArgumentRule( "initialValue"           , Real::getClassTypeSpec()              , "The first value in the MRF.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "increments"             , ModelVector<Real>::getClassTypeSpec() , "The increments of the process, assumed to be on the log-scale.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "beta"                   , Real::getClassTypeSpec()              , "The correlation coefficient for the predictor variables.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "predictors"             , ModelVector<Real>::getClassTypeSpec() , "The predictor variables, assumed to be.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "initialValueIsLogScale" , RlBoolean::getClassTypeSpec()         , "Is valueInitial on the log-scale, like the increments?", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "order"                  , Natural::getClassTypeSpec()           , "The order of the MRF (1 or 2).", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_assembleContinuousMRF::getClassType(void)
{

    static std::string rev_type = "Func_assembleContinuousMRF";

	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_assembleContinuousMRF::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_assembleContinuousMRF::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnassembleContinuousMRF";

    return f_name;
}


const TypeSpec& Func_assembleContinuousMRF::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
