#include "AssembleOrder1ContinuousMRFFunction.h"
#include "AssembleOrder2ContinuousMRFFunction.h"
#include "Func_assembleContinuousMRF.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

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
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* increments = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    bool theta1_is_log_scale                                                = static_cast<const RlBoolean &>( args[2].getVariable()->getRevObject() ).getValue();
    size_t order                                                            = static_cast<const Natural &>( args[3].getVariable()->getRevObject() ).getValue();

    if (order == 1)
    {
      RevBayesCore::AssembleOrder1ContinuousMRFFunction* f = new RevBayesCore::AssembleOrder1ContinuousMRFFunction( theta1, increments, theta1_is_log_scale );
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

    if ( !rules_set )
    {


        argumentRules.push_back( new ArgumentRule( "initialValue"           , Real::getClassTypeSpec()              , "The first value in the MRF.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "increments"             , ModelVector<Real>::getClassTypeSpec() , "The increments of the process, assumed to be on the log-scale.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
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
