#include "PosteriorPredictiveProbabilityFunction.h"
#include "Func_posteriorPredictiveProbability.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_posteriorPredictiveProbability::Func_posteriorPredictiveProbability( void ) : TypedFunction<Simplex>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_posteriorPredictiveProbability* Func_posteriorPredictiveProbability::clone( void ) const
{
    
    return new Func_posteriorPredictiveProbability( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::Simplex >* Func_posteriorPredictiveProbability::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* v = static_cast<const ModelVector<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* x = static_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::PosteriorPredictiveProbabilityFunction* f = new RevBayesCore::PosteriorPredictiveProbabilityFunction( v, x );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_posteriorPredictiveProbability::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "v", ModelVector<Real>::getClassTypeSpec(), "The simulated statistics.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec()             , "The opbserved value."     , ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_posteriorPredictiveProbability::getClassType(void)
{
    
    static std::string rev_type = "Func_posteriorPredictiveProbability";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_posteriorPredictiveProbability::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_posteriorPredictiveProbability::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "posteriorPredictiveProbability";
    
    return f_name;
}


const TypeSpec& Func_posteriorPredictiveProbability::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
