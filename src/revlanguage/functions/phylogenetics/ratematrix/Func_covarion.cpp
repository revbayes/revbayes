#include "Func_covarion.h"
#include "CovarionFunction.h"
#include "Natural.h"
#include "RateMatrix_JC.h"
#include "Real.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_covarion::Func_covarion( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_covarion* Func_covarion::clone( void ) const
{
    
    return new Func_covarion( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_covarion::createFunction( void ) const
{
    
    bool rescale = static_cast<const RlBoolean &>( this->args[3].getVariable()->getRevObject() ).getDagNode()->getValue();
    RevBayesCore::CovarionFunction* f = new RevBayesCore::CovarionFunction( rescale );
    
    // set transition rates for observed character states
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    f->setRateMatrices( rm );
    
    // set transition rates for hidden character states
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* rs = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    f->setRateScalars( rs );
    
    // set transition rates for hidden character states
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* sr = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    f->setSwitchRates( sr );
    
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_covarion::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "RateMatrices", ModelVector<RateGenerator>::getClassTypeSpec(), "Rate matrices for the characters per state.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "RateScalars", ModelVector<RealPos>::getClassTypeSpec(),         "Rate multipliers per state.",   ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "SwitchRates", ModelVector<ModelVector<RealPos> >::getClassTypeSpec(),         "Rates between state.",   ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rescaled",                    RlBoolean::getClassTypeSpec(), "Should the matrix be normalized?",          ArgumentRule::BY_VALUE,              ArgumentRule::ANY, new RlBoolean(true) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_covarion::getClassType(void)
{
    
    static std::string rev_type = "Func_covarion";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_covarion::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_covarion::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnCovarion";
    
    return f_name;
}


const TypeSpec& Func_covarion::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
