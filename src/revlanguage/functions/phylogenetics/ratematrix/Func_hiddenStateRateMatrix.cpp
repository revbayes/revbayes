#include "Func_hiddenStateRateMatrix.h"

#include "HiddenStateRateMatrixFunction.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlRateGenerator.h"
#include "TypeSpec.h"

namespace RevLanguage { class Real; }

using namespace RevLanguage;

/** default constructor */
Func_hiddenStateRateMatrix::Func_hiddenStateRateMatrix( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_hiddenStateRateMatrix* Func_hiddenStateRateMatrix::clone( void ) const
{
    
    return new Func_hiddenStateRateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_hiddenStateRateMatrix::createFunction( void ) const
{
    
    bool rescale = static_cast<const RlBoolean &>( this->args[2].getVariable()->getRevObject() ).getDagNode()->getValue();
    RevBayesCore::HiddenStateRateMatrixFunction* f = new RevBayesCore::HiddenStateRateMatrixFunction( rescale );
    
    // set transition rates for observed character states
    if ( this->args[0].getVariable()->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* tr = static_cast<const ModelVector<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
        f->setObservedTransitionRates( tr );
    }
    else if ( this->args[0].getVariable()->getRevObject().isType( RateGenerator::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RateGenerator >* tr = static_cast<const RateGenerator &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
        f->setObservedTransitionRates( tr );
    }
    else
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* tr = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
        f->setObservedTransitionRates( tr );
    }
    
    // set transition rates for hidden character states
    if ( this->args[1].getVariable()->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* tr = static_cast<const ModelVector<Real> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
        f->setHiddenTransitionRates( tr );
    }
    else if ( this->args[1].getVariable()->getRevObject().isType( RateGenerator::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RateGenerator >* tr = static_cast<const RateGenerator &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
        f->setHiddenTransitionRates( tr );
    }
    else
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* tr = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
        f->setHiddenTransitionRates( tr );
    }
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_hiddenStateRateMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        std::vector<TypeSpec> transition_rate_types;
        transition_rate_types.push_back( ModelVector<ModelVector<RealPos> >::getClassTypeSpec() );
        transition_rate_types.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        transition_rate_types.push_back( RateGenerator::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "observered_transition_rates", transition_rate_types,         "Transition rates between observed states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "hidden_transition_rates",     transition_rate_types,         "Transition rates between hidden states.",   ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rescaled",                    RlBoolean::getClassTypeSpec(), "Should the matrix be normalized?",          ArgumentRule::BY_VALUE,              ArgumentRule::ANY, new RlBoolean(true) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_hiddenStateRateMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_hiddenStateRateMatrix";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_hiddenStateRateMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_hiddenStateRateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnHiddenStateRateMatrix";
    
    return f_name;
}


const TypeSpec& Func_hiddenStateRateMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
