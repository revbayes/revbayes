#include "Func_FlowT2Populations.h"
#include "FlowT2PopulationsRateMatrixFunction.h"
#include "RateMatrix_FlowT2Populations.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_FlowT2Populations::Func_FlowT2Populations( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_FlowT2Populations* Func_FlowT2Populations::clone( void ) const
{
    
    return new Func_FlowT2Populations( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_FlowT2Populations::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* n = static_cast<const ModelVector<RealPos> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* m = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* l = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* s = static_cast<const ModelVector<RealPos> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();

    if ( n->getValue().size() != 2 )
    {
        throw RbException("Make sure each vector has two entries.");
    }

    RevBayesCore::FlowT2PopulationsRateMatrixFunction* f = new RevBayesCore::FlowT2PopulationsRateMatrixFunction( n, m, l, s );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_FlowT2Populations::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "nu"          , ModelVector<Real>::getClassTypeSpec(), "The population sizes: nu=(nu_1,nu_2).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "mu"          , ModelVector<Real>::getClassTypeSpec(), "The mutation rates: mu=(mu_A,mu_a).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lambda"      , ModelVector<Real>::getClassTypeSpec(), "The allele flow rates: lambda=(lambda_12,lambda_21).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "sigma"       , ModelVector<Real>::getClassTypeSpec(), "The selection coefficients: sigma=(sigma_A,sigma_a).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_FlowT2Populations::getClassType(void)
{
    
    static std::string rev_type = "Func_FlowT2Populations";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_FlowT2Populations::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_FlowT2Populations::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnFlowT2Populations";
    
    return f_name;
}


const TypeSpec& Func_FlowT2Populations::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



