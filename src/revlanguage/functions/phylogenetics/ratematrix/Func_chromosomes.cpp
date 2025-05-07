#include <iosfwd>
#include <string>
#include <vector>

#include "ChromosomesRateMatrixFunction.h"
#include "Func_chromosomes.h"
#include "Natural.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RateGenerator.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

using namespace RevLanguage;

/** default constructor */
Func_chromosomes::Func_chromosomes( void ) : TypedFunction<RateMatrix>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_chromosomes* Func_chromosomes::clone( void ) const
{
    
    return new Func_chromosomes( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_chromosomes::createFunction( void ) const {
    
    
    RevBayesCore::TypedDagNode< std::int64_t >* n           = static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* gamma    = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* delta    = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* rho      = static_cast<const RealPos &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* eta      = static_cast<const RealPos &>( this->args[4].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* gamma_l  = static_cast<const Real &>( this->args[5].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* delta_l  = static_cast<const Real &>( this->args[6].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::ChromosomesRateMatrixFunction* f = new RevBayesCore::ChromosomesRateMatrixFunction( n, gamma, delta, rho, eta, gamma_l, delta_l );
    
	return f;
}


/* Get argument rules */
const ArgumentRules& Func_chromosomes::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
      
        argumentRules.push_back( new ArgumentRule( "maxChromosomes", Natural::getClassTypeSpec(), "Maximum number of chromosomes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "gamma"         , RealPos::getClassTypeSpec(), "Rate of gain of a single chromosome.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "delta"         , RealPos::getClassTypeSpec(), "Rate of loss of a single chromosome.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "rho"           , RealPos::getClassTypeSpec(), "Rate of polyploidization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "eta"           , RealPos::getClassTypeSpec(), "Rate of demi-polyploidization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "gamma_l"       , Real::getClassTypeSpec(), "Rate modifier for chromosome gains.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "delta_l"       , Real::getClassTypeSpec(), "Rate modifier for chromosome losses.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_chromosomes::getClassType(void)
{
    
    static std::string rev_type = "Func_chromosomes";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_chromosomes::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_chromosomes::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnChromosomes";
    
    return f_name;
}


const TypeSpec& Func_chromosomes::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
