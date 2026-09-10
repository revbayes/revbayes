#include <iosfwd>
#include <string>
#include <vector>

#include "ChromosomesPloidyRateMatrixFunction.h"
#include "Func_chromosomesPloidy.h"
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
Func_chromosomesPloidy::Func_chromosomesPloidy( void ) : TypedFunction<RateMatrix>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_chromosomesPloidy* Func_chromosomesPloidy::clone( void ) const
{
    
    return new Func_chromosomesPloidy( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_chromosomesPloidy::createFunction( void ) const {
    
    
    RevBayesCore::TypedDagNode< std::int64_t >* n           = static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* gamma_d  = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* delta_d  = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* rho_d    = static_cast<const RealPos &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* eta_d    = static_cast<const RealPos &>( this->args[4].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* gamma_dl = static_cast<const Real &>( this->args[5].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* delta_dl = static_cast<const Real &>( this->args[6].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* gamma_p  = static_cast<const RealPos &>( this->args[7].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* delta_p  = static_cast<const RealPos &>( this->args[8].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* rho_p    = static_cast<const RealPos &>( this->args[9].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* eta_p    = static_cast<const RealPos &>( this->args[10].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* gamma_pl = static_cast<const Real &>( this->args[11].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* delta_pl = static_cast<const Real &>( this->args[12].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* beta     = static_cast<const Real &>( this->args[13].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::ChromosomesPloidyRateMatrixFunction* f = new RevBayesCore::ChromosomesPloidyRateMatrixFunction( n, gamma_d, delta_d, rho_d, eta_d, gamma_dl, delta_dl, gamma_p, delta_p, rho_p, eta_p, gamma_pl, delta_pl, beta);
    
	return f;
}


/* Get argument rules */
const ArgumentRules& Func_chromosomesPloidy::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
      
        argumentRules.push_back( new ArgumentRule( "maxChromosomes", Natural::getClassTypeSpec(), "Maximum number of chromosomes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "gamma_d"       , RealPos::getClassTypeSpec(), "Diploid rate of gain of a single chromosome.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "delta_d"       , RealPos::getClassTypeSpec(), "Diploid rate of loss of a single chromosome.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "rho_d"         , RealPos::getClassTypeSpec(), "Diploid rate of polyploidization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "eta_d"         , RealPos::getClassTypeSpec(), "Disploid rate of demi-polyploidization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "gamma_dl"      , Real::getClassTypeSpec(), "Diploid rate modifier for chromosome gains.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "delta_dl"      , Real::getClassTypeSpec(), "Diploid rate modifier for chromosome losses.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "gamma_p"       , RealPos::getClassTypeSpec(), "Polyploid rate of gain of a single chromosome.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "delta_p"       , RealPos::getClassTypeSpec(), "Polyploid rate of loss of a single chromosome.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "rho_p"         , RealPos::getClassTypeSpec(), "Polyploid rate of polyploidization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "eta_p"         , RealPos::getClassTypeSpec(), "Polyploid rate of demi-polyploidization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "gamma_pl"      , Real::getClassTypeSpec(), "Polyploid rate modifier for chromosome gains.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "delta_pl"      , Real::getClassTypeSpec(), "Polyploid rate modifier for chromosome losses.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "beta"          , Real::getClassTypeSpec(), "Rate of rediploidization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_chromosomesPloidy::getClassType(void)
{
    
    static std::string rev_type = "Func_chromosomesPloidy";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_chromosomesPloidy::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_chromosomesPloidy::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnChromosomesPloidy";
    
    return f_name;
}


const TypeSpec& Func_chromosomesPloidy::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
