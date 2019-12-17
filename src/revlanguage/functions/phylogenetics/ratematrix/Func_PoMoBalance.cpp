#include "Func_PoMoBalance.h"
#include "PoMoBalanceRateMatrixFunction.h"
#include "RateMatrix_PoMoBalance.h"
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
Func_PoMoBalance::Func_PoMoBalance( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PoMoBalance* Func_PoMoBalance::clone( void ) const
{
    
    return new Func_PoMoBalance( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_PoMoBalance::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< double                        >* n = static_cast<const RealPos              &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex          >* p = static_cast<const Simplex              &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* r = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* s = static_cast<const ModelVector<RealPos> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double                        >* b = static_cast<const RealPos              &>( this->args[4].getVariable()->getRevObject() ).getDagNode();

    if ( ((p->getValue().size())*(p->getValue().size())-(p->getValue().size()))/2 !=  r->getValue().size() )
    {
        throw RbException("The number of alleles and exchangeabilities do not match.");
    }

    RevBayesCore::PoMoBalanceRateMatrixFunction* f = new RevBayesCore::PoMoBalanceRateMatrixFunction( n, p, r, s, b );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_PoMoBalance::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "N"          , RealPos::getClassTypeSpec(), "Effective population size.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "pi"         , Simplex::getClassTypeSpec(), "Allele base frequencies: pi=(pi_A,pi_C,pi_G,pi_T).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rho"        , ModelVector<Real>::getClassTypeSpec(), "Exchangeabilities: rho=(rho_AC,rho_AG,rho_AT,rho_CG,rho_CT,rho_GT).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "sigma"      , ModelVector<Real>::getClassTypeSpec(), "Selection coefficients: sigma=(sigma_A,sigma_C,sigma_G,sigma_T).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "beta"       , RealPos::getClassTypeSpec(), "Balancing selection coefficient.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_PoMoBalance::getClassType(void)
{
    
    static std::string rev_type = "Func_PoMoBalance";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PoMoBalance::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PoMoBalance::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPoMoBalance";
    
    return f_name;
}


const TypeSpec& Func_PoMoBalance::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



