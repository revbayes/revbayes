#include "Func_revPoMoBalance4N.h"
#include "revPoMoBalance4NRateMatrixFunction.h"
#include "RateMatrix_revPoMoBalance4N.h"
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
Func_revPoMoBalance4N::Func_revPoMoBalance4N( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_revPoMoBalance4N* Func_revPoMoBalance4N::clone( void ) const
{
    
    return new Func_revPoMoBalance4N( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_revPoMoBalance4N::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< long                          >* n  = static_cast<const Natural              &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex          >* p  = static_cast<const Simplex              &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* r  = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* s  = static_cast<const ModelVector<Real>    &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* b  = static_cast<const ModelVector<RealPos> &>( this->args[4].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<long>   >* Bf = static_cast<const ModelVector<Natural> &>( this->args[5].getVariable()->getRevObject() ).getDagNode();

    if ( p->getValue().size() !=  4 )
    {
        throw RbException("The number of base frequencies pi does not match the number of alleles: 4.");
    }

    if ( r->getValue().size() !=  6 )
    {
        throw RbException("The number of exchangeabilities rho does not match the number of pairwise combinations of alleles: 6.");
    }

    if ( s->getValue().size() !=  4 )
    {
        throw RbException("The number of fitness coefficients phi does not match the number of alleles: 4.");
    }

    if ( b->getValue().size() !=  6 )
    {
        throw RbException("The number of balacing selection coefficients beta does not match the number of pairwise combinations of alleles: 6.");
    }

    if ( Bf->getValue().size() !=  6 )
    {
        throw RbException("The number of perfered frequencies B does not match the number of pairwise combinations of alleles: 6.");
    }

    RevBayesCore::revPoMoBalance4NRateMatrixFunction* f = new RevBayesCore::revPoMoBalance4NRateMatrixFunction( n, p, r, s, b, Bf );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_revPoMoBalance4N::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "N"          , Natural::getClassTypeSpec(), "Population size", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "pi"         , Simplex::getClassTypeSpec(), "Allele base frequencies: pi=(pi_A,pi_C,pi_G,pi_T)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rho"        , ModelVector<RealPos>::getClassTypeSpec(), "Exchangeabilities: rho=(rho_AC,rho_AG,rho_AT,rho_CG,rho_CT,rho_GT)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "phi"        , ModelVector<Real>   ::getClassTypeSpec(), "Fitness coefficients: phi=(phi_A,phi_C,phi_G,phi_T)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "beta"       , ModelVector<RealPos>::getClassTypeSpec(), "Balancing selection coefficients: beta=(beta_AC,beta_AG,beta_AT,beta_CG,beta_CT,beta_GT)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "B"          , ModelVector<Natural>::getClassTypeSpec(), "Preferred frequencies: B=(B_Ac,B_Ag,B_At,B_Cg,B_Ct,B_Gt)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_revPoMoBalance4N::getClassType(void)
{
    
    static std::string rev_type = "Func_revPoMoBalance4N";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_revPoMoBalance4N::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_revPoMoBalance4N::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnReversiblePoMoBalance4N";
    
    return f_name;
}


const TypeSpec& Func_revPoMoBalance4N::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



