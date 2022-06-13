#include "Func_revPoMoM2N.h"
#include "Func_revPoMoM2N.h"
#include "revPoMoM2NRateMatrixFunction.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RateMatrix_revPoMoM2N.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_revPoMoM2N::Func_revPoMoM2N( void ) : TypedFunction<RateMatrix>( )
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_revPoMoM2N* Func_revPoMoM2N::clone( void ) const
{

    return new Func_revPoMoM2N( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_revPoMoM2N::createFunction( void ) const
{
    long                                                         v  = static_cast<const Natural               &>( this->args[0].getVariable()->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode< double >*                        n  = static_cast<const RealPos               &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* m  = static_cast<const ModelVector<RealPos>  &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >*                        g  = static_cast<const RealPos               &>( this->args[3].getVariable()->getRevObject() ).getDagNode();

    if ( m->getValue().size() !=  2 )
    {
        throw RbException("The number of mutation rates does not match the number of pairwise combinations of alleles: 2.");
    }
    
    RevBayesCore::revPoMoM2NRateMatrixFunction* f = new RevBayesCore::revPoMoM2NRateMatrixFunction( v, n, m, g );

    return f;
}


/* Get argument rules */
const ArgumentRules& Func_revPoMoM2N::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {

        argumentRules.push_back( new ArgumentRule( "M"          , Natural::getClassTypeSpec(), "Virtual population size", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "N"          , RealPos::getClassTypeSpec(), "Effective Population size", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "mu"         , ModelVector<RealPos>::getClassTypeSpec(), "Vector of mutation rates: mu=(mu_a0a1,mu_a1a0)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "g"          , RealPos::getClassTypeSpec(), "Generation time", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_revPoMoM2N::getClassType(void)
{

    static std::string rev_type = "Func_revPoMoM2N";

	return rev_type;

}


/* Get class type spec describing type of object */
const TypeSpec& Func_revPoMoM2N::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;

}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_revPoMoM2N::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnReversiblePoMoM2N";

    return f_name;
}


const TypeSpec& Func_revPoMoM2N::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
