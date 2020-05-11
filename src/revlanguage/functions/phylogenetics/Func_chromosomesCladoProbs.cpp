#include "Func_chromosomesCladoProbs.h"

#include "CladogeneticProbabilityMatrix.h"
#include "ChromosomesCladogenicStateFunction.h"
#include "Natural.h"
#include "RlCladogeneticProbabilityMatrix.h"
#include "RlDeterministicNode.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Simplex; }

using namespace RevLanguage;

/** default constructor */
Func_chromosomesCladoProbs::Func_chromosomesCladoProbs( void ) : TypedFunction<CladogeneticProbabilityMatrix>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_chromosomesCladoProbs* Func_chromosomesCladoProbs::clone( void ) const {
    
    return new Func_chromosomesCladoProbs( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CladogeneticProbabilityMatrix >* Func_chromosomesCladoProbs::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* ep = static_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    int ns = (int)static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::ChromosomesCladogenicStateFunction* f = new RevBayesCore::ChromosomesCladogenicStateFunction( ep, ns );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_chromosomesCladoProbs::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "eventProbs", Simplex::getClassTypeSpec(), "The probabilities of the different event types.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "maxChromosomes", Natural::getClassTypeSpec(), "The maximum number of chromosomes.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_chromosomesCladoProbs::getClassType(void)
{
    
    static std::string rev_type = "Func_chromosomesCladoProbs";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_chromosomesCladoProbs::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_chromosomesCladoProbs::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnChromosomesCladoProbs";
    
    return f_name;
}


const TypeSpec& Func_chromosomesCladoProbs::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
