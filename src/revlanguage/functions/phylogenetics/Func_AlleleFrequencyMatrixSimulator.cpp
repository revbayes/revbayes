#include "Func_AlleleFrequencyMatrixSimulator.h"

#include <map>
#include <string>
#include <vector>

#include "AlleleFrequencySimulator.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlTaxon.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "HomologousDiscreteCharacterData.h"
#include "ModelObject.h"
#include "RbIterator.h"
#include "RbIteratorImpl.h"
#include "RbVector.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlMatrixRealPos.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/** default constructor */
Func_AlleleFrequencyMatrixSimulator::Func_AlleleFrequencyMatrixSimulator( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_AlleleFrequencyMatrixSimulator* Func_AlleleFrequencyMatrixSimulator::clone( void ) const
{
    
    return new Func_AlleleFrequencyMatrixSimulator( *this );
}


RevPtr<RevVariable> Func_AlleleFrequencyMatrixSimulator::execute()
{
    
    RevBayesCore::Tree* tree = NULL;
    std::vector<long> ps;
    long num_sites = 0;
    std::vector<long> samples;
    double root_branch = -1;
    
    long population_sizes                   = static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getValue();
    double generation_time                  = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getValue();
    bool moran_generations                  = static_cast<const RlBoolean &>( this->args[2].getVariable()->getRevObject() ).getValue();
    std::vector<double> mutation_rates      = static_cast<const ModelVector<RealPos> &>( this->args[3].getVariable()->getRevObject() ).getValue();
    double time                             = static_cast<const RealPos &>( this->args[4].getVariable()->getRevObject() ).getValue();
    long reps                               = static_cast<const Natural &>( this->args[5].getVariable()->getRevObject() ).getValue();


    RevBayesCore::AlleleFrequencySimulator sim = RevBayesCore::AlleleFrequencySimulator(tree, ps, generation_time, num_sites, mutation_rates, samples, root_branch, moran_generations );
    RevBayesCore::MatrixReal* m = sim.simulateAlleleFrequenciesMatrix( time, population_sizes, reps );
    
    return new RevVariable( new MatrixRealPos( m ) );
}


/* Get argument rules */
const ArgumentRules& Func_AlleleFrequencyMatrixSimulator::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "populationSize"    , Natural::getClassTypeSpec(), "The population size.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "generationTime"    , RealPos::getClassTypeSpec(), "The generation time for the simulations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "moranGenerations"  , RlBoolean::getClassTypeSpec(), "Is the generation time in Moran generation, i.e., scaled by N.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "mutationRates"     , ModelVector<RealPos>::getClassTypeSpec(), "The mutation rates from 0 to 1 and 1 to 0.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "time"              , RealPos::getClassTypeSpec(), "The time for the simulations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "reps"              , Natural::getClassTypeSpec(), "The number of replicate to simulate the frequencies.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argument_rules;
}


const std::string& Func_AlleleFrequencyMatrixSimulator::getClassType(void)
{
    
    static std::string rev_type = "Func_AlleleFrequencyMatrixSimulator";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_AlleleFrequencyMatrixSimulator::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_AlleleFrequencyMatrixSimulator::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "AlleleFrequencyMatrixSimulator";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_AlleleFrequencyMatrixSimulator::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RevNullObject::getClassTypeSpec();
    return return_typeSpec;
}


const TypeSpec& Func_AlleleFrequencyMatrixSimulator::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
