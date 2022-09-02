#include "AlleleFrequencySimulator.h"

#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Delimiter.h"
#include "RevObject.h"
#include "Real.h"
#include "RlAlleleFrequencySimulator.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlMatrixRealPos.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevLanguage { class Argument; }


using namespace RevLanguage;

AlleleFrequencySimulator::AlleleFrequencySimulator() : WorkspaceToCoreWrapperObject<RevBayesCore::AlleleFrequencySimulator>()
{

    // simulating a transition probability matrix
    
    ArgumentRules* tpm_arg_rules = new ArgumentRules();
    
    tpm_arg_rules->push_back( new ArgumentRule( "populationSize"    , Natural::getClassTypeSpec(), "The population size.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    tpm_arg_rules->push_back( new ArgumentRule( "time"              , RealPos::getClassTypeSpec(), "The time for the simulations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    tpm_arg_rules->push_back( new ArgumentRule( "reps"              , Natural::getClassTypeSpec(), "The number of replicate to simulate the frequencies.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

    methods.addFunction(new MemberProcedure( "transitionProbabilityMatrix", MatrixRealPos::getClassTypeSpec(), tpm_arg_rules) );

    
    
    // simulating a transition probability vector
    
    ArgumentRules* tpv_arg_rules = new ArgumentRules();
    
    tpv_arg_rules->push_back( new ArgumentRule( "populationSize"    , Natural::getClassTypeSpec(), "The population size.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    tpv_arg_rules->push_back( new ArgumentRule( "startState"        , Natural::getClassTypeSpec(), "The start state for the transition probability simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    tpv_arg_rules->push_back( new ArgumentRule( "time"              , RealPos::getClassTypeSpec(), "The time for the simulations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    tpv_arg_rules->push_back( new ArgumentRule( "reps"              , Natural::getClassTypeSpec(), "The number of replicate to simulate the frequencies.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

    methods.addFunction(new MemberProcedure( "transitionProbabilityVector", ModelVector< Probability >::getClassTypeSpec(), tpv_arg_rules) );

    
    
    // simulating a transition probability matrix
    
    ArgumentRules* epoch_arg_rules = new ArgumentRules();
    
    epoch_arg_rules->push_back( new ArgumentRule( "populationSize"    , ModelVector<Natural>::getClassTypeSpec(), "The vector of population sizes per epoch.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    epoch_arg_rules->push_back( new ArgumentRule( "startState"        , Natural::getClassTypeSpec(), "The start state for the allele frequency simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    epoch_arg_rules->push_back( new ArgumentRule( "sampleSize"        , Natural::getClassTypeSpec(), "The sample size, i.e., number of individuals, at present.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    epoch_arg_rules->push_back( new ArgumentRule( "time"              , ModelVector<RealPos>::getClassTypeSpec(), "The vector of epochs times for the simulations. The process starts at time 0 and goes forward to the present.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    epoch_arg_rules->push_back( new ArgumentRule( "reps"              , Natural::getClassTypeSpec(), "The number of replicate to simulate the frequencies.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

    methods.addFunction(new MemberProcedure( "epochSFS", ModelVector< Probability >::getClassTypeSpec(), epoch_arg_rules) );

    
    
    // simulating a data matrix
    
    ArgumentRules* data_matrix_arg_rules = new ArgumentRules();
    
    data_matrix_arg_rules->push_back( new ArgumentRule( "tree"              , TimeTree::getClassTypeSpec(), "The tree along which we want to simulate.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    data_matrix_arg_rules->push_back( new ArgumentRule( "populationSizes"   , ModelVector<Natural>::getClassTypeSpec(), "The population sizes for all branches, including the root branch.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    data_matrix_arg_rules->push_back( new ArgumentRule( "numSites"          , Natural::getClassTypeSpec(), "The number of sites to simulate.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    data_matrix_arg_rules->push_back( new ArgumentRule( "samplesPerSpecies" , ModelVector<Natural>::getClassTypeSpec(), "The observed number of samples per species.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    data_matrix_arg_rules->push_back( new ArgumentRule( "rootBranch"        , RealPos::getClassTypeSpec(), "The length of the root branch for simulating polymorphisms at the root.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    data_matrix_arg_rules->push_back( new ArgumentRule( "variable"          , RlBoolean::getClassTypeSpec(), "Do we condition on only observing variable sites?", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    data_matrix_arg_rules->push_back( new ArgumentRule( "filename"          , RlString::getClassTypeSpec(), "The filename for the counts file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

    methods.addFunction(new MemberProcedure( "dataMatrix", RlUtils::Void, data_matrix_arg_rules) );
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
AlleleFrequencySimulator* AlleleFrequencySimulator::clone(void) const
{
    
    return new AlleleFrequencySimulator(*this);
}


void AlleleFrequencySimulator::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // get the parameter values
    double gt               = static_cast<const RealPos &>( generation_time->getRevObject() ).getValue();
    bool mg                 = static_cast<const RlBoolean &>( moran_generations->getRevObject() ).getValue();
    std::vector<double> mr  = static_cast<const ModelVector<RealPos> &>( mutation_rates->getRevObject() ).getValue();
    
    value = new RevBayesCore::AlleleFrequencySimulator(gt, mr, mg);

    
}


/* Map calls to member methods */
RevPtr<RevVariable> AlleleFrequencySimulator::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "dataMatrix")
    {
        found = true;
        
        RevBayesCore::Tree* tree                        = static_cast<const TimeTree&>( args[0].getVariable()->getRevObject() ).getValue().clone();
        const std::vector<long>& population_sizes       = static_cast<const ModelVector<Natural> &>( args[1].getVariable()->getRevObject() ).getValue();
        long num_sites                                  = static_cast<const Natural &>( args[2].getVariable()->getRevObject() ).getValue();
        const std::vector<long>& samples_per_species    = static_cast<const ModelVector<Natural> &>( args[3].getVariable()->getRevObject() ).getValue();
        double root_branch                              = static_cast<const RealPos &>( args[4].getVariable()->getRevObject() ).getValue();
        bool variable                                   = static_cast<const RlBoolean &>( args[5].getVariable()->getRevObject() ).getValue();
        const std::string& fn                           = static_cast<const RlString &>( args[6].getVariable()->getRevObject() ).getValue();

        value->simulateAlleleFrequencies(tree, population_sizes, num_sites, samples_per_species, root_branch, fn, variable);
        
        delete tree;
        
    }
    else if (name == "transitionProbabilityMatrix")
    {
        found = true;
        
        long population_sizes                   = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        double time                             = static_cast<const RealPos &>( args[1].getVariable()->getRevObject() ).getValue();
        long reps                               = static_cast<const Natural &>( args[2].getVariable()->getRevObject() ).getValue();

        RevBayesCore::MatrixReal* m = value->simulateAlleleFrequenciesMatrix(time, population_sizes, reps);
        
        return new RevVariable( new MatrixRealPos( m ) );
    }
    else if (name == "transitionProbabilityVector")
    {
        found = true;
            
        long population_sizes                   = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        long start                              = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue() - 1;
        double time                             = static_cast<const RealPos &>( args[2].getVariable()->getRevObject() ).getValue();
        long reps                               = static_cast<const Natural &>( args[3].getVariable()->getRevObject() ).getValue();

        RevBayesCore::RbVector<double>* m = value->simulateAlleleFrequenciesVector(time, population_sizes, reps, start);
        
        return new RevVariable( new ModelVector<Probability>( *m ) );
    }
    else if (name == "epochSFS")
    {
        found = true;
        
        const std::vector<long>& population_sizes   = static_cast<const ModelVector<Natural> &>( args[0].getVariable()->getRevObject() ).getValue();
        long start                                  = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue() - 1;
        long sample_size                            = static_cast<const Natural &>( args[2].getVariable()->getRevObject() ).getValue();
        const std::vector<double>& times            = static_cast<const ModelVector<RealPos> &>( args[3].getVariable()->getRevObject() ).getValue();
        long reps                                   = static_cast<const Natural &>( args[4].getVariable()->getRevObject() ).getValue();

        RevBayesCore::RbVector<double>* m = value->simulateAlleleFrequenciesVectorEpoch(times, population_sizes, reps, start, sample_size);
        
        return new RevVariable( new ModelVector<Probability>( *m ) );
    }
    
    return WorkspaceToCoreWrapperObject<RevBayesCore::AlleleFrequencySimulator>::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& AlleleFrequencySimulator::getClassType(void)
{
    
    static std::string rev_type = "AlleleFrequencySimulator";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& AlleleFrequencySimulator::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::AlleleFrequencySimulator>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string AlleleFrequencySimulator::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "AlleleFrequencySimulator";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& AlleleFrequencySimulator::getParameterRules(void) const
{
    
    static MemberRules argument_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "generationTime"    , RealPos::getClassTypeSpec(), "The generation time for the simulations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "moranGenerations"  , RlBoolean::getClassTypeSpec(), "Is the generation time in Moran generation, i.e., scaled by N.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "mutationRates"     , ModelVector<RealPos>::getClassTypeSpec(), "The mutation rates from 0 to 1 and 1 to 0.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argument_rules;
}


/** Get type spec */
const TypeSpec& AlleleFrequencySimulator::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void AlleleFrequencySimulator::printValue(std::ostream &o) const
{
    
    o << "AlleleFrequencySimulator";
}


/** Set a member variable */
void AlleleFrequencySimulator::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "generationTime")
    {
        generation_time = var;
    }
    else if ( name == "moranGenerations")
    {
        moran_generations = var;
    }
    else if ( name == "mutationRates")
    {
        mutation_rates = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}
