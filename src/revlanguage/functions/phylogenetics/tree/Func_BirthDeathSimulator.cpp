#include "Argument.h"
#include "ArgumentRule.h"
#include "BirthDeathForwardSimulator.h"
#include "Func_BirthDeathSimulator.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RbException.h"
#include "RealPos.h"
#include "RlSimplex.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "RlUtils.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"

#include <math.h>       /* log2 */

using namespace RevLanguage;

/** Default constructor */
Func_BirthDeathSimulator::Func_BirthDeathSimulator( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_BirthDeathSimulator* Func_BirthDeathSimulator::clone( void ) const
{
    
    return new Func_BirthDeathSimulator( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_BirthDeathSimulator::execute( void )
{
    
//    int num_taxa                = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
//    const std::string& type     = static_cast<const RlString &>( args[1].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::BirthDeathForwardSimulator simulator;
    
    size_t arg_index = 0;
    
    const std::vector<double> &timeline = static_cast<const ModelVector<RealPos> &>( args[arg_index].getVariable()->getRevObject() ).getValue();
    simulator.setTimeline( timeline );
    
    ++arg_index;
    std::vector< std::vector<double> > speciation_rates;
    if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< ModelVector< RealPos > >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector< RevBayesCore::RbVector<double> > &tmp_speciation_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        for (size_t i=0; i<tmp_speciation_rates.size(); ++i)
        {
            speciation_rates.push_back( tmp_speciation_rates[i] );
        }
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< RealPos >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector<double> &tmp_speciation_rates = static_cast<const ModelVector<RealPos> &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        speciation_rates.push_back( tmp_speciation_rates );
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( RealPos::getClassTypeSpec() ) )
    {
        double tmp_speciation_rates = static_cast<const RealPos &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        speciation_rates.push_back( std::vector<double>(1, tmp_speciation_rates) );
    }
    simulator.setSpeciationRate( speciation_rates );
    
    ++arg_index;
    std::vector< std::vector<double> > extinction_rates;
    if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< ModelVector< RealPos > >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector< RevBayesCore::RbVector<double> > &tmp_extinction_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        for (size_t i=0; i<tmp_extinction_rates.size(); ++i)
        {
            extinction_rates.push_back( tmp_extinction_rates[i] );
        }
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< RealPos >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector<double> &tmp_extinction_rates = static_cast<const ModelVector<RealPos> &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        extinction_rates.push_back( tmp_extinction_rates );
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( RealPos::getClassTypeSpec() ) )
    {
        double tmp_extinction_rates = static_cast<const RealPos &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        extinction_rates.push_back( std::vector<double>(1, tmp_extinction_rates) );
    }
    simulator.setExtinctionRate( extinction_rates );
    
    
    ++arg_index;
    std::vector< std::vector<double> > sampling_rates;
    if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< ModelVector< RealPos > >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector< RevBayesCore::RbVector<double> > &tmp_sampling_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        for (size_t i=0; i<tmp_sampling_rates.size(); ++i)
        {
            sampling_rates.push_back( tmp_sampling_rates[i] );
        }
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< RealPos >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector<double> &tmp_sampling_rates = static_cast<const ModelVector<RealPos> &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        sampling_rates.push_back( tmp_sampling_rates );
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( RealPos::getClassTypeSpec() ) )
    {
        double tmp_sampling_rates = static_cast<const RealPos &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        sampling_rates.push_back( std::vector<double>(1, tmp_sampling_rates) );
    }
    simulator.setSamplingRate( sampling_rates );
    
    
    ++arg_index;
    std::vector< std::vector<double> > sampling_extinction_rates;
    if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< ModelVector< Probability > >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector< RevBayesCore::RbVector<double> > &tmp_sampling_extinction_rates = static_cast<const ModelVector< ModelVector<Probability> > &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        for (size_t i=0; i<tmp_sampling_extinction_rates.size(); ++i)
        {
            sampling_extinction_rates.push_back( tmp_sampling_extinction_rates[i] );
        }
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< Probability >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector<double> &tmp_sampling_extinction_rates = static_cast<const ModelVector<Probability> &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        sampling_extinction_rates.push_back( tmp_sampling_extinction_rates );
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( Probability::getClassTypeSpec() ) )
    {
        double tmp_sampling_extinction_rates = static_cast<const Probability &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        sampling_extinction_rates.push_back( std::vector<double>(1, tmp_sampling_extinction_rates) );
    }
    simulator.setSamplingExtinctionRate( sampling_extinction_rates );
    
    
    ++arg_index;
    std::vector< std::vector<double> > speciation_probability;
    if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< ModelVector< Probability > >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector< RevBayesCore::RbVector<double> > &tmp_speciation_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        for (size_t i=0; i<tmp_speciation_probs.size(); ++i)
        {
            speciation_probability.push_back( tmp_speciation_probs[i] );
        }
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< Probability >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector<double> &tmp_speciation_probs = static_cast<const ModelVector<Probability> &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        speciation_probability.push_back( tmp_speciation_probs );
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( Probability::getClassTypeSpec() ) )
    {
        double tmp_speciation_probs = static_cast<const Probability &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        speciation_probability.push_back( std::vector<double>(1, tmp_speciation_probs) );
    }
    simulator.setBurstProbability( speciation_probability );
    
    
    ++arg_index;
    std::vector< std::vector<double> > extinction_probability;
    if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< ModelVector< Probability > >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector< RevBayesCore::RbVector<double> > &tmp_extinction_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        for (size_t i=0; i<tmp_extinction_probs.size(); ++i)
        {
            extinction_probability.push_back( tmp_extinction_probs[i] );
        }
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< Probability >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector<double> &tmp_extinction_probs = static_cast<const ModelVector<Probability> &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        extinction_probability.push_back( tmp_extinction_probs );
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( Probability::getClassTypeSpec() ) )
    {
        double tmp_extinction_probs = static_cast<const Probability &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        extinction_probability.push_back( std::vector<double>(1, tmp_extinction_probs) );
    }
    simulator.setMassExtinctionProbability( extinction_probability );
    
    
    ++arg_index;
    std::vector< std::vector<double> > sampling_probability;
    if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< ModelVector< Probability > >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector< RevBayesCore::RbVector<double> > &tmp_sampling_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        for (size_t i=0; i<tmp_sampling_probs.size(); ++i)
        {
            sampling_probability.push_back( tmp_sampling_probs[i] );
        }
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< Probability >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector<double> &tmp_sampling_probs = static_cast<const ModelVector<Probability> &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        sampling_probability.push_back( tmp_sampling_probs );
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( Probability::getClassTypeSpec() ) )
    {
        double tmp_sampling_probs = static_cast<const Probability &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        sampling_probability.push_back( std::vector<double>(1, tmp_sampling_probs) );
    }
    simulator.setSamplingProbability( sampling_probability );
    
    
    ++arg_index;
    std::vector< std::vector<double> > sampling_extinction_probability;
    if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< ModelVector< Probability > >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector< RevBayesCore::RbVector<double> > &tmp_sampling_extinction_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        for (size_t i=0; i<tmp_sampling_extinction_probs.size(); ++i)
        {
            sampling_extinction_probability.push_back( tmp_sampling_extinction_probs[i] );
        }
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( ModelVector< Probability >::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector<double> &tmp_sampling_extinction_probs = static_cast<const ModelVector<Probability> &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        sampling_extinction_probability.push_back( tmp_sampling_extinction_probs );
    }
    else if ( args[arg_index].getVariable()->getRevObject().isType( Probability::getClassTypeSpec() ) )
    {
        double tmp_sampling_extinction_probs = static_cast<const Probability &>( args[arg_index].getVariable()->getRevObject() ).getValue();
        sampling_extinction_probability.push_back( std::vector<double>(1, tmp_sampling_extinction_probs) );
    }
    simulator.setSamplingExtinctionProbability( sampling_extinction_probability );
    
    
    ++arg_index;
    const std::vector<double> &root_probs = static_cast<const Simplex &>( args[arg_index].getVariable()->getRevObject() ).getValue();
    simulator.setRootCategoryProbabilities( root_probs );
    
    
    ++arg_index;
    double time = static_cast<const RealPos &>( args[arg_index].getVariable()->getRevObject() ).getValue();
    
    // the time tree object (topology + times)
    RevBayesCore::Tree *my_tree = simulator.simulateTreeConditionTime( time, RevBayesCore::BirthDeathForwardSimulator::ROOT );
    
    return new RevVariable( new TimeTree( my_tree ) );
}


/** Get argument rules */
const ArgumentRules& Func_BirthDeathSimulator::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        std::vector<TypeSpec> rate_options;
        rate_options.push_back( ModelVector< ModelVector< RealPos > >::getClassTypeSpec() );
        rate_options.push_back( ModelVector< RealPos >::getClassTypeSpec() );
        rate_options.push_back( RealPos::getClassTypeSpec() );
        
        std::vector<TypeSpec> prob_options;
        prob_options.push_back( ModelVector< ModelVector< Probability > >::getClassTypeSpec() );
        prob_options.push_back( ModelVector< Probability >::getClassTypeSpec() );
        prob_options.push_back( Probability::getClassTypeSpec() );
        
        argumentRules.push_back( new ArgumentRule( "timeline", ModelVector< RealPos >::getClassTypeSpec(), "The endpoints of the time intervals (episodes). You should include 0 at the end. We use ages before the present.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>( RevBayesCore::RbVector<double>(1,0) ) ) );
        argumentRules.push_back( new ArgumentRule( "lambda", rate_options, "The speciation rates for each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "mu", rate_options, "The extinction rates for each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos( 0.0 ) ) );
        argumentRules.push_back( new ArgumentRule( "phi", rate_options, "The sampling rates for each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos( 0.0 ) ) );
        argumentRules.push_back( new ArgumentRule( "r", prob_options, "The extinction probability when rate-sampling happens for each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability( 0.0 ) ) );

        argumentRules.push_back( new ArgumentRule( "Lambda", prob_options, "The burst probability at the end of each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability( 0.0 ) ) );
        argumentRules.push_back( new ArgumentRule( "Mu", prob_options, "The (mass) extinction probability at the end of each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability( 0.0 ) ) );
        argumentRules.push_back( new ArgumentRule( "Phi", prob_options, "The sampling probability at the end of each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability( 0.0 ) ) );
        argumentRules.push_back( new ArgumentRule( "R", prob_options, "The extinction probability when event-sampling happens for each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability( 0.0 ) ) );
        
        argumentRules.push_back( new ArgumentRule( "rootCategory", Simplex::getClassTypeSpec(), "The probabilities of the categories for the root.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Simplex( RevBayesCore::Simplex(1,1) ) ) );

        argumentRules.push_back( new ArgumentRule( "time", RealPos::getClassTypeSpec(), "The time/age before the present.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        
//        argumentRules.push_back( new ArgumentRule( "num_taxa", Natural::getClassTypeSpec(), "How many taxa this tree has.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
//
//        std::vector<std::string> optionsCondition;
//        optionsCondition.push_back( "balanced" );
//        optionsCondition.push_back( "caterpillar" );
//        //        optionsCondition.push_back( "random" );
//        argumentRules.push_back( new OptionRule( "type"    , new RlString("balanced"), optionsCondition, "The type of the shape of the topology." ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_BirthDeathSimulator::getClassType(void)
{
    
    static std::string rev_type = "Func_BirthDeathSimulator";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_BirthDeathSimulator::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_BirthDeathSimulator::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "simBirthDeath";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_BirthDeathSimulator::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_BirthDeathSimulator::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = TimeTree::getClassTypeSpec();
    
    return return_typeSpec;
}
