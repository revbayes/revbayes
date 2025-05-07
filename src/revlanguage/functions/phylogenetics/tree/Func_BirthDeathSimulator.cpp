#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "BirthDeathForwardSimulator.h"
#include "Func_BirthDeathSimulator.h"
#include "ModelVector.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlSimplex.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "Procedure.h"
#include "OptionRule.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "RlString.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Tree; }

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
    else
    {
        sampling_extinction_probability = sampling_extinction_rates;
    }

    simulator.setSamplingExtinctionProbability( sampling_extinction_probability );

    ++arg_index;
    const std::vector<double> &root_probs = static_cast<const Simplex &>( args[arg_index].getVariable()->getRevObject() ).getValue();
    simulator.setRootCategoryProbabilities( root_probs );

    ++arg_index;
    double time = static_cast<const RealPos &>( args[arg_index].getVariable()->getRevObject() ).getValue();
    
    ++arg_index;
    RevBayesCore::BirthDeathForwardSimulator::SIM_CONDITION cdt = RevBayesCore::BirthDeathForwardSimulator::TIME;
    const std::string& cdt_str = static_cast<const RlString &>( args[arg_index].getVariable()->getRevObject() ).getValue();

    ++arg_index;
    long max_lineages = static_cast<const Natural &>( args[arg_index].getVariable()->getRevObject() ).getValue();
    simulator.setMaxNumLineages(max_lineages);

    ++arg_index;
    bool complete_tree = static_cast<const RlBoolean &>( args[arg_index].getVariable()->getRevObject() ).getValue();
    simulator.setCompleteTree(complete_tree);

    if ( cdt_str == "time" )
    {
        cdt = RevBayesCore::BirthDeathForwardSimulator::TIME;
    }
    else if ( cdt_str == "root" )
    {
        cdt = RevBayesCore::BirthDeathForwardSimulator::ROOT;
    }
    else if ( cdt_str == "survival" )
    {
        cdt = RevBayesCore::BirthDeathForwardSimulator::SURVIVAL;
    }

    // the time tree object (topology + times)
    RevBayesCore::Tree *my_tree = simulator.simulateTreeConditionTime( time, cdt );

    return new RevVariable( new TimeTree( my_tree ) );
}


/** Get argument rules */
const ArgumentRules& Func_BirthDeathSimulator::getArgumentRules( void ) const
{

    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;

    if ( rules_set == false )
    {
        std::vector<TypeSpec> rate_options;
        rate_options.push_back( ModelVector< ModelVector< RealPos > >::getClassTypeSpec() );
        rate_options.push_back( ModelVector< RealPos >::getClassTypeSpec() );
        rate_options.push_back( RealPos::getClassTypeSpec() );

        std::vector<TypeSpec> prob_options;
        prob_options.push_back( ModelVector< ModelVector< Probability > >::getClassTypeSpec() );
        prob_options.push_back( ModelVector< Probability >::getClassTypeSpec() );
        prob_options.push_back( Probability::getClassTypeSpec() );

        argument_rules.push_back( new ArgumentRule( "timeline", ModelVector< RealPos >::getClassTypeSpec(), "The endpoints of the time intervals (episodes). You should include 0 at the end. We use ages before the present.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>( RevBayesCore::RbVector<double>(1,0) ) ) );
        argument_rules.push_back( new ArgumentRule( "lambda", rate_options, "The speciation rates for each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "mu", rate_options, "The extinction rates for each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos( 0.0 ) ) );
        argument_rules.push_back( new ArgumentRule( "phi", rate_options, "The sampling rates for each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos( 0.0 ) ) );
        argument_rules.push_back( new ArgumentRule( "r", prob_options, "The extinction probability when rate-sampling happens for each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability( 0.0 ) ) );

        argument_rules.push_back( new ArgumentRule( "Lambda", prob_options, "The burst probability at the end of each interval (first value is ignored).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability( 0.0 ) ) );
        argument_rules.push_back( new ArgumentRule( "Mu", prob_options, "The (mass) extinction probability at the end of each interval (first value is ignored).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability( 0.0 ) ) );
        argument_rules.push_back( new ArgumentRule( "Phi", prob_options, "The sampling probability at the end of each interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability( 0.0 ) ) );
        argument_rules.push_back( new ArgumentRule( "R", prob_options, "The extinction probability when event-sampling happens for each interval (first value is ignored). If NULL, r is used instead.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        argument_rules.push_back( new ArgumentRule( "rootCategory", Simplex::getClassTypeSpec(), "The probabilities of the categories for the root.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Simplex( RevBayesCore::Simplex(1,1) ) ) );

        argument_rules.push_back( new ArgumentRule( "time", RealPos::getClassTypeSpec(), "The time/age before the present.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        std::vector<std::string> options_condition;
        options_condition.push_back( "time" );
        options_condition.push_back( "root" );
        options_condition.push_back( "survival" );
        argument_rules.push_back( new OptionRule( "condition", new RlString("root"), options_condition, "What outcome should we condition on?" ) );

        argument_rules.push_back( new ArgumentRule( "maxNumLineages", Natural::getClassTypeSpec(), "The maximum number of lineages allowed by the simulator. Simulations that reach this size will be aborted and re-started.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Natural( 100000 ) ) );
        argument_rules.push_back( new ArgumentRule( "completeTree", RlBoolean::getClassTypeSpec(), "Should the tree include all lineages, even those that went extinct?", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean( false ) ) );

        rules_set = true;
    }

    return argument_rules;
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
    std::string f_name = "simForwardBirthDeath";

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
