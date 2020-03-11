#include <stddef.h>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_GLHBDSP.h"
#include "GeneralizedLineageHeterogeneousBirthDeathSamplingProcess.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RlCladogeneticProbabilityMatrix.h"
#include "RlRateGenerator.h"
#include "RealPos.h"
#include "RlSimplex.h"
#include "RlStochasticMatrix.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_GLHBDSP::Dist_GLHBDSP() : TypedDistribution<TimeTree>()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_GLHBDSP* Dist_GLHBDSP::clone( void ) const
{
    return new Dist_GLHBDSP(*this);
}


/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::TypedDistribution<RevBayesCore::Tree>* Dist_GLHBDSP::createDistribution( void ) const
{

	std::cout << "Calling constructor for the GLHBDSP." << std::endl;

    // the start age
    RevBayesCore::TypedDagNode<double>* sa = static_cast<const RealPos &>( start_age->getRevObject() ).getDagNode();

    // the start condition
    bool uo = ( start_condition == "originAge" ? true : false );

    // sampling condition
    const std::string& cond = static_cast<const RlString &>( condition->getRevObject() ).getValue();

    // regular parameters
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* l_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( lambda->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           l_times = static_cast<const ModelVector< RealPos> &>( lambda_times->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* m_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( mu->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           m_times = static_cast<const ModelVector< RealPos> &>( mu_times->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* p_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( phi->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           p_times = static_cast<const ModelVector< RealPos> &>( phi_times->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* d_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( delta->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           d_times = static_cast<const ModelVector< RealPos> &>( delta_times->getRevObject() ).getDagNode();

    // mass-event parameters
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* u_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( upsilon->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           u_times = static_cast<const ModelVector< RealPos> &>( upsilon_times->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* g_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( gamma->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           g_times = static_cast<const ModelVector< RealPos> &>( gamma_times->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* r_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( rho->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           r_times = static_cast<const ModelVector< RealPos> &>( rho_times->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* x_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( xi->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           x_times = static_cast<const ModelVector< RealPos> &>( xi_times->getRevObject() ).getDagNode();

    // state-change events
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RateGenerator > >*                 h_mats  = static_cast<const ModelVector< RateGenerator > &>( eta->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                                        h_times = static_cast<const ModelVector< RealPos> &>( eta_times->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::CladogeneticProbabilityMatrix > >* w_mats  = static_cast<const ModelVector< CladogeneticProbabilityMatrix > &>( omega->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                                        w_times = static_cast<const ModelVector< RealPos> &>( omega_times->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::MatrixReal > >*                    z_mats  = static_cast<const ModelVector< StochasticMatrix > &>( zeta->getRevObject() ).getDagNode();

    // root frequency
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex >* rf = static_cast<const Simplex &>( root_frequencies->getRevObject() ).getDagNode();;

    // checks for safety
    std::cout << "Checking argument sizes." << std::endl;
    if ( l_rates->getValue().size() != l_times->getValue().size() + 1 )
    {
    	throw RbException( "Number of lambda vectors does not match the number of intervals." );
    }
    if ( m_rates->getValue().size() != m_times->getValue().size() + 1 )
    {
    	throw RbException( "Number of mu vectors does not match the number of intervals." );
    }
    if ( p_rates->getValue().size() != p_times->getValue().size() + 1 )
    {
    	throw RbException( "Number of phi vectors does not match the number of intervals." );
    }
    if ( d_rates->getValue().size() != d_times->getValue().size() + 1 )
    {
    	throw RbException( "Number of delta vectors does not match the number of intervals." );
    }
    if ( u_probs->getValue().size() != u_times->getValue().size() )
    {
    	throw RbException( "Number of mass-speciation events does not match the number of times." );
    }
    if ( g_probs->getValue().size() != g_times->getValue().size() )
    {
    	throw RbException( "Number of mass-extinction events does not match the number of times." );
    }
    if ( z_mats->getValue().size() != g_times->getValue().size() )
    {
    	throw RbException( "Number of mass-extinction associated state-change matrices does not match the number of times." );
    }
    if ( r_probs->getValue().size() != r_times->getValue().size() )
    {
    	throw RbException( "Number of mass-sampling events does not match the number of times." );
    }
    if ( x_probs->getValue().size() != x_times->getValue().size() )
    {
    	throw RbException( "Number of mass-destructive-sampling events does not match the number of times." );
    }
    if ( h_mats->getValue().size() != h_times->getValue().size() + 1 )
    {
    	throw RbException( "Number of eta matrices does not match the number of intervals." );
    }
    if ( w_mats->getValue().size() != w_times->getValue().size() + 1 )
    {
    	throw RbException( "Number of omega matrices does not match the number of intervals." );
    }

    // make the distribution
    std::cout << "Creating the distribution object." << std::endl;
    RevBayesCore::GeneralizedLineageHeterogeneousBirthDeathSamplingProcess* d = new RevBayesCore::GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(
    	sa, cond, rf,
		l_rates, l_times,
		m_rates, m_times,
		p_rates, p_times,
		d_rates, d_times,
		u_probs, u_times,
		g_probs, g_times,
		r_probs, r_times,
		x_probs, x_times,
		h_mats,  h_times,
		w_mats,  w_times,
		z_mats,
		uo);

    // return the distribution
    std::cout << "Done." << std::endl;
    return d;

}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_GLHBDSP::getClassType( void )
{

    static std::string rev_type = "Dist_GLHBDSP";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_GLHBDSP::getClassTypeSpec( void )
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<TimeTree>::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_GLHBDSP::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "GLHBDSP" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_GLHBDSP::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "GeneralizedLineageHeterogeneousBirthDeathProcess";

    return d_name;
}

MethodTable Dist_GLHBDSP::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution<TimeTree>::getDistributionMethods();

    return methods;
}



/**
 * Get the member rules used to create the constructor of this object.
 * \return The member rules.
 */
const MemberRules& Dist_GLHBDSP::getParameterRules(void) const
{

    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {
    	// the start of the process
    	std::vector<std::string> aliases;
        aliases.push_back("rootAge");
        aliases.push_back("originAge");
        dist_member_rules.push_back( new ArgumentRule( aliases, RealPos::getClassTypeSpec(), "The start time of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        // regular events (and times)
        dist_member_rules.push_back( new ArgumentRule( "lambda",       ModelVector< ModelVector<RealPos> >::getClassTypeSpec(),         "The vector of speciation rates for each time interval.",                                     ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "lambdaTimes",  ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which speciation rates change.",                                                ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );
        dist_member_rules.push_back( new ArgumentRule( "mu",           ModelVector< ModelVector<RealPos> >::getClassTypeSpec(),         "The vector of extinction rates for each time interval.",                                     ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "muTimes",      ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which extinction rates change.",                                                ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );
        dist_member_rules.push_back( new ArgumentRule( "phi",          ModelVector< ModelVector<RealPos> >::getClassTypeSpec(),         "The vector of sampling rates for each time interval.",                                       ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "phiTimes",     ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which sampling rates change.",                                                  ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );
        dist_member_rules.push_back( new ArgumentRule( "delta",        ModelVector< ModelVector<RealPos> >::getClassTypeSpec(),         "The vector of destructive-sampling rates for each time interval.",                           ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "deltaTimes",   ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which destructive-sampling rates change.",                                      ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );

        // mass events
        dist_member_rules.push_back( new ArgumentRule( "upsilon",      ModelVector< ModelVector<Probability> >::getClassTypeSpec(),     "The vector of speciation probabilities for each mass-speciation event.",                     ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "upsilonTimes", ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which mass-speciation events occur.",                                           ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "gamma",        ModelVector< ModelVector<Probability> >::getClassTypeSpec(),     "The vector of extinction probabilities for each mass-extinction event.",                     ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "gammaTimes",   ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which mass-extinction events occur.",                                           ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "rho",          ModelVector< ModelVector<Probability> >::getClassTypeSpec(),     "The vector of sampling probabilities for each mass-sampling event.",                         ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "rhoTimes",     ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which mass-sampling events occur.",                                             ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "xi",           ModelVector< ModelVector<Probability> >::getClassTypeSpec(),     "The vector of destructive-sampling probabilities for each mass-destructive-sampling event.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "xiTimes",      ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which mass-destructive-sampling events occur.",                                 ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        // state changes
        dist_member_rules.push_back( new ArgumentRule( "eta",          ModelVector< RateGenerator >::getClassTypeSpec(),                "The anagenetic rates of change for each time interval.",                                     ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "etaTimes",     ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which the anagenetic rates change.",                                            ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );
        dist_member_rules.push_back( new ArgumentRule( "omega",        ModelVector< CladogeneticProbabilityMatrix>::getClassTypeSpec(), "The cladogenetic event probabilities for each time interval.",                               ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "omegaTimes",   ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which the cladogenetic rates change.",                                          ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );
        dist_member_rules.push_back( new ArgumentRule( "zeta",         ModelVector<StochasticMatrix>::getClassTypeSpec(),               "The probabilities of change for each mass-extinction event.",                                ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        // the root frequency
        dist_member_rules.push_back( new ArgumentRule( "rootFreq",     Simplex::getClassTypeSpec(),                                     "Frequencies of each state at the beginning of the process.",                                 ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        // conditioning
		std::vector<std::string> options_condition;
		options_condition.push_back( "time" );
		options_condition.push_back( "survival" );
		dist_member_rules.push_back( new OptionRule( "condition", new RlString("time"), options_condition, "The condition of the process." ) );

        rules_set = true;
    }

    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_GLHBDSP::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();

    return ts;
}


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Dist_GLHBDSP::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

	if ( name == "rootAge" || name == "originAge" )
	{
		start_condition = name;
		start_age = var;
	}
	else if ( name == "lambda" )
	{
		lambda = var;
	}
	else if ( name == "lambdaTimes" )
	{
		lambda_times = var;
	}
	else if ( name == "mu" )
	{
		mu = var;
	}
	else if ( name == "muTimes" )
	{
		mu_times = var;
	}
	else if ( name == "phi" )
	{
		phi = var;
	}
	else if ( name == "phiTimes" )
	{
		phi_times = var;
	}
	else if ( name == "delta" )
	{
		delta = var;
	}
	else if ( name == "deltaTimes" )
	{
		delta_times = var;
	}
	else if ( name == "upsilon" )
	{
		upsilon = var;
	}
	else if ( name == "upsilonTimes" )
	{
		upsilon_times = var;
	}
	else if ( name == "gamma" )
	{
		gamma = var;
	}
	else if ( name == "gammaTimes" )
	{
		gamma_times = var;
	}
	else if ( name == "rho" )
	{
		rho = var;
	}
	else if ( name == "rhoTimes" )
	{
		rho_times = var;
	}
	else if ( name == "xi" )
	{
		xi = var;
	}
	else if ( name == "xiTimes" )
	{
		xi_times = var;
	}
	else if ( name == "eta" )
	{
		eta = var;
	}
	else if ( name == "etaTimes" )
	{
		eta_times = var;
	}
	else if ( name == "omega" )
	{
		omega = var;
	}
	else if ( name == "omegaTimes" )
	{
		omega_times = var;
	}
	else if ( name == "zeta" )
	{
		zeta = var;
	}
	else if ( name == "rootFreq" )
	{
		root_frequencies = var;
	}
	else if ( name == "condition" )
	{
		condition = var;
	}
	else
	{
		Distribution::setConstParameter(name, var);
	}

}
