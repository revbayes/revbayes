#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_GLHBDSP.h"

#include "RlDiscretizedContinuousCharacterData.h"
#include "GeneralizedLineageHeterogeneousBirthDeathSamplingProcess.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlCladogeneticProbabilityMatrix.h"
#include "RlRateGenerator.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
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

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_GLHBDSP::Dist_GLHBDSP() : TypedDistribution<TimeTree>()
{

}

Dist_GLHBDSP::~Dist_GLHBDSP()
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
    // the start age
    RevBayesCore::TypedDagNode<double>* sa = static_cast<const RealPos &>( start_age->getRevObject() ).getDagNode();

    // the start type
    bool use_origin = start_type == "originAge" ? true : false;

    // sampling condition
    const std::string& cond = static_cast<const RlString &>( condition->getRevObject() ).getValue();

    if ( cond != "time" )
    {
        if ( start_type == "originAge" )
        {
        	if ( cond != "survival" && cond != "sampled" && cond != "sampledExtant" && cond != "tree" && cond != "treeExtant" )
        	{
        		throw RbException() << "Cannot condition on " << cond << " when starting from the origin." ; 
        	}
        }
        else if ( start_type == "rootAge" )
        {
        	if ( cond != "survival" && cond != "sampled" && cond != "sampledExtant" && cond != "sampledMRCA" )
        	{
        		throw RbException() << "Cannot condition on " << cond << " when starting from the root." ; 
        	}
        }
        else
        {
        	throw RbException("How did you get here???");
        }
    }

    // taxa
    std::vector<RevBayesCore::Taxon> tax = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    // root frequency
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* root_freq = static_cast<const Simplex &>( pi->getRevObject() ).getDagNode();;

    // number of states
    size_t num_states = size_t(static_cast<const Natural &>( n_states->getRevObject() ).getValue());
    if ( num_states < 2 )
    {
    	throw RbException("nStates must be 2 or greater.");
    }

    // number of processors
    size_t num_procs = size_t(static_cast<const Natural &>( n_proc->getRevObject() ).getValue());

    // tolerances
    double absolute_tolerance = (double)static_cast<const RealPos &>(abs_tol->getRevObject()).getValue();
    double relative_tolerance = (double)static_cast<const RealPos &>(rel_tol->getRevObject()).getValue();

    // indexing
    bool zi = static_cast<const RlBoolean &>( zero_index->getRevObject() ).getValue();

    // make the distribution
    RevBayesCore::GeneralizedLineageHeterogeneousBirthDeathSamplingProcess* d = new RevBayesCore::GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(tax, sa, cond, root_freq, num_states, use_origin, zi, num_procs, absolute_tolerance, relative_tolerance);

    // speciation rates
    if ( lambda->getRevObject() != RevNullObject::getInstance() )
    {
    	if ( lambda->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    	{ // case 1: time homogeneous
    		RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* l_rates = static_cast<const ModelVector< RealPos > &>( lambda->getRevObject() ).getDagNode();
    	    d->setLambda(l_rates);
    	}
    	else if ( lambda->getRevObject().isType( ModelVector< ModelVector<RealPos> >::getClassTypeSpec() ) )
		{ // case 2: time heterogeneous
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* l_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( lambda->getRevObject() ).getDagNode();
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           l_times = static_cast<const ModelVector< RealPos > &>( lambda_times->getRevObject() ).getDagNode();
    	    d->setLambda(l_rates, l_times);
		}
    	else
    	{
    		throw RbException("How did you get here?");
    	}
    }
    else
    {
    	throw RbException("Must provide a vector of speciation rates.");
    }

    // extinction rates
    if ( mu->getRevObject() != RevNullObject::getInstance() )
    {
    	if ( mu->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    	{ // case 1: time homogeneous
    		RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* m_rates = static_cast<const ModelVector< RealPos > &>( mu->getRevObject() ).getDagNode();
    	    d->setMu(m_rates);
    	}
    	else if ( mu->getRevObject().isType( ModelVector< ModelVector<RealPos> >::getClassTypeSpec() ) )
		{ // case 2: time heterogeneous
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* m_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( mu->getRevObject() ).getDagNode();
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           m_times = static_cast<const ModelVector< RealPos > &>( mu_times->getRevObject() ).getDagNode();
    	    d->setMu(m_rates, m_times);
		}
    	else
    	{
    		throw RbException("How did you get here?");
    	}
    }

    // sampling rates
    if ( phi->getRevObject() != RevNullObject::getInstance() )
    {
    	if ( phi->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    	{ // case 1: time homogeneous
    		RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* p_rates = static_cast<const ModelVector< RealPos > &>( phi->getRevObject() ).getDagNode();
    	    d->setPhi(p_rates);
    	}
    	else if ( phi->getRevObject().isType( ModelVector< ModelVector<RealPos> >::getClassTypeSpec() ) )
		{ // case 2: time heterogeneous
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* p_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( phi->getRevObject() ).getDagNode();
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           p_times = static_cast<const ModelVector< RealPos > &>( phi_times->getRevObject() ).getDagNode();
    	    d->setPhi(p_rates, p_times);
		}
    	else
    	{
    		throw RbException("How did you get here?");
    	}
    }

    // destructive-sampling rates
    if ( delta->getRevObject() != RevNullObject::getInstance() )
    {
    	if ( delta->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    	{ // case 1: time homogeneous
    		RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* d_rates = static_cast<const ModelVector< RealPos > &>( delta->getRevObject() ).getDagNode();
    	    d->setDelta(d_rates);
    	}
    	else if ( delta->getRevObject().isType( ModelVector< ModelVector<RealPos> >::getClassTypeSpec() ) )
		{ // case 2: time heterogeneous
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* d_rates = static_cast<const ModelVector< ModelVector<RealPos> > &>( delta->getRevObject() ).getDagNode();
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           d_times = static_cast<const ModelVector< RealPos > &>( delta_times->getRevObject() ).getDagNode();
    	    d->setDelta(d_rates, d_times);
		}
    	else
    	{
    		throw RbException("How did you get here?");
    	}
    }

    // mass-speciation event
    if ( upsilon->getRevObject() != RevNullObject::getInstance() )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* u_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( upsilon->getRevObject() ).getDagNode();
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           u_times = static_cast<const ModelVector< RealPos> &>( upsilon_times->getRevObject() ).getDagNode();
        d->setUpsilon(u_probs, u_times);
    }

    // mass-extinction event
    if ( gamma->getRevObject() != RevNullObject::getInstance() )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* g_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( gamma->getRevObject() ).getDagNode();
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           g_times = static_cast<const ModelVector< RealPos> &>( gamma_times->getRevObject() ).getDagNode();
        d->setGamma(g_probs, g_times);
    }

    // mass-sampling event
    if ( rho->getRevObject() != RevNullObject::getInstance() )
    {
    	if ( rho->getRevObject().isType( Probability::getClassTypeSpec() ) )
    	{ // simple case
    		RevBayesCore::TypedDagNode< double >* r_probs = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();
    		d->setRho(r_probs);
    	}
    	else
    	{ // state-depenent case
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* r_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( rho->getRevObject() ).getDagNode();
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           r_times = static_cast<const ModelVector< RealPos> &>( rho_times->getRevObject() ).getDagNode();
            d->setRho(r_probs, r_times);
    	}
    }

    // mass-destructive-sampling event
    if ( xi->getRevObject() != RevNullObject::getInstance() )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RbVector<double> > >* x_probs = static_cast<const ModelVector< ModelVector<Probability> > &>( xi->getRevObject() ).getDagNode();
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                           x_times = static_cast<const ModelVector< RealPos> &>( xi_times->getRevObject() ).getDagNode();
        d->setXi(x_probs, x_times);
    }

    // anagenetic changes
    if ( eta->getRevObject() != RevNullObject::getInstance() )
    {
    	if ( eta->getRevObject().isType( RealPos::getClassTypeSpec() ) )
    	{
    		RevBayesCore::TypedDagNode< double >* h_rate = static_cast<const RealPos &>( eta->getRevObject() ).getDagNode();
    	    d->setEta(h_rate);
    	}
    	else if ( eta->getRevObject().isType( RateGenerator::getClassTypeSpec() ) )
    	{ // case 2: time homogeneous
    		RevBayesCore::TypedDagNode< RevBayesCore::RateGenerator >* h_mats = static_cast<const RateGenerator &>( eta->getRevObject() ).getDagNode();
    	    d->setEta(h_mats);
    	}
    	else if ( eta->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
		{ // case 3: time heterogeneous
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::RateGenerator > >* h_mats  = static_cast<const ModelVector<RateGenerator> &>( eta->getRevObject() ).getDagNode();
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                        h_times = static_cast<const ModelVector< RealPos > &>( eta_times->getRevObject() ).getDagNode();
    	    d->setEta(h_mats, h_times);
		}
    	else
    	{
    		throw RbException("How did you get here?");
    	}
    }

    // cladogenetic changes
    if ( omega->getRevObject() != RevNullObject::getInstance() )
    {
    	if ( omega->getRevObject().isType( CladogeneticProbabilityMatrix::getClassTypeSpec() ) )
    	{ // case 1: time homogeneous
    		RevBayesCore::TypedDagNode< RevBayesCore::CladogeneticProbabilityMatrix >* w_mats = static_cast<const CladogeneticProbabilityMatrix &>( omega->getRevObject() ).getDagNode();
    	    d->setOmega(w_mats);
    	}
    	else if ( omega->getRevObject().isType( ModelVector<CladogeneticProbabilityMatrix>::getClassTypeSpec() ) )
		{ // case 2: time heterogeneous
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::CladogeneticProbabilityMatrix > >* w_mats  = static_cast<const ModelVector<CladogeneticProbabilityMatrix> &>( omega->getRevObject() ).getDagNode();
    	    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*                                        w_times = static_cast<const ModelVector< RealPos > &>( omega_times->getRevObject() ).getDagNode();
    	    d->setOmega(w_mats, w_times);
		}
    	else
    	{
    		throw RbException("How did you get here?");
    	}
    }

    // mass-extinction mediated state-change event
    if ( zeta->getRevObject() != RevNullObject::getInstance() )
    {
    	RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::MatrixReal > >* z_mats = static_cast<const ModelVector< StochasticMatrix > &>( zeta->getRevObject() ).getDagNode();;
        d->setZeta(z_mats);
    }

    // refresh SSE process outcome given assigned parameter values
    d->redrawValue();
    
    // return the distribution
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
    a_names.push_back( "TimeHeterogeneousLineageHeterogeneousBirthDeathSamplingProcess" );

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

    ArgumentRules* clampCharDataArgRules = new ArgumentRules();
    std::vector<TypeSpec> data_types;
    data_types.push_back( AbstractHomologousDiscreteCharacterData::getClassTypeSpec() );
    data_types.push_back( DiscretizedContinuousCharacterData::getClassTypeSpec() );
    clampCharDataArgRules->push_back( new ArgumentRule( "value", data_types, "The observed value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "clampCharData", RlUtils::Void, clampCharDataArgRules ) );

    ArgumentRules* getCharDataArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "getCharData", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), getCharDataArgRules ) );

    ArgumentRules* dumpModelArgRules = new ArgumentRules();
    dumpModelArgRules->push_back( new ArgumentRule( "filename", RlString::getClassTypeSpec(), "The file to write the model in.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "dumpModel", RlUtils::Void, dumpModelArgRules ) );

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
    	std::vector<std::string> age_types;
        age_types.push_back("rootAge");
        age_types.push_back("originAge");
        dist_member_rules.push_back( new ArgumentRule( age_types, RealPos::getClassTypeSpec(), "The start time of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        // the root frequency
        dist_member_rules.push_back( new ArgumentRule( "pi",     Simplex::getClassTypeSpec(),                                     "Frequencies of each state at the beginning of the process.",                                 ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        // regular events (and times)
        std::vector<TypeSpec> rate_types;
        rate_types.push_back( ModelVector< RealPos>::getClassTypeSpec() );
        rate_types.push_back( ModelVector< ModelVector<RealPos> >::getClassTypeSpec() );

        dist_member_rules.push_back( new ArgumentRule( "lambda",       rate_types,                                                      "The vector of speciation rates for each time interval.",                                     ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "lambdaTimes",  ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which speciation rates change.",                                                ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );

        dist_member_rules.push_back( new ArgumentRule( "mu",           rate_types,                                                      "The vector of extinction rates for each time interval.",                                     ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "muTimes",      ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which extinction rates change.",                                                ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );

        dist_member_rules.push_back( new ArgumentRule( "phi",          rate_types,                                                      "The vector of sampling rates for each time interval.",                                       ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "phiTimes",     ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which sampling rates change.",                                                  ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );

        dist_member_rules.push_back( new ArgumentRule( "delta",        rate_types,                                                      "The vector of destructive-sampling rates for each time interval.",                           ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "deltaTimes",   ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which destructive-sampling rates change.",                                      ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );

        // mass events
        dist_member_rules.push_back( new ArgumentRule( "upsilon",      ModelVector< ModelVector<Probability> >::getClassTypeSpec(),     "The vector of speciation probabilities for each mass-speciation event.",                     ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "upsilonTimes", ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which mass-speciation events occur.",                                           ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "gamma",        ModelVector< ModelVector<Probability> >::getClassTypeSpec(),     "The vector of extinction probabilities for each mass-extinction event.",                     ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "gammaTimes",   ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which mass-extinction events occur.",                                           ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );


        std::vector<TypeSpec> rho_types;
        rho_types.push_back( Probability::getClassTypeSpec() );
        rho_types.push_back( ModelVector< ModelVector<Probability> >::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "rho",          rho_types,                                                       "The vector of sampling probabilities for each mass-sampling event.",                         ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "rhoTimes",     ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which mass-sampling events occur.",                                             ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        dist_member_rules.push_back( new ArgumentRule( "xi",           ModelVector< ModelVector<Probability> >::getClassTypeSpec(),     "The vector of destructive-sampling probabilities for each mass-destructive-sampling event.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "xiTimes",      ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which mass-destructive-sampling events occur.",                                 ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        // state changes
        std::vector<TypeSpec> eta_types;
        eta_types.push_back( RealPos::getClassTypeSpec() );
        eta_types.push_back( RateGenerator::getClassTypeSpec() );
        eta_types.push_back( ModelVector<RateGenerator>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "eta",          eta_types,                                                       "The anagenetic rates of change for each time interval.",                                     ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "etaTimes",     ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which the anagenetic rates change.",                                            ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );

        std::vector<TypeSpec> omega_types;
        omega_types.push_back( CladogeneticProbabilityMatrix::getClassTypeSpec() );
        omega_types.push_back( ModelVector<CladogeneticProbabilityMatrix>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "omega",        omega_types,                                                     "The cladogenetic event probabilities for each time interval.",                               ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "omegaTimes",   ModelVector< RealPos >::getClassTypeSpec(),                      "The times at which the cladogenetic rates change.",                                          ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );

        dist_member_rules.push_back( new ArgumentRule( "zeta",         ModelVector<StochasticMatrix>::getClassTypeSpec(),               "The probabilities of change for each mass-extinction event.",                                ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        // conditioning
		std::vector<std::string> options_condition;
		options_condition.push_back( "time" );
		options_condition.push_back( "survival" );
		options_condition.push_back( "sampled" );
		options_condition.push_back( "sampledExtant" );
		options_condition.push_back( "sampledMRCA" );
		options_condition.push_back( "tree" );
		options_condition.push_back( "treeExtant" );
		dist_member_rules.push_back( new OptionRule( "condition", new RlString("time"), options_condition, "The condition of the process." ) );

		// taxa
        dist_member_rules.push_back( new ArgumentRule( "taxa", ModelVector<Taxon>::getClassTypeSpec(), "The taxa in the tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

		// number of states
        dist_member_rules.push_back( new ArgumentRule( "nStates", Natural::getClassTypeSpec(), "The number of discrete states.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(2) ) );

		// number of processors
        dist_member_rules.push_back( new ArgumentRule( "nProc", Natural::getClassTypeSpec(), "The number of processors for parallel calculations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );

        // tolerances
        dist_member_rules.push_back( new ArgumentRule( "absTol", RealPos::getClassTypeSpec(), "The absolute tolerance of the numerical integrator.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1e-7) ) );
        dist_member_rules.push_back( new ArgumentRule( "relTol", RealPos::getClassTypeSpec(), "The relative tolerance of the numerical integrator.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1e-7) ) );

        // zero indexing
        dist_member_rules.push_back( new ArgumentRule( "zeroIndex", RlBoolean::getClassTypeSpec(), "Does the state space include zero?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );


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
		start_type = name;
		start_age  = var;
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
	else if ( name == "pi" )
	{
		pi = var;
	}
	else if ( name == "condition" )
	{
		condition = var;
	}
	else if ( name == "taxa" )
	{
		taxa = var;
	}
	else if ( name == "nStates" )
	{
		n_states = var;
	}
	else if ( name == "nProc" )
	{
		n_proc = var;
	}
	else if ( name == "absTol")
	{
		abs_tol = var;
	}
	else if ( name == "relTol")
	{
		rel_tol = var;
	}
	else if ( name == "zeroIndex" )
	{
		zero_index = var;
	}
	else
	{
		Distribution::setConstParameter(name, var);
	}

}
