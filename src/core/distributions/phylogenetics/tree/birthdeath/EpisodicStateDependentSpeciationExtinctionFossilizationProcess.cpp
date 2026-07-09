#include <boost/assign/list_of.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "SSE_ODE.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "DistributionExponential.h"
#include "HomologousDiscreteCharacterData.h"
#include "EpisodicStateDependentSpeciationExtinctionFossilizationProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateMatrix_JC.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RlString.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "AbstractDiscreteTaxonData.h"
#include "AbstractTaxonData.h"
#include "Cloneable.h"
#include "DiscreteCharacterState.h"
#include "DiscreteTaxonData.h"
#include "NaturalNumbersState.h"
#include "RateGenerator.h"
#include "RbBitSet.h"
#include "RbException.h"
#include "RbSettings.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "Simplex.h"
#include "StringUtilities.h"
#include "Taxon.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TreeDiscreteCharacterData.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "boost/numeric/odeint.hpp" // IWYU pragma: keep
#include "RlUserInterface.h"  // for RBOUT

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }


using namespace RevBayesCore;


/**
 * Constructor.
 *
 * The constructor connects the parameters of the birth-death process (DAG structure)
 * and initializes the probability density by computing the combinatorial constant of the tree structure.
 */
EpisodicStateDependentSpeciationExtinctionFossilizationProcess::EpisodicStateDependentSpeciationExtinctionFossilizationProcess(const TypedDagNode<double> *age,
                                                                                                                               const TypedDagNode< Simplex >* p,
                                                                                                                               const std::string &cdt,
                                                                                                                               bool uo,
                                                                                                                               size_t min_num_lineages,
                                                                                                                               size_t max_num_lineages,
                                                                                                                               size_t exact_num_lineages,
                                                                                                                               double max_t,
                                                                                                                               bool prune,
                                                                                                                               bool condition_on_tip_states,
                                                                                                                               bool condition_on_num_tips,
                                                                                                                               bool condition_on_tree,
                                                                                                                               std::int64_t pr) : TypedDistribution<Tree>( new TreeDiscreteCharacterData() ),
    condition( cdt ),
    active_likelihood( std::vector<bool>(5, 0) ),
    changed_nodes( std::vector<bool>(5, false) ),
    dirty_nodes( std::vector<bool>(5, true) ),
    node_partial_likelihoods( std::vector<std::vector<std::vector<double> > >(5, std::vector<std::vector<double> >(2,std::vector<double>(2*p->getValue().size(),0))) ),
    extinction_probabilities( std::vector<std::vector<double> >( 500.0, std::vector<double>( p->getValue().size(), 0) ) ),
    num_states( p->getValue().size() ),
    scaling_factors( std::vector<std::vector<double> >(5, std::vector<double>(2,0.0) ) ),
    use_cladogenetic_events( false ),
    use_origin( uo ),
    sample_character_history( false ),
    average_speciation( std::vector<double>(5, 0.0) ),
    average_extinction( std::vector<double>(5, 0.0) ),
    num_shift_events( std::vector<std::int64_t>(5, 0.0) ),
    time_in_states( std::vector<double>(p->getValue().size(), 0.0) ),
    simmap( "" ),
    cladogenesis_matrix( NULL ),
    process_age( age ),
    lambda_const(NULL),
    lambda_var(NULL),
    mu_const(NULL),
    mu_var(NULL),
    phi_const( NULL),
    phi_var( NULL),
    eta_const( NULL),
    eta_var( NULL),
    epoch_times_lambda( NULL ),
    epoch_times_mu( NULL ),
    epoch_times_phi( NULL ),
    epoch_times_gamma( NULL ),
    epoch_times_eta( NULL ),
    epoch_times_Q( NULL ),
    pi( p ),
    survival_probs( NULL ),
    rho( NULL ),
    rho_per_state( NULL ),
    Q_default( p->getValue().size() ),
    min_num_lineages( min_num_lineages ),
    max_num_lineages( max_num_lineages ),
    exact_num_lineages( exact_num_lineages ),
    max_time( max_t ),
    prune_extinct_lineages( prune ),
    use_episodic_model( false ),
    condition_on_tip_states( condition_on_tip_states ),
    condition_on_num_tips( condition_on_num_tips ),
    condition_on_tree( condition_on_tree ),
    age_check_precision( pr ),
    NUM_TIME_SLICES( 500.0 )
{
    addParameter( pi );
    addParameter( rho );
    addParameter( process_age );
    
    if ( min_num_lineages > max_num_lineages )
    {
        throw RbException("minNumLineages cannot be greater than maxNumLineages.");
    }
    
    // set the length of the time slices used by the ODE for numerical integration
    dt = process_age->getValue() / NUM_TIME_SLICES * 1.0;

    value->getTreeChangeEventHandler().addListener( this );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
EpisodicStateDependentSpeciationExtinctionFossilizationProcess* EpisodicStateDependentSpeciationExtinctionFossilizationProcess::clone( void ) const
{
    EpisodicStateDependentSpeciationExtinctionFossilizationProcess* tmp = new EpisodicStateDependentSpeciationExtinctionFossilizationProcess( *this );
    tmp->getValue().getTreeChangeEventHandler().addListener(tmp);
    return tmp;
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
EpisodicStateDependentSpeciationExtinctionFossilizationProcess::~EpisodicStateDependentSpeciationExtinctionFossilizationProcess( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!

    // remove myself from the tree listeners
    value->getTreeChangeEventHandler().removeListener( this );
}


/**
 * Adds parameter-specific timeline to the set
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::addTimesToGlobalTimeline(std::set<double> &event_times, const RbVector<double>& par_times) const
{
  
    for (size_t i = 0; i < par_times.size(); ++i)
    {
        event_times.insert( par_times[i] );
    }

}


std::vector<double> EpisodicStateDependentSpeciationExtinctionFossilizationProcess::calculateTotalSpeciationRatePerState( double a ) const
{
    std::vector<double> total_rates = std::vector<double>(num_states, 0);
    if ( use_cladogenetic_events == true )
    {
        std::map<std::vector<unsigned>, double> eventMap;
        std::map<std::vector<unsigned>, double>::iterator it;
        
        // get cladogenesis event map (sparse speciation rate matrix)
        eventMap = cladogenesis_matrix->getValue().getEventMap();
        // iterate over each cladogenetic event possible
        for (it = eventMap.begin(); it != eventMap.end(); it++)
        {
            const std::vector<unsigned>& states = it->first;
            total_rates[states[0]] += it->second;
        }
    }
    else
    {
        total_rates = computeSpeciationRateAtTime( a );

    }
    return total_rates;
}


std::vector<double> EpisodicStateDependentSpeciationExtinctionFossilizationProcess::calculateTotalAnageneticRatePerState( double age ) const
{
    std::vector<double> total_rates = std::vector<double>(num_states, 0);
    const RateGenerator& rate_matrix = getEventRateMatrix( age );
    for (size_t i = 0; i < num_states; i++)
    {
        for (size_t j = 0; j < num_states; j++)
        {
            if (i != j)
            {
                total_rates[i] += rate_matrix.getRate(i, j, 0.0, getEventRate( age ));
            }
        }
    }

    return total_rates;
}


std::vector<double> EpisodicStateDependentSpeciationExtinctionFossilizationProcess::calculateExtinctionRatePerState( double a ) const
{
    
    return computeExtinctionRateAtTime( a );
}


const RbVector<double>& EpisodicStateDependentSpeciationExtinctionFossilizationProcess::computeExtinctionRateAtTime( double a ) const
{
    
    if ( use_episodic_model )
    {
        // get the rates for this time from the variable/episodic rates
        size_t index = computeEpochIndex(a);
        
        if ( index >= mu.size() ) throw RbException("Didn't rescale vector mu correctly.");
        
        return mu[index];
    }
    else
    {
        // use the constant rates
        const RbVector<double> &ext_rates = mu_const->getValue();
        return ext_rates;
    }
}


const RbVector<double>& EpisodicStateDependentSpeciationExtinctionFossilizationProcess::computeFossilizationRateAtTime( double a ) const
{
    
    if ( use_episodic_model )
    {
        // get the rates for this time from the variable/episodic rates
        size_t index = computeEpochIndex(a);
        
        if ( index >= phi.size() ) throw RbException("Didn't rescale vector phi correctly.");

        
        const RbVector<double> &fos_rates = phi[index];

        return fos_rates;
    }
    else
    {
        // use the constant rates
        const RbVector<double> &fos_rates = phi_const->getValue();
        return fos_rates;
    }
}


const RbVector<double>& EpisodicStateDependentSpeciationExtinctionFossilizationProcess::computeSpeciationRateAtTime( double a ) const
{
    
    if ( use_episodic_model )
    {
        // get the rates for this time from the variable/episodic rates
        size_t index = computeEpochIndex(a);

        if ( index >= lambda.size() ) throw RbException("Didn't rescale vector lambda correctly.");

        
        return lambda[index];
    }
    else
    {
        // use the constant rates
        const RbVector<double> &spe_rates = lambda_const->getValue();
        return spe_rates;
    }
}


const RbVector<double>& EpisodicStateDependentSpeciationExtinctionFossilizationProcess::computeSurvivalProbabilitiesAtTime( double a ) const
{
    
    // get the rates for this time from the variable/episodic rates
    size_t index = computeEpochIndex(a);
    
    if ( index >= gamma.size() ) throw RbException("Didn't rescale vector gamma correctly.");

        
//    const RbVector<double> &sp = survival_probs->getValue()[index];
    const RbVector<double> &sp = gamma[index];

    return sp;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double EpisodicStateDependentSpeciationExtinctionFossilizationProcess::computeLnProbability( void )
{
    
    // prepare the timelines and parameter vectors
    prepareTimeline();
    
    // check that the ages are in correct chronological order
    // i.e., no child is older than its parent
    const std::vector<TopologyNode*>& nodes = value->getNodes();
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        
        const TopologyNode &the_node = *(*it);
        if ( the_node.isRoot() == false )
        {
            
            if ( (the_node.getAge() - (*it)->getParent().getAge()) > 0 && the_node.isSampledAncestorTip() == false )
            {
                return RbConstants::Double::neginf;
            }
            else if ( (the_node.getAge() - (*it)->getParent().getAge()) > 0 && the_node.isSampledAncestorTip() == true )
            {
                return RbConstants::Double::neginf;
            }
            
        }
        
    }
    
    // check that the sampled ancestor nodes have a zero branch length
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        
        const TopologyNode &the_node = *(*it);
        if ( the_node.isSampledAncestorTip() == true )
        {
            
            if ( the_node.isFossil() == false )
            {
                return RbConstants::Double::neginf;
            }
            else if ( the_node.getBranchLength() > 0 )
            {
                return RbConstants::Double::neginf;
            }
            
        }
        
    }
    
    double num_initial_lineages = 2; // this needs to be a double!
    const TopologyNode& root = value->getRoot();

    if ( use_origin == true )
    {
        // If we are conditioning on survival from the origin,
        // then we must divide by 2 the log survival probability computed by AbstractBirthDeathProcess
        num_initial_lineages = 1;
    }
    // if conditioning on root, root node must be a "true" bifurcation event
    else if (root.getChild(0).isSampledAncestorTip() || root.getChild(1).isSampledAncestorTip())
    {
        return RbConstants::Double::neginf;
    }

    // present time
    double ra = root.getAge();
    double process_time = getOriginAge();
    
    if ( ra > process_time || ra != getRootAge() )
    {
        return RbConstants::Double::neginf;
    }
    
    const std::vector<TopologyNode*> &c = root.getChildren();

    for (std::vector<TopologyNode*>::const_iterator it = c.begin(); it != c.end(); ++it)
    {
        if ( ra < (*it)->getAge() )
        {
            return RbConstants::Double::neginf;
        }
    }

    if ( value->getNumberOfNodes() != dirty_nodes.size() )
    {
        resizeVectors(value->getNumberOfNodes());
    }
    
    // variable declarations and initialization
    double lnProbTimes = 0;
    
    // conditioning on survival
    if ( condition == "survival" )
    {
        lnProbTimes = - log( pSurvival(0, process_time,num_initial_lineages>1) );
    }
    
    // multiply the probability of a descendant of the initial species
    lnProbTimes += computeRootLikelihood();
    
    return lnProbTimes + lnProbTreeShape();
}


size_t EpisodicStateDependentSpeciationExtinctionFossilizationProcess::computeEpochIndex(double a) const
{
    
    size_t index = 0;
    
    if ( use_episodic_model )
    {
        
        while ( index < global_timeline.size() && a > global_timeline[index] )
        {
            ++index;
        }
    }
    
    return index;
}


double EpisodicStateDependentSpeciationExtinctionFossilizationProcess::computeEpochEnd(size_t i) const
{
    if ( use_episodic_model == false )
    {
        return RbConstants::Double::inf;
    }

    if ( i >= global_timeline.size() )
    {
        return RbConstants::Double::inf;
    }
    else
    {
        return global_timeline[i];
    }
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::computeNodeProbability(const RevBayesCore::TopologyNode &node, size_t node_index) const
{
    std::vector<double> &node_likelihood  = node_partial_likelihoods[node_index][active_likelihood[node_index]];

    // check for recomputation
    if ( dirty_nodes[node_index] == true || sample_character_history == true )
    {
        // mark as computed
        dirty_nodes[node_index] = false;
        
        if ( node.isTip() == true )
        {
            // this is a tip node
            TreeDiscreteCharacterData* tree = static_cast<TreeDiscreteCharacterData*>( this->value );

            std::vector<double> sampling;
            std::vector<double> extinction;
            if ( rho != NULL && rho_per_state == NULL )
            {
                sampling   = std::vector<double>(num_states, rho->getValue());
                extinction = std::vector<double>(num_states, 1.0 - rho->getValue());
            }
            else if ( rho == NULL && rho_per_state != NULL )
            {
                sampling   = rho_per_state->getValue();
                extinction = std::vector<double>(num_states, 1.0);
                for (size_t i=0; i<num_states; ++i)
                {
                    extinction[i] = 1.0 - sampling[i];
                }
            }
            else
            {
                throw RbException("Either a global sampling fraction or state-specific sampling fraction needs to be set.");
            }
            
            
            if ( node.isFossil() && node.getAge() > 1E-4 )
            {
                if ( phi_const == NULL && phi_var == NULL )
                {
                    throw(RbException("Tree has serially sampled tips, but no serial sampling rate was provided."));
                }
                
                double node_age= node.getAge();
                sampling = computeFossilizationRateAtTime( node_age );
                extinction = pExtinction(0.0, node_age);
            }
            
            RbBitSet obs_state(num_states);
            obs_state.set();
            bool gap = true;

            if ( tree->hasCharacterData() == true )
            {
                const DiscreteCharacterState &state = tree->getCharacterData().getTaxonData( node.getTaxon().getName() )[0];
                obs_state = state.getState();
                gap = (state.isMissingState() == true || state.isGapState() == true);
            }

            if (obs_state.size() > num_states)
            {
                throw RbException()<<"SSE model has "<<num_states<<" states, but observed data set has "<<obs_state.size()<<" states!";
            }
            else if (obs_state.size() < num_states)
            {
                std::ostringstream o;
                o<<"Warning: SSE model has "<<num_states<<" states, but observed data set has only "<<obs_state.size()<<" states!";
                RBOUT(o.str());
            }

            double all_states_impossible = true;
            for (size_t j = 0; j < num_states; ++j)
            {

                node_likelihood[j] = extinction[j];

                if ( obs_state.test( j ) == true || gap == true )
                {
                    if ( node.isFossil() && node.getAge() > 1E-4 )
                    {
                        node_likelihood[num_states+j] = sampling[j] * extinction[j];
                    }
                    else
                    {
                        node_likelihood[num_states+j] = sampling[j];
                    }
                }
                else
                {
                    node_likelihood[num_states+j] = 0.0;
                }

                if (node_likelihood[num_states+j] > 0) all_states_impossible = false;
            }

            // Should we print something here?  Possibly this never happens.
            assert(not all_states_impossible);
        }
        else
        {
            
            // this is an internal node
            const TopologyNode          &left           = node.getChild(0);
            size_t                      left_index      = left.getIndex();
            computeNodeProbability( left, left_index );
            const TopologyNode          &right          = node.getChild(1);
            size_t                      right_index     = right.getIndex();
            computeNodeProbability( right, right_index );
            
            // get the likelihoods of descendant nodes
            const std::vector<double> &left_likelihoods  = node_partial_likelihoods[left_index][active_likelihood[left_index]];
            const std::vector<double> &right_likelihoods = node_partial_likelihoods[right_index][active_likelihood[right_index]];

            std::map<std::vector<unsigned>, double> eventMap;
            std::vector<double> speciation_rates;
            if ( use_cladogenetic_events == true )
            {
                // get cladogenesis event map (sparse speciation rate matrix)
                eventMap = cladogenesis_matrix->getValue().getEventMap();
            }
            else
            {
                speciation_rates = computeSpeciationRateAtTime( node.getAge() );
            }
            
            bool speciation_node = true;
            if ( left.isSampledAncestorTip() || right.isSampledAncestorTip() )
            {
                speciation_node = (phi_const == NULL && phi_var == NULL);
            }

            // merge descendant likelihoods
            for (size_t i=0; i<num_states; ++i)
            {
                node_likelihood[i] = left_likelihoods[i];

                if ( use_cladogenetic_events == true && speciation_node == true )
                {
                    
                    double like_sum = 0.0;
                    std::map<std::vector<unsigned>, double>::iterator it;
                    for (it = eventMap.begin(); it != eventMap.end(); it++)
                    {
                        const std::vector<unsigned>& states = it->first;
                        double speciation_rate = it->second;
                        if (i == states[0])
                        {
                            double likelihoods = left_likelihoods[num_states + states[1]] * right_likelihoods[num_states + states[2]];
                            like_sum += speciation_rate * likelihoods;
                        }
                    }
                    node_likelihood[num_states + i] = like_sum;
                    
                }
                else
                {
                    node_likelihood[num_states + i] = left_likelihoods[num_states + i] * right_likelihoods[num_states + i];
                    node_likelihood[num_states + i] *= speciation_node ? speciation_rates[i] : 1.0;
                }
            }
            
        }
        
        double begin_age = node.getAge();
        double end_age = node.getParent().getAge();
        
        if ( node.isSampledAncestorTip() == false )
        {
            // calculate likelihoods for this branch
            if ( sample_character_history == false )
            {
                // numerically integrate over the entire branch length
                numericallyIntegrateProcess(node_likelihood, begin_age, end_age, true, false);
            }
            else
            {
                // calculate the conditional likelihoods for each time slice moving
                // along this branch backwards in time from the tip towards the root

                std::vector<std::vector<double> > branch_likelihoods;
                size_t current_dt = 0;
                
                // calculate partial likelihoods for each time slice and store them in branch_likelihoods
                while ( (current_dt * dt) + begin_age < end_age )
                {

                    std::vector<double> dt_likelihood;

                    double current_dt_start = (current_dt * dt) + begin_age;
                    double current_dt_end = ((current_dt + 1) * dt) + begin_age;
                    if (current_dt_end > end_age)
                    {
                        current_dt_end = end_age;
                    }
                    numericallyIntegrateProcess(node_likelihood, current_dt_start, current_dt_end, true, false);

                    std::vector<double>::const_iterator first = node_likelihood.begin() + num_states;
                    std::vector<double>::const_iterator last = node_likelihood.begin() + (num_states * 2);
                    dt_likelihood = std::vector<double>(first, last);

                    branch_likelihoods.push_back(dt_likelihood);
                    current_dt++;

                }
                
                // save the branch conditional likelihoods
                branch_partial_likelihoods[node_index] = branch_likelihoods;
            }
        }
        
        if ( RbSettings::userSettings().getUseScaling() == true ) //&& node_index % RbSettings::userSettings().getScalingDensity() == 0 )
        {
            // rescale the conditional likelihoods at the "end" of the branch
            double max = 0.0;
            for (size_t i=0; i<num_states; ++i)
            {
                if ( node_likelihood[num_states+i] > max )
                {
                    max = node_likelihood[num_states+i];
                }
            }
//            max *= num_states;

            if (max > 0)
            {
                assert(std::isfinite(max) and std::isfinite(1/max));
                for (size_t i=0; i<num_states; ++i)
                {
                    node_likelihood[num_states+i] /= max;
                }

                scaling_factors[node_index][active_likelihood[node_index]] = log(max);

                if ( node.isTip() == false )
                {
                    const TopologyNode          &left           = node.getChild(0);
                    size_t                      left_index      = left.getIndex();
                    const TopologyNode          &right          = node.getChild(1);
                    size_t                      right_index     = right.getIndex();
                    scaling_factors[node_index][active_likelihood[node_index]] += scaling_factors[left_index][active_likelihood[left_index]] + scaling_factors[right_index][active_likelihood[right_index]];
                }
            }

        }
        
    }

}


double EpisodicStateDependentSpeciationExtinctionFossilizationProcess::computeRootLikelihood( void ) const
{
    // get the likelihoods of descendant nodes
    const TopologyNode     &root            = value->getRoot();
    size_t                  node_index      = root.getIndex();
    const TopologyNode     &left            = root.getChild(0);
    size_t                  left_index      = left.getIndex();
    computeNodeProbability( left, left_index );
    const TopologyNode     &right           = root.getChild(1);
    size_t                  right_index     = right.getIndex();
    computeNodeProbability( right, right_index );

    // get the likelihoods of descendant nodes
    const std::vector<double> &left_likelihoods  = node_partial_likelihoods[left_index][active_likelihood[left_index]];
    const std::vector<double> &right_likelihoods = node_partial_likelihoods[right_index][active_likelihood[right_index]];

    std::vector<double> &node_likelihood  = node_partial_likelihoods[node_index][active_likelihood[node_index]];

    std::map<std::vector<unsigned>, double> eventMap;
    std::vector<double> speciation_rates;
    if ( use_cladogenetic_events == true )
    {
        // get cladogenesis event map (sparse speciation rate matrix)
        eventMap = cladogenesis_matrix->getValue().getEventMap();
    }
    else
    {
        speciation_rates = computeSpeciationRateAtTime( root.getAge() );
    }

    bool speciation_node = true;
    if ( left.isSampledAncestorTip() || right.isSampledAncestorTip() )
    {
        speciation_node = (phi_const == NULL && phi_var == NULL);
    }

    // merge descendant likelihoods
    for (size_t i=0; i<num_states; ++i)
    {
        node_likelihood[i] = left_likelihoods[i];

        // MRM 03/23/2020: the root is a cladogenetic event, so the
        // cladogenetic events should be included at the root. right now,
        // I am just forcing the cladogenetic events at the root, but there
        // may be a better solution.
        if ( use_cladogenetic_events == true && speciation_node == true )
        {

            double like_sum = 0.0;
            std::map<std::vector<unsigned>, double>::iterator it;
            for (it = eventMap.begin(); it != eventMap.end(); it++)
            {
                const std::vector<unsigned>& states = it->first;
                double speciation_rate = it->second;
                if (i == states[0])
                {
                    double likelihoods = left_likelihoods[num_states + states[1]] * right_likelihoods[num_states + states[2]];
                    like_sum += speciation_rate * likelihoods;
                }
            }
            node_likelihood[num_states + i] = like_sum;

        }
        else
        {
            node_likelihood[num_states + i] = left_likelihoods[num_states + i] * right_likelihoods[num_states + i];
//            node_likelihood[num_states + i] *= (speciation_node ? speciation_rates[i] : 1.0);
        }
    }
    
    // calculate likelihoods for the root branch
    if ( use_origin == true )
    {
        double begin_age = getRootAge();
        double end_age = getOriginAge();

        if ( sample_character_history == false )
        {
            // numerically integrate over the entire branch length
            numericallyIntegrateProcess(node_likelihood, begin_age, end_age, true, false);
        }
        else
        {
            // calculate the conditional likelihoods for each time slice moving
            // along this branch backwards in time from the tip towards the root

            std::vector<std::vector<double> > branch_likelihoods;
            size_t current_dt = 0;

            // calculate partial likelihoods for each time slice and store them in branch_likelihoods
            while ( (current_dt * dt) + begin_age < end_age )
            {

                std::vector<double> dt_likelihood;

                double current_dt_start = (current_dt * dt) + begin_age;
                double current_dt_end = ((current_dt + 1) * dt) + begin_age;
                if (current_dt_end > end_age)
                {
                    current_dt_end = end_age;
                }
                numericallyIntegrateProcess(node_likelihood, current_dt_start, current_dt_end, true, false);

                std::vector<double>::const_iterator first = node_likelihood.begin() + num_states;
                std::vector<double>::const_iterator last = node_likelihood.begin() + (num_states * 2);
                dt_likelihood = std::vector<double>(first, last);

                branch_likelihoods.push_back(dt_likelihood);
                current_dt++;

            }

            // save the branch conditional likelihoods
            branch_partial_likelihoods[node_index] = branch_likelihoods;
        }
    }

    // sum the root likelihoods
    const RbVector<double> &freqs = getRootFrequencies();
    double prob = 0.0;

    for (size_t i = 0; i < num_states; ++i)
    {
        prob += freqs[i] * node_likelihood[num_states + i];
    }

    scaling_factors[node_index][active_likelihood[node_index]] = scaling_factors[left_index][active_likelihood[left_index]] + scaling_factors[right_index][active_likelihood[right_index]];
    
    return log(prob) + scaling_factors[node_index][active_likelihood[node_index]];
}


/**
 * Takes a par.size() < global_timeline.size() vector and makes it the correct size to work with our global timeline.
 * The parameter has its own reference timeline, which we use to find the rate in the global intervals.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::expandNonGlobalProbabilityParameterVector(std::vector<RbVector<double> > &par, const std::vector<double> &par_times, const RbVector<double>& default_prob) const
{
    // @TODO @efficiency: this works but it would be faster to auto-advance indices rather than have an internal loop
    // Store the original values so we can overwrite the vector
    std::vector<RbVector<double> > old_par = par;
    par.resize( global_timeline.size() );

    // For each time in the global timeline, find the rate according to this variable's own timeline
    for (size_t i=0; i<global_timeline.size(); ++i)
    {
        bool global_time_is_variable_time = false;
        for (size_t j=0; j<par_times.size(); ++j)
        {
            if ( fabs(par_times[j] - global_timeline[i]) < DBL_EPSILON )
            {
                // time is in variable's timeline
                par[i] = old_par[j];
                global_time_is_variable_time = true;
                break;
            }
        }

        // Time is not in variable's own timeline, probability of event here is 0
        if ( !global_time_is_variable_time )
        {
            par[i] = default_prob;
        }
    }

}

/**
 * Takes a par.size() < global_timeline.size() vector and makes it the correct size to work with our global timeline.
 * The parameter has its own reference timeline, which we use to find the rate in the global intervals.
 * This works only for parameters (lambda,mu,phi,r), where the global timeline is simply a finer grid than the variable-specific timelines.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::expandNonGlobalRateParameterVector(std::vector<RbVector<double> > &par, const std::vector<double> &par_times) const
{
    // Store the original values so we can overwrite the vector
    std::vector<RbVector<double> > old_par = par;

    // For each time in the global timeline, find the rate according to this variable's own timeline
    par.clear();
    for (size_t i=0; i<global_timeline.size(); ++i)
    {
      // Where is this global time interval in the variable's timeline?
      size_t idx = findIndex(global_timeline[i],par_times);
      par.push_back( old_par[idx] );
    }
    par.push_back( old_par[old_par.size()-1] );

}



/**
 * Takes a par.size() < global_timeline.size() vector and makes it the correct size to work with our global timeline.
 * The parameter has its own reference timeline, which we use to find the rate in the global intervals.
 * This works only for parameters (lambda,mu,phi,r), where the global timeline is simply a finer grid than the variable-specific timelines.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::expandNonGlobalRateParameterVector(std::vector<double> &par, const std::vector<double> &par_times) const
{
    // Store the original values so we can overwrite the vector
    std::vector<double> old_par = par;

    // For each time in the global timeline, find the rate according to this variable's own timeline
    par.clear();
    for (size_t i=0; i<global_timeline.size(); ++i)
    {
      // Where is this global time interval in the variable's timeline?
      size_t idx = findIndex(global_timeline[i],par_times);
      par.push_back( old_par[idx] );
    }
    par.push_back( old_par[old_par.size()-1] );

}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::fireTreeChangeEvent( const RevBayesCore::TopologyNode &n, const unsigned& m )
{
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );

}


const RevBayesCore::AbstractHomologousDiscreteCharacterData& EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getCharacterData() const
{
    return static_cast<TreeDiscreteCharacterData*>(this->value)->getCharacterData();
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::drawJointConditionalAncestralStates(std::vector<size_t>& startStates, std::vector<size_t>& endStates)
{
    // now begin the root-to-tip pass, drawing ancestral states conditional on the start states
    
    // get the likelihoods of descendant nodes
    const TopologyNode          &root               = value->getRoot();
    size_t                       node_index         = root.getIndex();
    const TopologyNode          &left               = root.getChild(0);
    size_t                       left_index         = left.getIndex();
    const std::vector< double > &left_likelihoods   = node_partial_likelihoods[left_index][active_likelihood[left_index]];
    const TopologyNode          &right              = root.getChild(1);
    size_t                       right_index        = right.getIndex();
    const std::vector< double > &right_likelihoods  = node_partial_likelihoods[right_index][active_likelihood[right_index]];
    
    
    std::map<std::vector<unsigned>, double> eventMap;
    std::vector<double> speciation_rates;
    if ( use_cladogenetic_events == true )
    {
        // get cladogenesis event map (sparse speciation rate matrix)
        eventMap = cladogenesis_matrix->getValue().getEventMap();
    }
    else
    {
        speciation_rates = computeSpeciationRateAtTime( root.getAge() );
    }
    
    // get root frequencies
    const RbVector<double> &freqs = getRootFrequencies();
    
    std::map<std::vector<unsigned>, double> sample_probs;
    double sample_probs_sum = 0.0;
    std::map<std::vector<unsigned>, double>::iterator it;
    
    // calculate probabilities for each state
    if ( use_cladogenetic_events == true )
    {
        // iterate over each cladogenetic event possible
        // and initialize probabilities for each clado event
        for (it = eventMap.begin(); it != eventMap.end(); it++)
        {
            const std::vector<unsigned>& states = it->first;
            double speciation_rate = it->second;
            
            // we need to sample from the ancestor, left, and right states jointly,
            // so keep track of the probability of each clado event
            double prob = left_likelihoods[num_states + states[1]] * right_likelihoods[num_states + states[2]];
            prob *= freqs[states[0]] * speciation_rate;
            sample_probs[ states ] = prob;
            sample_probs_sum += prob;
        }
    }
    else
    {
        for (size_t i = 0; i < num_states; i++)
        {
            double likelihood = left_likelihoods[num_states + i] * right_likelihoods[num_states + i] * speciation_rates[i];
            std::vector<unsigned> states = boost::assign::list_of(i)(i)(i);
            sample_probs[ states ] = likelihood * freqs[i];
            sample_probs_sum += likelihood * freqs[i];
        }
    }
    
    // sample ancestor, left, and right character states from probs
    size_t a = 0, l = 0, r = 0;
    
    if (sample_probs_sum == 0)
    {
        RandomNumberGenerator* rng = GLOBAL_RNG;
        size_t u = rng->uniform01() * sample_probs.size();
        size_t v = 0;
        for (it = sample_probs.begin(); it != sample_probs.end(); it++)
        {
            if (u < v)
            {
                const std::vector<unsigned>& states = it->first;
                a = states[0];
                l = states[1];
                r = states[2];
                endStates[node_index] = a;
                startStates[left_index] = l;
                startStates[right_index] = r;
                break;
             }
             v++;
         }
    }
    else
    {
        RandomNumberGenerator* rng = GLOBAL_RNG;
        double u = rng->uniform01() * sample_probs_sum;
       
        for (it = sample_probs.begin(); it != sample_probs.end(); it++)
        {
            u -= it->second;
            if (u < 0.0)
            {
                const std::vector<unsigned>& states = it->first;
                a = states[0];
                l = states[1];
                r = states[2];
                endStates[node_index] = a;
                startStates[node_index] = a;
                startStates[left_index] = l;
                startStates[right_index] = r;
                break;
            }
        }
    }
    
    // recurse towards tips
    recursivelyDrawJointConditionalAncestralStates(left, startStates, endStates);
    recursivelyDrawJointConditionalAncestralStates(right, startStates, endStates);
    
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::recursivelyDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<size_t>& startStates, std::vector<size_t>& endStates)
{
    
    size_t node_index = node.getIndex();
    
    if ( node.isTip() == true )
    {
        const AbstractHomologousDiscreteCharacterData& data = static_cast<TreeDiscreteCharacterData*>(this->value)->getCharacterData();
        const AbstractDiscreteTaxonData& taxon_data = data.getTaxonData( node.getName() );
        
        const DiscreteCharacterState &char_state = taxon_data.getCharacter(0);
        
        // get the observed state at the tip if it is known, otherwise simulate it
        if ( char_state.isAmbiguous() == false && char_state.isMissingState() == false )
        {
            endStates[node_index] = char_state.getStateIndex();
        }
        else
        {
            // initialize the conditional likelihoods for this branch
            std::vector< double > branch_conditional_probs = std::vector<double>(2 * num_states, 0);
            size_t start_state = startStates[node_index];
            branch_conditional_probs[ num_states + start_state ] = 1.0;
            
            // first calculate extinction likelihoods via a backward time pass
            double end_age = node.getParent().getAge();
            numericallyIntegrateProcess(branch_conditional_probs, 0, end_age, true, true);
            
            // now calculate conditional likelihoods along branch in forward time
            end_age        = node.getParent().getAge() - node.getAge();
            numericallyIntegrateProcess(branch_conditional_probs, 0, end_age, false, false);
            
            double total_prob = 0.0;
            for (size_t i = 0; i < num_states; ++i)
            {
                if ( char_state.isMissingState() == true || char_state.isGapState() == true || char_state.isStateSet(i) == true )
                {
                    total_prob += branch_conditional_probs[ num_states + i ];
                }
            }
            
            RandomNumberGenerator* rng = GLOBAL_RNG;
            double u = rng->uniform01() * total_prob;
            
            for (size_t i = 0; i < num_states; ++i)
            {
                
                if ( char_state.isMissingState() == true || char_state.isGapState() == true || char_state.isStateSet(i) == true )
                {
                    u -= branch_conditional_probs[ num_states + i ];
                    if ( u <= 0.0 )
                    {
                        endStates[node_index] = i;
                        break;
                    }
                    
                }
                
            }
            
        }
    }
    else
    {
        // sample characters by their probability conditioned on the branch's start state going to end states
        
        // initialize the conditional likelihoods for this branch
        std::vector< double > branch_conditional_probs = std::vector<double>(2 * num_states, 0);
        size_t start_state = startStates[node_index];
        branch_conditional_probs[ num_states + start_state ] = 1.0;

        // first calculate extinction likelihoods via a backward time pass
        double end_age = node.getParent().getAge();
        numericallyIntegrateProcess(branch_conditional_probs, 0, end_age, true, true);
        
        // now calculate conditional likelihoods along branch in forward time
        end_age        = node.getParent().getAge() - node.getAge();
        numericallyIntegrateProcess(branch_conditional_probs, 0, end_age, false, false);
        
        std::map<std::vector<unsigned>, double> event_map;
        std::vector<double> speciation_rates;
        if ( use_cladogenetic_events == true )
        {
            // get cladogenesis event map (sparse speciation rate matrix)
            event_map = cladogenesis_matrix->getValue().getEventMap();
        }
        else
        {
            speciation_rates = computeSpeciationRateAtTime( node.getAge() );
        }
        
        // get likelihoods of descendant nodes
        const TopologyNode &left = node.getChild(0);
        size_t left_index = left.getIndex();
        std::vector< double > left_likelihoods = node_partial_likelihoods[left_index][active_likelihood[left_index]];
        const TopologyNode &right = node.getChild(1);
        size_t right_index = right.getIndex();
        std::vector< double > right_likelihoods = node_partial_likelihoods[right_index][active_likelihood[right_index]];
        
        std::map<std::vector<unsigned>, double> sample_probs;
        double sample_probs_sum = 0.0;
        std::map<std::vector<unsigned>, double>::iterator it;

        // calculate probabilities for each state
        if ( use_cladogenetic_events == true )
        {
            // iterate over each cladogenetic event possible
            // and initialize probabilities for each clado event
            for (it = event_map.begin(); it != event_map.end(); it++)
            {
                const std::vector<unsigned>& states = it->first;
                double speciation_rate = it->second;
                
                // we need to sample from the ancestor, left, and right states jointly,
                // so keep track of the probability of each clado event
                double prob = left_likelihoods[num_states + states[1]] * right_likelihoods[num_states + states[2]];
                prob *= speciation_rate * branch_conditional_probs[num_states + states[0]];
                sample_probs[ states ] = prob;
                sample_probs_sum += prob;
            }
        }
        else
        {
            for (size_t i = 0; i < num_states; i++)
            {
                double prob = left_likelihoods[num_states + i] * right_likelihoods[num_states + i] * speciation_rates[i];
                prob *= branch_conditional_probs[num_states + i];
                std::vector<unsigned> states = boost::assign::list_of(i)(i)(i);
                sample_probs[ states ] = prob;
                sample_probs_sum += prob;
            }
        }
        
        // finally, sample ancestor, left, and right character states from probs
        size_t a = 0, l = 0, r = 0;

        if (sample_probs_sum == 0)
        {
            RandomNumberGenerator* rng = GLOBAL_RNG;
            size_t u = rng->uniform01() * sample_probs.size();
            size_t v = 0;
            for (it = sample_probs.begin(); it != sample_probs.end(); it++)
            {
                if (u < v)
                {
                    const std::vector<unsigned>& states = it->first;
                    a = states[0];
                    l = states[1];
                    r = states[2];
                    endStates[node_index] = a;
                    startStates[left_index] = l;
                    startStates[right_index] = r;
                    break;
                 }
                 v++;
             }
        }
        else
        {
            RandomNumberGenerator* rng = GLOBAL_RNG;
            double u = rng->uniform01() * sample_probs_sum;
            
            for (it = sample_probs.begin(); it != sample_probs.end(); it++)
            {
                u -= it->second;
                if (u < 0.0)
                {
                    const std::vector<unsigned>& states = it->first;
                    a = states[0];
                    l = states[1];
                    r = states[2];
                    endStates[node_index] = a;
                    startStates[left_index] = l;
                    startStates[right_index] = r;
                    break;
                }
            }
        }
        
        // recurse towards tips
        recursivelyDrawJointConditionalAncestralStates(left, startStates, endStates);
        recursivelyDrawJointConditionalAncestralStates(right, startStates, endStates);
    }
    
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::recursivelyFlagNodeDirty( const RevBayesCore::TopologyNode &n ) {

    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();

    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( dirty_nodes[index] == false )
    {
        // the root doesn't have an ancestor
        if ( n.isRoot() == false )
        {
            recursivelyFlagNodeDirty( n.getParent() );
        }

        // set the flag
        dirty_nodes[index] = true;

        // if we previously haven't touched this node, then we need to change the active likelihood pointer
        if ( changed_nodes[index] == false )
        {
            active_likelihood[index] = (active_likelihood[index] == 0 ? 1 : 0);
            changed_nodes[index] = true;
        }

    }

}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::drawStochasticCharacterMap(std::vector<std::string>& character_histories, bool set_amb_char_data, bool use_simmap_default)
{
    // first populate partial likelihood vectors along all the branches
    sample_character_history = true;
    computeLnProbability();

    size_t attempts = 0;
    bool success = false;
    while (success == false)
    {
        if (attempts == 100000)
        {
            throw RbException("After 100000 attempts a character history could not be sampled with a non-zero probability. Try increasing nTimeSlices.");
        }

        for (size_t i = 0; i < num_states; i++)
        {
            time_in_states[i] = 0.0;
        }
        
        // get the likelihoods of descendant nodes
        const TopologyNode          &root               = value->getRoot();
        size_t                       node_index         = root.getIndex();
        const TopologyNode          &left               = root.getChild(0);
        size_t                       left_index         = left.getIndex();
        const std::vector< double > &left_likelihoods   = node_partial_likelihoods[left_index][active_likelihood[left_index]];
        const TopologyNode          &right              = root.getChild(1);
        size_t                       right_index        = right.getIndex();
        const std::vector< double > &right_likelihoods  = node_partial_likelihoods[right_index][active_likelihood[right_index]];
        

        // now begin the root-to-tip pass, drawing ancestral states for each time slice conditional on the start states
        std::map<std::vector<unsigned>, double> eventMap;
        std::vector<double> speciation_rates;
        if ( use_cladogenetic_events == true )
        {
            // get cladogenesis event map (sparse speciation rate matrix)
            eventMap = cladogenesis_matrix->getValue().getEventMap();
        }
        else
        {
            speciation_rates = computeSpeciationRateAtTime( root.getAge() );
        }
        
        
        // get root frequencies
        const RbVector<double> &freqs = getRootFrequencies();
        
        std::map<std::vector<unsigned>, double> sample_probs;
        double sample_probs_sum = 0.0;
        std::map<std::vector<unsigned>, double>::iterator it;
        
        // calculate probabilities for each state
        if ( use_cladogenetic_events == true )
        {
            // iterate over each cladogenetic event possible
            // and initialize probabilities for each clado event
            for (it = eventMap.begin(); it != eventMap.end(); it++)
            {
                const std::vector<unsigned>& states = it->first;
                double speciation_rate = it->second;
                
                // we need to sample from the ancestor, left, and right states jointly,
                // so keep track of the probability of each clado event
                double prob = left_likelihoods[num_states + states[1]] * right_likelihoods[num_states + states[2]];
                prob *= freqs[states[0]] * speciation_rate;
                sample_probs[ states ] = prob;
                sample_probs_sum += prob;
            }
        }
        else
        {
            for (size_t i = 0; i < num_states; i++)
            {
                double likelihood = left_likelihoods[num_states + i] * right_likelihoods[num_states + i] * speciation_rates[i];
                std::vector<unsigned> states = boost::assign::list_of(i)(i)(i);
                sample_probs[ states ] = likelihood * freqs[i];
                sample_probs_sum += likelihood * freqs[i];
            }
        }
        
        // sample ancestor, left, and right character states from probs
        size_t a = 0, l = 0, r = 0;
        
        if (sample_probs_sum != 0)
        {
            RandomNumberGenerator* rng = GLOBAL_RNG;
            double u = rng->uniform01() * sample_probs_sum;
            
            for (it = sample_probs.begin(); it != sample_probs.end(); it++)
            {
                u -= it->second;
                if (u < 0.0)
                {
                    const std::vector<unsigned>& states = it->first;
                    a = states[0];
                    l = states[1];
                    r = states[2];
                    break;
                }
            }
        
            // save the character history for the root
            std::string simmap_string = "{" + StringUtilities::toString(a) + "," + StringUtilities::toString( root.getBranchLength() ) + "}";
            character_histories[node_index] = simmap_string;
            
            // recurse towards tips
            bool success_l = recursivelyDrawStochasticCharacterMap(left, l, character_histories, set_amb_char_data, use_simmap_default);
            bool success_r = recursivelyDrawStochasticCharacterMap(right, r, character_histories, set_amb_char_data, use_simmap_default);
            success = success_l && success_r;
        }
        
        ++attempts;
    }

    Tree t = Tree(*value);
    t.clearNodeParameters();
    t.addNodeParameter( "character_history", character_histories, false );
    simmap = t.getSimmapNewickRepresentation();
    
    // turn off sampling until we need it again
    sample_character_history = false;

}


bool EpisodicStateDependentSpeciationExtinctionFossilizationProcess::recursivelyDrawStochasticCharacterMap(const TopologyNode &node, size_t start_state, std::vector<std::string>& character_histories, bool set_amb_char_data, bool use_simmap_default)
{
    size_t node_index = node.getIndex();
    std::vector<double> speciation_rates = calculateTotalSpeciationRatePerState( node.getAge() );
    std::vector<double> extinction_rates = calculateExtinctionRatePerState( node.getAge() );
    
    // reset the number of rate-shift events
    num_shift_events[node_index] = 0;
    
    // sample characters by their probability conditioned on the branch's start state going to end states
    
    // initialize the conditional likelihoods for this branch
    std::vector< double > branch_conditional_probs = std::vector<double>(2 * num_states, 0);
    branch_conditional_probs[ num_states + start_state ] = 1.0;
    
    // first calculate extinction likelihoods via a backward time pass
    double start_time = node.getParent().getAge();
    numericallyIntegrateProcess(branch_conditional_probs, 0, start_time, true, true);
    
    // now calculate conditional likelihoods along branch in forward time
    double branch_length = node.getParent().getAge() - node.getAge();
    size_t current_dt = 0;
    double current_dt_start = 0;
    double current_dt_end = 0;
    
    size_t current_state = start_state;
    
    // set up vectors to hold the transition events
    std::vector<size_t> transition_states;
    std::vector<double> transition_times;
    transition_states.push_back(current_state);
    
    int downpass_dt = int( branch_partial_likelihoods[node_index].size() ) - 1;
   
    // keep track of rates in each time interval so we can calculate per branch averages of each rate
    double total_speciation_rate = 0.0;
    double total_extinction_rate = 0.0;
    double num_dts = 0.0;

    // loop over every time slice, stopping before the last time slice
    while ( downpass_dt >= 0 && ((current_dt + 1) * dt) < branch_length)
    {
        current_dt_start = (current_dt * dt);
        current_dt_end = ((current_dt + 1) * dt);
        
        numericallyIntegrateProcess(branch_conditional_probs, current_dt_start, current_dt_end, false, false);

        // draw state for this time slice
        size_t new_state = current_state;
        double probs_sum = 0.0;
        for (size_t i = 0; i < num_states; i++)
        {
            probs_sum += branch_conditional_probs[i + num_states] * branch_partial_likelihoods[node_index][downpass_dt][i];
        }
        if ( probs_sum == 0.0 )
        {
            return false;
        }
        else
        {
            RandomNumberGenerator* rng = GLOBAL_RNG;
            double u = rng->uniform01() * probs_sum;

            for (size_t i = 0; i < num_states; i++)
            {
                u -= branch_conditional_probs[i + num_states] * branch_partial_likelihoods[node_index][downpass_dt][i];
                if (u < 0.0)
                {
                    new_state = i;
                    break;
                }
            }
        }
        
        // check if there was a character state transition
        if (new_state != current_state)
        {
            double time_since_last_transition = 0.0;
            double transition_times_sum = 0.0;
            for (size_t j = 0; j < transition_times.size(); j++)
            {
                transition_times_sum += transition_times[j];
            }
            time_since_last_transition = current_dt_end - transition_times_sum;

            transition_times.push_back(time_since_last_transition);
            transition_states.push_back(new_state);
            current_state = new_state;
            
            ++num_shift_events[node_index];
        }
        
        // condition branch_conditional_probs on the sampled state
        for (size_t i = 0; i < num_states; i++)
        {
            if (i == current_state)
            {
                branch_conditional_probs[ num_states + i ] = 1.0;
            }
            else
            {
                branch_conditional_probs[ num_states + i ] = 0.0;
            }
        }
        
        current_dt++;
        downpass_dt--;
        
        // keep track of rates in this interal so we can calculate per branch averages of each rate
        total_speciation_rate += speciation_rates[current_state];
        total_extinction_rate += extinction_rates[current_state];
        time_in_states[current_state] += dt;
        num_dts += 1;
    }
    
    if ( node.isTip() == true )
    {
        // the last time slice of the branch will be the observed state
        
        AbstractHomologousDiscreteCharacterData& data = static_cast<TreeDiscreteCharacterData*>(this->value)->getCharacterData();
        AbstractDiscreteTaxonData& taxon_data = data.getTaxonData( node.getName() );
        
        DiscreteCharacterState &char_state = taxon_data.getCharacter(0);
        size_t new_state = current_state;
        
        if ( char_state.isAmbiguous() == false )
        {
            new_state = char_state.getStateIndex();
        }
        else
        {
            // use the simulated state
            if (set_amb_char_data == true)
            {
                // overwrite the character data
                char_state.setMissingState(false);
                char_state.setStateByIndex(new_state);
            }
        }
        
        // keep track of rates in this interval so we can calculate per branch averages of each rate
        total_speciation_rate += speciation_rates[new_state];
        total_extinction_rate += extinction_rates[new_state];
        time_in_states[new_state] += dt;
        num_dts += 1;
        
        // check if there was a character state transition
        if (new_state != current_state)
        {
            double time_since_last_transition = 0.0;
            double transition_times_sum = 0.0;
            for (size_t j = 0; j < transition_times.size(); j++)
            {
                transition_times_sum += transition_times[j];
            }
            time_since_last_transition = current_dt_end - transition_times_sum;
            
            transition_times.push_back(time_since_last_transition);
            transition_states.push_back(new_state);
            ++num_shift_events[node_index];
        }
        
        // add the length of the final character state
        double time_since_last_transition = 0.0;
        double transition_times_sum = 0.0;
        for (size_t j = 0; j < transition_times.size(); j++)
        {
            transition_times_sum += transition_times[j];
        }
        time_since_last_transition = branch_length - transition_times_sum;
        transition_times.push_back(time_since_last_transition);
        
        // make SIMMAP string
        std::string simmap_string = "{";

        if (use_simmap_default == true)
        {
            for (size_t i = transition_times.size(); i > 0; i--)
            {
                simmap_string = simmap_string + StringUtilities::toString(transition_states[i - 1]) + "," + StringUtilities::toString(transition_times[i - 1]);
                if (i != 1)
                {
                    simmap_string = simmap_string + ":";
                }
            }
        }
        else
        {
            for (size_t i = 0; i < transition_times.size(); i++)
            {
                if (i != 0)
                {
                    simmap_string = simmap_string + ":";
                }
                simmap_string = simmap_string + StringUtilities::toString(transition_states[i]) + "," + StringUtilities::toString(transition_times[i]);
            }
        }
        simmap_string = simmap_string + "}";
        
        // calculate average diversification rates on this branch
        average_speciation[node_index] = total_speciation_rate / num_dts;
        average_extinction[node_index] = total_extinction_rate / num_dts;

        // save the character history for this branch
        character_histories[node_index] = simmap_string;
        
    }
    else
    {
        // the last time slice of the branch will be the state of the node before any cladogenetic events
        
        std::map<std::vector<unsigned>, double> event_map;
        if ( use_cladogenetic_events == true )
        {
            // get cladogenesis event map (sparse speciation rate matrix)
            event_map = cladogenesis_matrix->getValue().getEventMap();
        }
        
        // get likelihoods of descendant nodes
        const TopologyNode     &left                = node.getChild(0);
        size_t                  left_index          = left.getIndex();
        std::vector< double >   left_likelihoods    = node_partial_likelihoods[left_index][active_likelihood[left_index]];
        const TopologyNode     &right               = node.getChild(1);
        size_t                  right_index         = right.getIndex();
        std::vector< double >   right_likelihoods   = node_partial_likelihoods[right_index][active_likelihood[right_index]];
        
        std::map<std::vector<unsigned>, double> sample_probs;
        double sample_probs_sum = 0.0;
        std::map<std::vector<unsigned>, double>::iterator it;
        
        // calculate probabilities for each state
        if ( use_cladogenetic_events == true )
        {
            // iterate over each cladogenetic event possible
            // and initialize probabilities for each clado event
            for (it = event_map.begin(); it != event_map.end(); it++)
            {
                const std::vector<unsigned>& states = it->first;
                double speciation_rate = it->second;
                
                // we need to sample from the ancestor, left, and right states jointly,
                // so keep track of the probability of each clado event
                double prob = left_likelihoods[num_states + states[1]] * right_likelihoods[num_states + states[2]];
                prob *= speciation_rate * branch_conditional_probs[num_states + states[0]];
                sample_probs[ states ] = prob;
                sample_probs_sum += prob;
            }
        }
        else
        {
            for (size_t i = 0; i < num_states; i++)
            {
                double prob = left_likelihoods[num_states + i] * right_likelihoods[num_states + i] * speciation_rates[i];
                prob *= branch_conditional_probs[num_states + i];
                std::vector<unsigned> states = boost::assign::list_of(i)(i)(i);
                sample_probs[ states ] = prob;
                sample_probs_sum += prob;
            }
        }
        
        // finally, sample ancestor, left, and right character states from probs
        size_t a = 0;
        size_t l = 0;
        size_t r = 0;
        
        if (sample_probs_sum == 0)
        {
            return false;
        }
        else
        {
            RandomNumberGenerator* rng = GLOBAL_RNG;
            double u = rng->uniform01() * sample_probs_sum;
            
            for (it = sample_probs.begin(); it != sample_probs.end(); it++)
            {
                u -= it->second;
                if (u < 0.0)
                {
                    const std::vector<unsigned>& states = it->first;
                    a = states[0];
                    l = states[1];
                    r = states[2];
                    break;
                }
            }
        }
        
        // keep track of rates in this interval so we can calculate per branch averages of each rate
        total_speciation_rate += speciation_rates[a];
        total_extinction_rate += extinction_rates[a];
        time_in_states[a] += dt;
        num_dts += 1;
        
        // check if there was a character state transition
        if (a != current_state)
        {
            double time_since_last_transition = 0.0;
            double transition_times_sum = 0.0;
            for (size_t j = 0; j < transition_times.size(); j++)
            {
                transition_times_sum += transition_times[j];
            }
            time_since_last_transition = current_dt_end - transition_times_sum;

            transition_times.push_back(time_since_last_transition);
            transition_states.push_back(a);
            ++num_shift_events[node_index];
        }
        
        // add the length of the final character state
        double time_since_last_transition = 0.0;
        double transition_times_sum = 0.0;
        for (size_t j = 0; j < transition_times.size(); j++)
        {
            transition_times_sum += transition_times[j];
        }
        time_since_last_transition = branch_length - transition_times_sum;

        transition_times.push_back(time_since_last_transition);
        
        // make SIMMAP string
        std::string simmap_string = "{";
        for (size_t i = transition_times.size(); i > 0; i--)
        {
            simmap_string = simmap_string + StringUtilities::toString(transition_states[i - 1]) + "," + StringUtilities::toString(transition_times[i - 1]);
            if (i != 1)
            {
                simmap_string = simmap_string + ":";
            }
        }
        simmap_string = simmap_string + "}";
        
        // save the character history for this branch
        character_histories[node_index] = simmap_string;
        
        // calculate average diversification rates on this branch
        average_speciation[node_index] = total_speciation_rate / num_dts;
        average_extinction[node_index] = total_extinction_rate / num_dts;
        
        // recurse towards tips
        bool success_l = recursivelyDrawStochasticCharacterMap(left, l, character_histories, set_amb_char_data, use_simmap_default);
        bool success_r = recursivelyDrawStochasticCharacterMap(right, r, character_histories, set_amb_char_data, use_simmap_default);
        return success_l && success_r;
    }
    return true;
}


RevLanguage::RevPtr<RevLanguage::RevVariable> EpisodicStateDependentSpeciationExtinctionFossilizationProcess::executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found)
{
    if (name == "clampCharData")
    {
        found = true;
        
        const AbstractHomologousDiscreteCharacterData& v = static_cast<const TypedDagNode<AbstractHomologousDiscreteCharacterData > *>( args[0] )->getValue();
    
        // check if the tip names match
        bool match = true;
        std::vector<string> tips = value->getTipNames();
        for (size_t i = 0; i < tips.size(); i++)
        {
            found = false;
            for (size_t j = 0; j < v.getNumberOfTaxa(); j++)
            {
                if (tips[i] == v[j].getTaxonName())
                {
                    found = true;
                    break;
                }
            }
            if (found == false)
            {
                match = false;
                break;
            }
        }
        if (match == false)
        {
            throw RbException("To clamp a character data object all taxa present in the tree must be present in the character data.");
        }
        
        static_cast<TreeDiscreteCharacterData*>(this->value)->setCharacterData( v.clone() );
   
        // Sebastian (20210519): We should not waste computations here if we actually don't need it. Try to do lazy evaluations.
        // I keep this here if we find out later that these were indeed.
        // simulate character history over the tree conditioned on the new tip data
//        size_t num_nodes = value->getNumberOfNodes();
//        std::vector<std::string> character_histories(num_nodes);
//        drawStochasticCharacterMap(character_histories);
//        static_cast<TreeDiscreteCharacterData*>(this->value)->setTimeInStates(time_in_states);

        return NULL;
    }
    
    if (name == "getCharData")
    {
        found = true;
        RevLanguage::AbstractHomologousDiscreteCharacterData *tip_states = new RevLanguage::AbstractHomologousDiscreteCharacterData( getCharacterData() );
        return new RevLanguage::RevVariable( tip_states );
    }
    if ( name == "getCharHistory" )
    {
        found = true;
        return new RevLanguage::RevVariable( new RlString( simmap ) );
    }
    return TypedDistribution<Tree>::executeProcedure( name, args, found );
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::executeMethod(const std::string &name, const std::vector<const DagNode *> &args, RbVector<std::int64_t> &rv) const
{
   
    if ( name == "numberEvents" )
    {
        rv = num_shift_events;
    }
    else
    {
        throw RbException() << "The state dependent birth-death process does not have a member method called '" << name << "'.";
    }

}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::executeMethod(const std::string &name, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
   
    if ( name == "averageSpeciationRate" )
    {
        rv = average_speciation;
    }
    else if ( name == "averageExtinctionRate" )
    {
        rv = average_extinction;
    }
    else if ( name == "getTimeInStates" )
    {
        rv = time_in_states;
    }
    else
    {
        throw RbException() << "The state dependent birth-death process does not have a member method called '" << name << "'.";
    }

}


/**
 * return the index i so that s_{i-1} <= t < s_i
 * where s_i is the global timeline of events
 * s_0 = 0.0
 * s_l = Inf
 */
size_t EpisodicStateDependentSpeciationExtinctionFossilizationProcess::findIndex(double t) const
{
    // @TODO @efficiency: this would be much faster if we can get std::lower_bound to work consistently
    // Linear search for interval because std::lower_bound is not cooperating
    if (global_timeline.size() == 1)
    {
        return (t <= global_timeline[0] ? 0 : 1);
    }
    else if ( t <= global_timeline[0] )
    {
        return 0;
    }
    else
    {
        for (size_t i=1; i < global_timeline.size(); ++i)
        {
            if (t > (global_timeline[i]-1E-5) && t <= (global_timeline[i+1]-1E-5))
            {
                return i-1;
            }
        }

        return global_timeline.size();
    }
}

/**
 * return the index i so that x_{i-1} <= t < x_i
 * where x is one of the input vector timelines
 */
size_t EpisodicStateDependentSpeciationExtinctionFossilizationProcess::findIndex(double t, const std::vector<double> &timeline) const
{

    // Linear search for interval because std::lower_bound is not cooperating
    if (timeline.size() == 1)
    {
        return (t <= timeline[0] ? 0 : 1);
    }
    else if ( t <= timeline[0] )
    {
        return 0;
    }
    else
    {
        for (size_t i=1; i < timeline.size(); ++i)
        {
            if (t > timeline[i-1] && t <= timeline[i])
            {
                return i;
            }
        }

        return timeline.size();
    }
}


/**
 * Get the affected nodes by a change of this node.
 * If the root age has changed than we need to call get affected again.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    if ( affecter == process_age )
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
}


/**
 * Get the event rate
 */
double EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getEventRate( double age ) const
{

    if ( use_episodic_model == true )
    {
        size_t index_epoch = computeEpochIndex(age);
        
        if ( index_epoch >= eta.size() ) throw RbException("Didn't rescale vector eta correctly.");

        return eta[index_epoch];
    }
    else if ( eta_const != NULL )
    {
        return eta_const->getValue();
    }
    else
    {
        return 1.0;
    }

}


/**
 * Get the event rate generator
 */
const RateGenerator& EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getEventRateMatrix(double age) const
{

    if ( use_episodic_model == true )
    {
        size_t index_epoch = computeEpochIndex(age);
        
        if ( index_epoch >= Q.size() ) throw RbException("Didn't rescale vector Q correctly.");

        return Q[index_epoch];
    }
    else
    {
        if ( Q_const != NULL )
        {
            return Q_const->getValue();
        }
        else
        {
            return Q_default;
        }
    }

}


double EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getOriginAge( void ) const
{

    return process_age->getValue();
}


std::vector<double> EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getAverageExtinctionRatePerBranch( void ) const
{
    return average_extinction;
}


std::vector<double> EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getAverageSpeciationRatePerBranch( void ) const
{
    return average_speciation;
}


std::vector<std::int64_t> EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getNumberOfShiftEventsPerBranch( void ) const
{
    return num_shift_events;
}


std::vector<double> EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getTimeInStates( void ) const
{
    return time_in_states;
}


/**
 * By default, the root age is assumed to be equal to the origin time.
 * This should be overridden if a distinct root age is needed
 */
double EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getRootAge( void ) const
{

    if (use_origin)
    {
        if (value->getNumberOfNodes() > 0)
        {
            return value->getRoot().getAge();
        }
        else
        {
            return 0;
        }
    }
    else
        return getOriginAge();
}


/**
 * Get the stationary root frequencies
 */
std::vector<double> EpisodicStateDependentSpeciationExtinctionFossilizationProcess::getRootFrequencies(void) const
{

    if ( pi != NULL )
    {
        return pi->getValue();
    }
    else
    {
        return std::vector<double>(num_states, 1.0/num_states);
    }

}

bool EpisodicStateDependentSpeciationExtinctionFossilizationProcess::isEpisodicModel(void) const
{
    bool has_interval_times = false;
    // For there to be no intervals, every timeline must either be NULL or have size 0
    if ( (epoch_times_lambda != NULL      && epoch_times_lambda->getValue().size() > 0 ) ||
         (epoch_times_mu != NULL          && epoch_times_mu->getValue().size() > 0 ) ||
         (epoch_times_phi != NULL         && epoch_times_phi->getValue().size() > 0 ) ||
         (epoch_times_gamma != NULL       && epoch_times_gamma->getValue().size() > 0 ) ||
         (epoch_times_eta != NULL         && epoch_times_eta->getValue().size() > 0 ) ||
         (epoch_times_Q != NULL           && epoch_times_Q->getValue().size() > 0 ) )
    {
        has_interval_times = true;
    }

    bool all_parameters_are_scalars = false;
    // For all parameters to be scalars,
    // 1) rate parameters must either be homogenous or they must have size <= 1 (1 for scalar, 0 if it's null)
    // 2) Lambda/Mu must be of size 0 or NULL
    // 3) Phi must be of size 1 or a scalar
    if ( (lambda_var == NULL     || lambda_var->getValue().size() <= 1)     &&
         (mu_var == NULL         || mu_var->getValue().size() <= 1)         &&
         (phi_var == NULL        || phi_var->getValue().size() <= 1)        &&
         (survival_probs == NULL || survival_probs->getValue().size() == 0) &&
         (eta_var     == NULL    || eta_var->getValue().size() == 0)        &&
         (Q_var == NULL          || Q_var->getValue().size() <= 1) )
    {
         all_parameters_are_scalars = true;
    }


    if (has_interval_times && all_parameters_are_scalars)
    {
        throw RbException("No timeline(s) was (were) provided but there are non-scalar parameters.");
    }

    return has_interval_times && !all_parameters_are_scalars;
}



/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::keepSpecialization(const DagNode *affecter)
{
    
    if ( affecter == process_age )
    {
        dag_node->keepAffected();
    }
    
    // reset all flags
    for (std::vector<bool>::iterator it = this->dirty_nodes.begin(); it != this->dirty_nodes.end(); ++it)
    {
        (*it) = false;
    }

    for (std::vector<bool>::iterator it = this->changed_nodes.begin(); it != this->changed_nodes.end(); ++it)
    {
        (*it) = false;
    }

}


double EpisodicStateDependentSpeciationExtinctionFossilizationProcess::lnProbTreeShape(void) const
{
    // the birth death divergence times density is derived for a (ranked) unlabeled oriented tree
    // so we convert to a (ranked) labeled non-oriented tree probability by multiplying by 2^{n+m-1} / n!
    // where n is the number of extant tips, m is the number of extinct tips

    int num_taxa = (int)value->getNumberOfTips();
    int num_extinct = (int)value->getNumberOfExtinctTips();
    int num_sa = (int)value->getNumberOfSampledAncestors();

    return (num_taxa - num_sa - 1) * RbConstants::LN2 - RbMath::lnFactorial(num_taxa - num_sa);
}


/*
 * Here wepopulate all parameter vectors with their final values.
 * This requires that we:
 *    1) Clear out old values of all parameter vectors
 *    2) Refill and sort vector-valued parameters (leaving scalar parameters alone) to go from present to past
 *    3) Sort (assemble first if needed) the global timeline, attach the first time (the offset)
 * Then we can fill in our final vector for each parameter, which will be a vector of the same size as the global timeline
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::prepareTimeline( void ) const
{
    // clean all the sets
    lambda.clear();
    mu.clear();
    phi.clear();
    eta.clear();
    gamma.clear();
    Q.clear();

    global_timeline.clear();

    // put in current values for vector parameters so we can re-order them as needed
    RbVector<double> empty_timeline;
    RbVector<double> lambda_times = ( epoch_times_lambda != NULL ? epoch_times_lambda->getValue() : empty_timeline );
    RbVector<double> mu_times     = ( epoch_times_mu     != NULL ? epoch_times_mu->getValue()     : empty_timeline );
    RbVector<double> phi_times    = ( epoch_times_phi    != NULL ? epoch_times_phi->getValue()    : empty_timeline );
    RbVector<double> gamma_times  = ( epoch_times_gamma  != NULL ? epoch_times_gamma->getValue()  : empty_timeline );
    RbVector<double> eta_times    = ( epoch_times_eta    != NULL ? epoch_times_eta->getValue()    : empty_timeline );
    RbVector<double> Q_times      = ( epoch_times_Q      != NULL ? epoch_times_Q->getValue()      : empty_timeline );

    // If it's a constant-rate process, make sure we only have scalars
    use_episodic_model = isEpisodicModel();
    if ( use_episodic_model == false )
    {
        global_timeline = std::vector<double>(0,0.0);
    }
    // We only need to assemble a global timeline if
    else
    {
        // check if correct number of speciation rates were provided
        // if provided as a vector, sort to the correct timescale
        if ( lambda_const == NULL && lambda_var == NULL)
        {
            throw RbException("Speciation rate must be of type RealPos or RealPos[]");
        }
        else if ( lambda_var != NULL )
        {
            if ( lambda_times.size() == 0 )
            {
                throw RbException("No time intervals provided for the piecewise constant speciation rates.");
            }
            if ( lambda_var->getValue().size() - lambda_times.size() != 1 )
            {
                throw RbException() << "Number of speciation rates (" << lambda_var->getValue().size() << ") does not match number of time intervals (" << lambda_times.size() << ")";
            }
        }

        // check if correct number of extinction rates were provided
        // if provided as a vector, sort to the correct timescale
        if ( mu_var != NULL )
        {
            if ( mu_times.size() == 0 )
            {
                throw RbException("No time intervals provided for the piecewise constant extinction rates.");
            }
            if ( mu_var->getValue().size() - mu_times.size() != 1 )
            {
                throw RbException() << "Number of extinction rates (" << mu_var->getValue().size() << ") does not match number of time intervals (" << mu_times.size() << ")";
            }
        }

        // check if correct number of fossilization rates were provided
        // if provided as a vector, sort to the correct timescale
        if ( phi_var != NULL )
        {
            if ( phi_times.size() == 0 )
            {
                throw RbException("No time intervals provided for the piecewise constant fossilization rates.");
            }
            if ( phi_var->getValue().size() - phi_times.size() != 1 )
            {
                throw RbException() << "Number of fossilization rates (" << phi_var->getValue().size() << ") does not match number of time intervals (" << phi_times.size() << ")";
            }
        }

        // check if correct number of mass extinction survival probabilities were provided
        // if provided as a vector, sort to the correct timescale
        if ( survival_probs != NULL )
        {
            if ( gamma_times.size() == 0 )
            {
                throw RbException("No time intervals provided for the mass extinction survival probabilities.");
            }
            if ( survival_probs->getValue().size() - gamma_times.size() != 0 )
            {
                throw RbException() << "Number of mass extinction survival probabilities (" << survival_probs->getValue().size() << ") does not match number of time intervals (" << gamma_times.size() << ")";
            }
        }

        // check if correct number of transition rates were provided
        // if provided as a vector, sort to the correct timescale
        if ( eta_var != NULL )
        {
            if ( eta_times.size() == 0 )
            {
                throw RbException("No time intervals provided for the transition rates.");
            }
            if ( eta_var->getValue().size() - eta_times.size() != 1 )
            {
                throw RbException() << "Number of transition rates (" << eta_var->getValue().size() << ") does not match number of time intervals (" << eta_times.size() << ")";
            }
        }

        // check if correct number of transition rate matrices were provided
        // if provided as a vector, sort to the correct timescale
        if ( Q_var != NULL )
        {
            if ( Q_times.size() == 0 )
            {
                throw RbException("No time intervals provided for the transition rate matrices.");
            }
            if ( Q_var->getValue().size() - Q_times.size() != 1 )
            {
                throw RbException() << "Number of transition rate matrices (" << Q_var->getValue().size() << ") does not match number of time intervals (" << Q_times.size() << ")";
            }
        }


        // now we start assembling the global timeline by finding the union of unique intervals for all parameters
        std::set<double> event_times;
        addTimesToGlobalTimeline(event_times, lambda_times);
        addTimesToGlobalTimeline(event_times, mu_times);
        addTimesToGlobalTimeline(event_times, phi_times);
        addTimesToGlobalTimeline(event_times, gamma_times);
        addTimesToGlobalTimeline(event_times, eta_times);
        addTimesToGlobalTimeline(event_times, Q_times);
        
        for (std::set<double>::const_iterator it = event_times.begin(); it != event_times.end(); ++it)
        {
            global_timeline.push_back( *it );
        }

        // we are done with setting up the timeline (i.e., using all the provided timelines) and checking all dimensions of parameters

    }

    // For each parameter vector, we now make sure that its size matches the size of the global vector
    // For a rate parameter, there are four cases
    //     1) It is a vector and it matches the size of the global timeline, in which case it is already sorted and we can use it
    //     2) It is a vector and it DOES NOT match the size of the global timeline, in which case we must expand it to match
    //     3) It is a scalar, in which case we simply populate a vector of the correct size with the value
    //     4) It is empty, in which case we simply populate a vector of the correct size with the default value

    // get vector of speciation rates
    if ( lambda_var != NULL )
    {
        lambda = lambda_var->getValue();
        sortNonGlobalTimesAndParameters(lambda,lambda_times);

        if ( lambda.size() != global_timeline.size() + 1)
        {
            expandNonGlobalRateParameterVector(lambda,lambda_times);
        } // else it matches in size and is already sorted and is thus ready to be used
    }
    else
    {
        lambda = std::vector<RbVector<double> >(global_timeline.size()+1,lambda_const->getValue());
    }

    // Get vector of death rates
    if ( mu_var != NULL )
    {
        mu = mu_var->getValue();
        sortNonGlobalTimesAndParameters(mu,mu_times);
        
        if ( mu.size() != global_timeline.size() + 1 )
        {
            expandNonGlobalRateParameterVector(mu,mu_times);
        } // else it matches in size and is already sorted and is thus ready to be used
    }
    else
    {
        mu = std::vector< RbVector<double> >(global_timeline.size()+1, ( mu_const == NULL ? RbVector<double>(num_states,0.0) : mu_const->getValue()) );
    }

    // Get vector of sampling rates
    if ( phi_var != NULL )
    {
        phi = phi_var->getValue();
        sortNonGlobalTimesAndParameters(phi,phi_times);
        if ( phi.size() != global_timeline.size() + 1)
        {
            expandNonGlobalRateParameterVector(phi,phi_times);
        } // else it matches in size and is already sorted and is thus ready to be used
    }
    else
    {
        RbVector<double> phi_val = ( phi_const != NULL ? phi_const->getValue() : RbVector<double>(num_states,0.0) );
        phi = std::vector< RbVector<double> >(global_timeline.size()+1, phi_val);
    }

    // For each parameter vector, we now make sure that its size matches the size of the global vector
    // For gamma, there are two cases
    //     1) It is a vector and is is of length global_timeline.size() - 1, in which case we add an event with probability 0.0 at the present, and it is ready to use
    //     2) It is a vector and it DOES NOT match the size of the global timeline, in which case we must expand it to match, which automatically adds an event of P=0.0 at the present

    // Get vector of burst birth probabilities
    if ( survival_probs != NULL )
    {
        gamma = survival_probs->getValue();
        sortNonGlobalTimesAndParameters(gamma,gamma_times);
        // Expand if needed
        if (gamma_times.size() != global_timeline.size())
        {
            expandNonGlobalProbabilityParameterVector(gamma, gamma_times, RbVector<double>(num_states,1.0));
        }
    }
    else
    {
        // User specified nothing, there are no birth bursts
        gamma = std::vector<RbVector<double> >(global_timeline.size(), RbVector<double>(num_states,1.0) );
    }

    // Get vector of transition rates
    if ( eta_var != NULL )
    {
        eta = eta_var->getValue();
        sortNonGlobalTimesAndParameters(eta,eta_times);
        // Expand if needed
        if (eta_times.size() != global_timeline.size() + 1)
        {
            expandNonGlobalRateParameterVector(eta,eta_times);
        }
    }
    else if ( eta_const != NULL )
    {
        // User specified nothing
         eta = std::vector<double>(global_timeline.size()+1,eta_const->getValue());
    }
    else
    {
        // User specified nothing
         eta = std::vector<double>(global_timeline.size()+1,1.0);
    }
    
    // Get vector of transition rates
    if ( Q_var != NULL )
    {
        Q = Q_var->getValue();
//        sortNonGlobalTimesAndParameters(Q,Q_times);
        // Expand if needed
        if (Q_times.size() != global_timeline.size() + 1)
        {
//            expandNonGlobalRateParameterVector(Q,Q_times);
        }
    }
    else if ( Q_const != NULL )
    {
        // User specified nothing
         Q = RbVector<RateGenerator>(global_timeline.size()+1,Q_const->getValue() );
    }
    else
    {
        // User specified nothing
         Q = RbVector<RateGenerator>(global_timeline.size()+1,Q_default );
    }

    
//    std::cerr << "Timeline:\t\t" << global_timeline << std::endl;
//    std::cerr << "Lambda:\t\t";
//    for (size_t i=0; i<lambda.size(); ++i) std::cerr << lambda[i][0] << " ";
//    std::cerr << std::endl;
//    std::cerr << "Mu:\t\t\t";
//    for (size_t i=0; i<mu.size(); ++i) std::cerr << mu[i][0] << " ";
//    std::cerr << std::endl;
//    std::cerr << "Phi:\t\t";
//    for (size_t i=0; i<phi.size(); ++i) std::cerr << phi[i][0] << " ";
//    std::cerr << std::endl;
}



std::vector<double> EpisodicStateDependentSpeciationExtinctionFossilizationProcess::pExtinction(double start, double end) const
{
    
    
    std::vector<double> sampling_probability;
    if ( rho != NULL && rho_per_state == NULL )
    {
        sampling_probability   = std::vector<double>(num_states, rho->getValue());
    }
    else if ( rho == NULL && rho_per_state != NULL )
    {
        sampling_probability   = rho_per_state->getValue();
    }
    else
    {
        throw RbException("Either a global sampling fraction or state-specific sampling fraction needs to be set.");
    }

    std::vector< double > initial_state = std::vector<double>(2*num_states,0);
    for (size_t i=0; i<num_states; ++i)
    {
        initial_state[i] = 1.0 - sampling_probability[i];
        initial_state[num_states + i] = sampling_probability[i];
    }
    
    numericallyIntegrateProcess(initial_state, start, end, true, false);
    
    return initial_state;
}


double EpisodicStateDependentSpeciationExtinctionFossilizationProcess::pSurvival(double start, double end) const
{

    // delegate to specific function that manages survival between origin and root
    return pSurvival(start, end, use_origin == false);
}


double EpisodicStateDependentSpeciationExtinctionFossilizationProcess::pSurvival(double start, double end, bool speciation) const
{
    

    std::vector< double > initial_state = pExtinction(start,end);
    std::vector<double>   speciation_rates = calculateTotalSpeciationRatePerState( start );

    double prob = 0.0;
    const RbVector<double> &freqs = getRootFrequencies();
    for (size_t i=0; i<num_states; ++i)
    {
        // we need to check if we should condition on survival of the speciation event
        if ( speciation == true )
        {
            prob += freqs[i]*(1.0-initial_state[i])*(1.0-initial_state[i])*speciation_rates[i];
        }
        else
        {
            prob += freqs[i]*(1.0-initial_state[i]);
        }
        
    }
    
    return prob;
}



/**
 * Redraw the current value. We delegate this to the simulate method.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::redrawValue( void )
{

    size_t attempts = 0;
    //while (attempts < 100000)
    while (attempts < 10000)
    {
        bool success = false;

        if ( condition_on_tree == true && value->getNumberOfTips() > 0 )
        {
            // simulate a character history conditioned on the observed tree

            // make character data objects -- all unknown/missing
            std::vector<string> tips = value->getTipNames();
            HomologousDiscreteCharacterData<NaturalNumbersState> *tip_data = new HomologousDiscreteCharacterData<NaturalNumbersState>();
            for (size_t i = 0; i < tips.size(); i++)
            {
                DiscreteTaxonData<NaturalNumbersState> this_tip_data = DiscreteTaxonData<NaturalNumbersState>(tips[i]);
                NaturalNumbersState state = NaturalNumbersState(0, num_states);
                state.setState("?");
                this_tip_data.addCharacter(state);
                tip_data->addTaxonData(this_tip_data);
            }
            static_cast<TreeDiscreteCharacterData*>(this->value)->setCharacterData(tip_data);
           
            // simulate character history over the new tree
            size_t num_nodes = value->getNumberOfNodes();
            if (num_nodes > 2)
            {
                std::vector<std::string> character_histories(num_nodes);
                drawStochasticCharacterMap(character_histories, true);
            }
            static_cast<TreeDiscreteCharacterData*>(this->value)->setTimeInStates(time_in_states);

            success = true;
        }
        else if ( condition_on_tip_states == true && static_cast<TreeDiscreteCharacterData *>(this->value)->hasCharacterData() == true )
        {
            success = simulateTreeConditionedOnTips(attempts);
        }
        else
        {
            success = simulateTree(attempts);
        }

        if (success == true)
        {
            return;
        }
        ++attempts;
    }
    throw RbException("After 100000 attempts a character-dependent birth death tree could not be simulated. Try changing minNumLineages or maxNumLineages.");
}



/**
 * Restore the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::restoreSpecialization(const DagNode *affecter)
{
    
    if ( affecter == process_age )
    {
        if ( use_origin == false )
        {
            value->getRoot().setAge( process_age->getValue() );
        }

        if ( dag_node != NULL )
        {
            dag_node->restoreAffected();
        }
    }
    
    // reset the flags
    for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
    {
        (*it) = false;
    }

    // restore the active likelihoods vector
    for (size_t index = 0; index < changed_nodes.size(); ++index)
    {
        // we have to restore, that means if we have changed the active likelihood vector
        // then we need to revert this change
        if ( changed_nodes[index] == true )
        {
            active_likelihood[index] = (active_likelihood[index] == 0 ? 1 : 0);
        }

        // set all flags to false
        changed_nodes[index] = false;
    }

}



void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setCladogenesisMatrix(const TypedDagNode< CladogeneticSpeciationRateMatrix >* cm)
{
    
    // remove the old parameter first
    this->removeParameter( cladogenesis_matrix );
    
    // set the value
    cladogenesis_matrix = cm;
    
    // should we use the event map for the speciation rates?
    use_cladogenetic_events = true;
    
    // add the new parameter
    this->addParameter( cladogenesis_matrix );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setExtinctionRates(const TypedDagNode< RbVector<double> >* r)
{

    // remove the old parameter first
    this->removeParameter( mu_const );
    this->removeParameter( mu_var );
    this->removeParameter( epoch_times_mu );

    // set the value
    mu_const = r;
    mu_var   = NULL;
    epoch_times_mu = NULL;

    // add the new parameter
    this->addParameter( mu_const );
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setExtinctionRates(const TypedDagNode< RbVector< RbVector<double> > >* r, const TypedDagNode<RbVector<double> >* t)
{

    // remove the old parameter first
    this->removeParameter( mu_const );
    this->removeParameter( mu_var );
    this->removeParameter( epoch_times_mu );

    // set the value
    mu_var   = r;
    mu_const = NULL;
    epoch_times_mu = t;
    
    // add the new parameter
    this->addParameter( mu_var );
    this->addParameter( epoch_times_mu );
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setFossilizationRates(const TypedDagNode< RbVector<double> >* r)
{

    // remove the old parameter first
    this->removeParameter( phi_const );
    this->removeParameter( phi_var );
    this->removeParameter( epoch_times_phi );

    // set the value
    phi_const = r;
    phi_var   = NULL;
    epoch_times_phi = NULL;

    // add the new parameter
    this->addParameter( phi_const );
    this->addParameter( epoch_times_phi );

}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setFossilizationRates(const TypedDagNode< RbVector< RbVector<double> > >* r, const TypedDagNode<RbVector<double> >* t)
{

    // remove the old parameter first
    this->removeParameter( phi_const );
    this->removeParameter( phi_var );
    this->removeParameter( epoch_times_phi );

    // set the value
    phi_var   = r;
    phi_const = NULL;
    epoch_times_phi = t;

    use_episodic_model = true;

    // add the new parameter
    this->addParameter( phi_var );
    this->addParameter( epoch_times_phi );

}

void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setMassExtinctionSurvivalProbabilities(const TypedDagNode<RbVector<RbVector<double> > > *p, const TypedDagNode<RbVector<double> >* t)
{
    // remove the old parameter first
    this->removeParameter( survival_probs );
    this->removeParameter( epoch_times_gamma );

    // set the value
    survival_probs = p;
    epoch_times_gamma = t;
    
    // add the new parameter
    this->addParameter( survival_probs );
    this->addParameter( epoch_times_gamma );

}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setSampleCharacterHistory(bool sample_history)
{
    sample_character_history = sample_history;
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setSamplingFraction(const TypedDagNode<double> *f)
{
    
    // remove the old parameter first
    this->removeParameter( rho );
    this->removeParameter( rho_per_state );

    
    // set the value
    rho = f;
    rho_per_state = NULL;
    
    // add the new parameter
    this->addParameter( rho );
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setSamplingFraction(const TypedDagNode< RbVector<double> > *f)
{
    
    // remove the old parameter first
    this->removeParameter( rho );
    this->removeParameter( rho_per_state );

    
    // set the value
    rho_per_state = f;
    rho = NULL;
    
    // add the new parameter
    this->addParameter( rho_per_state );
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setSpeciationRates(const TypedDagNode< RbVector<double> >* r)
{
    
    // remove the old parameter first
    this->removeParameter( lambda_const );
    this->removeParameter( lambda_var );
    this->removeParameter( epoch_times_lambda );

    // set the value
    lambda_const = r;
    lambda_var   = NULL;
    epoch_times_lambda = NULL;

    // should we use the event map for the speciation rates?
    use_cladogenetic_events = false;
    
    // add the new parameter
    this->addParameter( lambda_const );
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setSpeciationRates(const TypedDagNode< RbVector< RbVector<double> > >* r, const TypedDagNode<RbVector<double> >* t)
{
    
    // remove the old parameter first
    this->removeParameter( lambda_const );
    this->removeParameter( lambda_var );
    removeParameter( epoch_times_lambda );

    // set the value
    lambda_var   = r;
    lambda_const = NULL;
    epoch_times_lambda = t;

    // should we use the event map for the speciation rates?
    use_cladogenetic_events = false;
    
    // add the new parameter
    this->addParameter( lambda_var );
    this->addParameter( epoch_times_lambda );

}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setTransitionRate(const TypedDagNode<double> *r)
{
    
    // remove the old parameter first
    this->removeParameter( Q_const );
    this->removeParameter( Q_var );
    this->removeParameter( eta_const );
    this->removeParameter( eta_var );
    this->removeParameter( epoch_times_Q );
    this->removeParameter( epoch_times_eta );

    // set the value
    eta_const = r;
    eta_var   = NULL;
    Q_const   = NULL;
    Q_var     = NULL;
    epoch_times_eta = NULL;
    epoch_times_Q = NULL;

    // add the new parameter
    this->addParameter( eta_const );
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setTransitionRate(const TypedDagNode< RbVector<double> > *r, const TypedDagNode<RbVector<double> >* t)
{
    
    // remove the old parameter first
    this->removeParameter( Q_const );
    this->removeParameter( Q_var );
    this->removeParameter( eta_const );
    this->removeParameter( eta_var );
    this->removeParameter( epoch_times_Q );
    this->removeParameter( epoch_times_eta );

    // set the value
    eta_const = NULL;
    eta_var   = r;
    Q_const   = NULL;
    Q_var     = NULL;
    epoch_times_eta = t;
    epoch_times_Q = NULL;

    // add the new parameter
    this->addParameter( eta_var );
    this->addParameter( epoch_times_eta );
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setTransitionRateMatrix(const TypedDagNode<RateGenerator> *m)
{
    
    // remove the old parameter first
    this->removeParameter( Q_const );
    this->removeParameter( Q_var );
    this->removeParameter( eta_const );
    this->removeParameter( eta_var );
    this->removeParameter( epoch_times_Q );
    this->removeParameter( epoch_times_eta );

    // set the value
    eta_const = NULL;
    eta_var   = NULL;
    Q_const   = m;
    Q_var     = NULL;
    epoch_times_eta = NULL;
    epoch_times_Q = NULL;

    // add the new parameter
    this->addParameter( Q_const );
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setTransitionRateMatrix(const TypedDagNode<RbVector<RateGenerator> > *m, const TypedDagNode<RbVector<double> >* t)
{
    
    // remove the old parameter first
    this->removeParameter( Q_const );
    this->removeParameter( Q_var );
    this->removeParameter( eta_const );
    this->removeParameter( eta_var );
    this->removeParameter( epoch_times_Q );
    this->removeParameter( epoch_times_eta );

    // set the value
    eta_const = NULL;
    eta_var   = NULL;
    Q_const   = NULL;
    Q_var     = m;
    epoch_times_eta = NULL;
    epoch_times_Q = t;
    
    // add the new parameter
    this->addParameter( Q_var );
    this->addParameter( epoch_times_Q );
}


void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setNumberOfTimeSlices( double n )
{
    
    NUM_TIME_SLICES = n;
    dt = process_age->getValue() / NUM_TIME_SLICES;
    
}


/**
 * Set the current value.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::setValue(Tree *v, bool f )
{
    
    std::vector<Taxon> taxa = v->getTaxa();
    RevBayesCore::Tree *newv = TreeUtilities::startingTreeInitializer(*v, taxa, age_check_precision);
//    AbstractRootedTreeDistribution::setValue(newv, f);
        
    if (newv->isBinary() == false)
    {
        throw RbException("The character-dependent birth death process is only implemented for binary trees.");
    }

    value->getTreeChangeEventHandler().removeListener( this );

    // delegate to super class
    static_cast<TreeDiscreteCharacterData *>(this->value)->setTree( *newv );

    resizeVectors(newv->getNumberOfNodes());
    
    // clear memory
    delete v;
    v = NULL;
    
    value->getTreeChangeEventHandler().addListener( this );
    
    if ( process_age != NULL && use_origin == false )
    {
        const StochasticNode<double> *stoch_process_age = dynamic_cast<const StochasticNode<double>* >(process_age);
        if ( stoch_process_age != NULL )
        {
            const_cast<StochasticNode<double> *>(stoch_process_age)->setValue( new double( value->getRoot().getAge() ), f);
        }
        else
        {
            value->getRoot().setAge( process_age->getValue() );
        }
        
    }

    // make character data objects -- all unknown/missing
    std::vector<string> tips = value->getTipNames();
    HomologousDiscreteCharacterData<NaturalNumbersState> *tip_data = new HomologousDiscreteCharacterData<NaturalNumbersState>();
    for (size_t i = 0; i < tips.size(); i++)
    {
        DiscreteTaxonData<NaturalNumbersState> this_tip_data = DiscreteTaxonData<NaturalNumbersState>(tips[i]);
        NaturalNumbersState state = NaturalNumbersState(0, num_states);
        state.setState("?");
        this_tip_data.addCharacter(state);
        tip_data->addTaxonData(this_tip_data);
    }
    static_cast<TreeDiscreteCharacterData*>(this->value)->setCharacterData(tip_data);
    
    // Sebastian (20210519): We should not waste computations here if we actually don't need it. Try to do lazy evaluations.
    // I keep this here if we find out later that these were indeed.
    // simulate character history over the new tree
//    size_t num_nodes = value->getNumberOfNodes();
//    if (num_nodes > 2)
//    {
//        std::vector<std::string> character_histories(num_nodes);
//        drawStochasticCharacterMap(character_histories);
//    }
    static_cast<TreeDiscreteCharacterData*>(this->value)->setTimeInStates(time_in_states);
}


bool EpisodicStateDependentSpeciationExtinctionFossilizationProcess::simulateTreeConditionedOnTips( size_t attempts )
{

    if ( use_cladogenetic_events == true )
    {
        throw RbException("Simulations conditioned on the tip states are not yet implemented for cladogenetic SSE models.");
    }
    if ( prune_extinct_lineages == false )
    {
        throw RbException("Simulations conditioned on the tip states are currently implemented only when pruneExtinctLineages is set to true.");
    }
    if ( use_episodic_model == true )
    {
        throw RbException("Simulations conditioned on the tip states are currently implemented only when rates are constant over time.");
    }
    
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // a vector keeping track of the lineages currently surviving in each state
    // as we simulate forward in time
    std::vector< std::vector<size_t> > lineages_in_state = std::vector< std::vector<size_t> >(num_states, std::vector<size_t>());
    std::vector< std::vector<size_t> > extinct_lineages_in_state = std::vector< std::vector<size_t> >(num_states, std::vector<size_t>());

    // CharacterData object to hold the tip states
    const AbstractHomologousDiscreteCharacterData& tip_data = static_cast<TreeDiscreteCharacterData*>(this->value)->getCharacterData();
    if ( tip_data.getNumberOfTaxa() < 2 )
    {
        throw RbException("Simulations conditioned on the tip states require at least two extant lineages.");
    }

    // vectors keeping track of the total rate of all
    // speciation/anagenetic/extinction events for each state
    const RateGenerator *rate_matrix = &getEventRateMatrix( 0 );
    std::vector<double> extinction_rates = calculateExtinctionRatePerState( 0.0 );
    std::vector<double> total_speciation_rates = calculateTotalSpeciationRatePerState( 0.0 );
    std::vector<double> total_anagenetic_rates = calculateTotalAnageneticRatePerState( 0.0 );
    std::vector<double> r = std::vector<double>(num_states, 0);

    // create a vector of nodes for our simulated tree
    std::vector<TopologyNode*> nodes;
    Tree *sim_tree = new Tree();
    
    // make nodes for each observed tip state
    double t = 0.0;
    for (size_t i = 0; i < tip_data.getNumberOfTaxa(); ++i)
    {
        TopologyNode* tip_node = new TopologyNode(i);
        std::string tip_name = tip_data.getTaxa()[i].getName();
        tip_node->setName(tip_name);
        size_t state_index = 0;
        if (tip_data.getTaxonData(tip_name)[0].isAmbiguous() == false)
        {
            state_index = tip_data.getTaxonData(tip_name)[0].getStateIndex();
        }
        else
        {
            // state is ambigious so sample one of the observed states randomly
            double num_observed_states = tip_data.getTaxonData(tip_name)[0].getNumberObservedStates();
            if (num_observed_states > 0)
            {
                double u = rng->uniform01() * num_observed_states;
                for (size_t j = 0; j < num_states; ++j)
                {
                    if (tip_data.getTaxonData(tip_name)[0].isStateSet(j) == true)
                    {
                        --u;
                        if (u < 0)
                        {
                            state_index = j;
                            break;
                        }
                    }
                }
            }
            else
            {
                double u = rng->uniform01() * num_states;
                for (size_t j = 0; j < num_states; ++j)
                {
                    --u;
                    if (u < 0)
                    {
                        state_index = j;
                        break;
                    }
                }
            }
        }
        
        tip_node->setAge(t);
        tip_node->setTimeInStates(std::vector<double>(num_states, 0.0));
        tip_node->setNumberOfShiftEvents( 0 );
        lineages_in_state[state_index].push_back(i);
        nodes.push_back(tip_node);
    }
    
    // simulate moving backwards in time
    while (true) {

        // scale extinction as a function of time so simulations can't go back in time forever....
        if (true)
        {
            for (size_t i = 0; i < num_states; ++i)
            {
                extinction_rates[i] -= t/10; // TODO make this a user option
                if (extinction_rates[i] < 0)
                {
                    extinction_rates[i] = 0.0;
                }
            }
        }

        // calculate c and g from Hua and Bromham 2016
        for (size_t i = 0; i < num_states; ++i)
        {
            r[i] = 0.0;
            if (lineages_in_state[i].size() > 1)
            {
                r[i] += total_speciation_rates[i];
            }
            if (lineages_in_state[i].size() > 0)
            {
                r[i] += extinction_rates[i] + total_anagenetic_rates[i];
            }
        }
        double g = 0;
        double c = 0;
        for (size_t i = 0; i < num_states; ++i)
        {
            if (lineages_in_state[i].size() > 0)
            {
                g += r[i] * (lineages_in_state[i].size() - 1);
            }
            c += r[i] * (lineages_in_state[i].size());
        }
        if (g == 0)
        {
            g = c;
        }
        
        // use rejection sampling to sample a time for the next event
        double dt = 0.0;
        std::vector<double> prob_speciation = std::vector<double>(num_states, 0);
        std::vector<double> prob_extinction = std::vector<double>(num_states, 0);
        std::vector< std::vector<double> > prob_transition = std::vector< std::vector<double> >(num_states, std::vector<double>(num_states, 0));
        std::vector<double> prob_transition_sum = std::vector<double>(num_states, 0);
        std::vector<double> prob_state = std::vector<double>(num_states, 0);
        double prob_sum = 0.0;
        size_t tries = 0;
        while (true)
        {
        
            // propose a new time from proposal distribution g
            dt = RbStatistics::Exponential::rv( g, *rng );

            // calculate probability for the new time
            for (size_t i = 0; i < num_states; ++i)
            {
                double total_rate_spec = 0.0;
                double total_rate_ext = 0.0;
                std::vector<double> total_rate_ana = std::vector<double>(num_states, 0);
                for (size_t j = 0; j < num_states; ++j)
                {
                    if (i == j)
                    {
                        total_rate_spec += r[j] * (lineages_in_state[j].size() - 1);
                        total_rate_ext += r[j] * (lineages_in_state[j].size() + 1);
                    }
                    else
                    {
                        total_rate_spec += r[j] * lineages_in_state[j].size();
                        total_rate_ext += r[j] * lineages_in_state[j].size();

                        // the total rates for the transition from i into j
                        for (size_t k = 0; k < num_states; ++k)
                        {
                            if (k == i)
                            {
                                total_rate_ana[j] += r[k] * (lineages_in_state[k].size() + 1);
                            }
                            else if (k == j)
                            {
                                total_rate_ana[j] += r[k] * (lineages_in_state[k].size() - 1);
                            }
                            else
                            {
                                total_rate_ana[j] += r[k] * lineages_in_state[k].size();
                            }
                        }
                    }
                }

                if (lineages_in_state[i].size() > 1)
                {
                    prob_speciation[i] = total_speciation_rates[i] * (lineages_in_state[i].size() - 1) * exp(-1 * dt * total_rate_spec);
                }
                prob_extinction[i] = extinction_rates[i] * (lineages_in_state[i].size() + 1) * exp(-1 * dt * total_rate_ext);

                for (size_t j = 0; j < num_states; ++j)
                {
                    if (i != j && lineages_in_state[j].size() > 0)
                    {
                        prob_transition[i][j] = rate_matrix->getRate(i, j, 0.0, getEventRate( t )) * (lineages_in_state[i].size() + 1) * exp(-1 * dt * total_rate_ana[j]);
                        prob_transition_sum[i] += prob_transition[i][j];
                    }
                }

                prob_state[i] = prob_speciation[i] + prob_extinction[i] + prob_transition_sum[i];
                prob_sum += prob_state[i];
            }

            // check if we accept the new time
            double u = rng->uniform01();
            if (u <= prob_sum / (c * exp(-1 * dt * g)))
            {
                break;
            }

            // otherwise reinitialize and try again
            prob_speciation = std::vector<double>(num_states, 0);
            prob_extinction = std::vector<double>(num_states, 0);
            prob_transition = std::vector< std::vector<double> >(num_states, std::vector<double>(num_states, 0));
            prob_transition_sum = std::vector<double>(num_states, 0);
            prob_state = std::vector<double>(num_states, 0);
            prob_sum = 0.0;
       
            tries++;
            if (tries == 100)
            {
                return false;
            }
        }

        t = t + dt;

        // stop and retry if lineages didn't coalesce in time
        if (t > max_time)
        {
            delete sim_tree;
            nodes.clear();
            return false;
        }
      
        // extend all current branches to the new time
        for (size_t i = 0; i < num_states; ++i)
        {
            for (size_t j = 0; j < lineages_in_state[i].size(); ++j)
            {
                size_t idx = lineages_in_state[i][j];
                std::vector<double> state_times = nodes[idx]->getTimeInStates();
                state_times[i] += dt;
                nodes[idx]->setTimeInStates(state_times);
            }
        }
        
        // determine the state for the event that occurred
        size_t event_state = 0;
        double u = rng->uniform01() * prob_sum;
        for (size_t i = 0; i < num_states; ++i)
        {
            u -= prob_state[i];
            if (u < 0)
            {
                event_state = i;
                break;
            }
        }

        // determine the type of event
        std::string event_type = "";
        u = rng->uniform01() * prob_state[event_state];
        while (true) {
            u = u - prob_extinction[event_state];
            if (u < 0)
            {
                event_type = "extinction";
                break;
            }
            u = u - prob_speciation[event_state];
            if (u < 0)
            {
                event_type = "speciation";
                break;
            }
            u = u - prob_transition_sum[event_state];
            if (u < 0)
            {
                event_type = "anagenetic";
                break;
            }
        }

        if (event_type == "extinction")
        {

            size_t node_index = nodes.size();
            TopologyNode* e = new TopologyNode(node_index);
            e->setAge(t);
            e->setTimeInStates(std::vector<double>(num_states, 0.0));
            std::stringstream ss;
            ss << "ex" << node_index;
            std::string name = ss.str();
            e->setName(name);
            extinct_lineages_in_state[event_state].push_back(node_index);
            lineages_in_state[event_state].push_back(node_index);
            nodes.push_back(e);

        }
        
        if (event_type == "anagenetic")
        {
            // sample new state to transition to
            size_t new_state = 0;
            double u = rng->uniform01() * prob_transition_sum[event_state];
            for (size_t i = 0; i < num_states; i++)
            {
                u -= prob_transition[event_state][i];
                if (u < 0)
                {
                    new_state = i;
                    break;
                }
            }

            // determine which lineage gets the event
            size_t node_index = 0;
            u = rng->uniform01() * static_cast<double>(lineages_in_state[new_state].size());
            node_index = lineages_in_state[new_state][floor(u)];
            
            // remove this lineage from the new state and add it to old state
            lineages_in_state[new_state].erase(std::remove(lineages_in_state[new_state].begin(), lineages_in_state[new_state].end(), node_index), lineages_in_state[new_state].end());
            lineages_in_state[event_state].push_back(node_index);
            
            // increment the shift counter
            nodes[node_index]->setNumberOfShiftEvents( nodes[node_index]->getNumberOfShiftEvents() + 1 );

        }
        
        if (event_type == "speciation")
        {
            // pick two daughter lineages
            size_t daughter1 = 0;
            double u = rng->uniform01() * static_cast<double>(lineages_in_state[event_state].size());
            daughter1 = lineages_in_state[event_state][floor(u)];
            size_t daughter2 = daughter1;
            while (daughter1 == daughter2)
            {
                u = rng->uniform01() * static_cast<double>(lineages_in_state[event_state].size());
                daughter2 = lineages_in_state[event_state][floor(u)];
            }
            
            // check to see if this is the root
            size_t num_lineages = 0;
            for (size_t i = 0; i < num_states; ++i)
            {
                num_lineages += lineages_in_state[i].size();
            }
            bool is_root = false;
            if (num_lineages == 2)
            {
                is_root = true;
            }

            // make node for parent
            size_t parent_index = nodes.size();
            TopologyNode* p = new TopologyNode(parent_index);
            p->setAge(t);
            p->setTimeInStates(std::vector<double>(num_states, 0.0));
            p->setNumberOfShiftEvents(0);
            p->addChild(nodes[daughter1]);
            p->addChild(nodes[daughter2]);
            nodes[daughter1]->setParent(p);
            nodes[daughter2]->setParent(p);
            lineages_in_state[event_state].push_back(parent_index);
            nodes.push_back(p);

            // remove the children nodes from the vector of current lineages
            lineages_in_state[event_state].erase(std::remove(lineages_in_state[event_state].begin(), lineages_in_state[event_state].end(), daughter1), lineages_in_state[event_state].end());
            lineages_in_state[event_state].erase(std::remove(lineages_in_state[event_state].begin(), lineages_in_state[event_state].end(), daughter2), lineages_in_state[event_state].end());

            if (is_root == true)
            {
                sim_tree->setRoot(p, true);
                sim_tree->setRooted(true);
                break;
            }
        }
    }
  
    // prune extinct lineage if necessary
    if (prune_extinct_lineages == true)
    {
        for (size_t i = 0; i < num_states; ++i)
        {
            for (size_t j = 0; j < extinct_lineages_in_state[i].size(); ++j)
            {
                size_t this_node = extinct_lineages_in_state[i][j];
                if (nodes[this_node]->isTip() == true)
                {
                    sim_tree->dropTipNodeWithName( nodes[this_node]->getName() );
                }
            }
        }
    }
    
    // update character history vectors
    resizeVectors(sim_tree->getNumberOfNodes());
    simmap = "";
    for (size_t i = 0; i < sim_tree->getNumberOfNodes(); ++i)
    {
        double branch_total_speciation = 0.0;
        double branch_total_extinction = 0.0;
        for (size_t j = 0; j < num_states; ++j)
        {
            time_in_states[j] += sim_tree->getNodes()[i]->getTimeInStates()[j];
            branch_total_speciation += sim_tree->getNodes()[i]->getTimeInStates()[j] * total_speciation_rates[j];
            branch_total_extinction += sim_tree->getNodes()[i]->getTimeInStates()[j] * extinction_rates[j];
        }
        if (sim_tree->getNodes()[i]->getBranchLength() > 0)
        {
            average_speciation[i] = branch_total_speciation/sim_tree->getNodes()[i]->getBranchLength();
            average_extinction[i] = branch_total_extinction/sim_tree->getNodes()[i]->getBranchLength();
            num_shift_events[i]   = sim_tree->getNodes()[i]->getNumberOfShiftEvents();
        }
    }
    
    // set the simulated values
    value->getTreeChangeEventHandler().removeListener( this );
    static_cast<TreeDiscreteCharacterData *>(this->value)->setTree( *sim_tree );
    delete sim_tree;
    nodes.clear();
    value->getTreeChangeEventHandler().addListener( this );
    static_cast<TreeDiscreteCharacterData*>(this->value)->setTimeInStates(time_in_states);
    return true;
    
}


/**
 *
 */
bool EpisodicStateDependentSpeciationExtinctionFossilizationProcess::simulateTree( size_t attempts )
{

    // prepare the timelines for simulation
    prepareTimeline();
    
    if ( use_origin == true && condition_on_num_tips == false )
    {
        // if originAge is set we start with one lineage
        // if rootAge is set we start with two lineages and their speciation event
        throw RbException("Simulations are currently only implemented when rootAge is set. You set the originAge.");
    }
    
    if (exact_num_lineages < 2 && condition_on_num_tips == true)
    {
        throw RbException("When simulating conditioned on the number of tips exactNumLineages must be 2 or more.");
    }
    
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // a vector keeping track of the lineages currently surviving in each state
    // as we simulate forward in time
    std::vector< std::vector<size_t> > lineages_in_state = std::vector< std::vector<size_t> >(num_states, std::vector<size_t>());
    std::vector< std::vector<size_t> > extinct_lineages_in_state = std::vector< std::vector<size_t> >(num_states, std::vector<size_t>());

    // CharacterData object to hold the tip states
    HomologousDiscreteCharacterData<NaturalNumbersState> *tip_data = new HomologousDiscreteCharacterData<NaturalNumbersState>();

    // vectors keeping track of the total rate of all
    // cladogenetic/anagenetic/extinction events for each state
    std::vector<double> extinction_rates = calculateExtinctionRatePerState( process_age->getValue() );
    std::vector<double> total_speciation_rates = calculateTotalSpeciationRatePerState( process_age->getValue() );
    std::vector<double> total_anagenetic_rates = calculateTotalAnageneticRatePerState( process_age->getValue() );
    std::vector<double> total_rate_for_state = std::vector<double>(num_states, 0);
    for (size_t i = 0; i < num_states; i++)
    {
        total_rate_for_state[i] = extinction_rates[i] + total_speciation_rates[i] + total_anagenetic_rates[i];
    }

    // get the speciation rates, extinction rates, and Q matrix
    std::map<std::vector<unsigned>, double> eventMap;
    std::vector<double> speciation_rates;
    std::map<std::vector<unsigned>, double>::iterator it;
    if ( use_cladogenetic_events == true )
    {
        eventMap = cladogenesis_matrix->getValue().getEventMap();
    }
    else
    {
        speciation_rates = calculateTotalSpeciationRatePerState( process_age->getValue() );
    }
    const RateGenerator *rate_matrix = &getEventRateMatrix( process_age->getValue() );

    // a vector of all nodes in our simulated tree
    std::vector<TopologyNode*> nodes;

    // initialize the root node
    TopologyNode* root = new TopologyNode();
    double t = process_age->getValue();
    if (condition_on_num_tips == true)
    {
        t = 0.0;
    }
    root->setAge(t);
    root->setTimeInStates(std::vector<double>(num_states, 0.0));
    root->setNumberOfShiftEvents(0);
    nodes.push_back(root);

    // now draw a state for the root cladogenetic event
    
    // get root frequencies
    const RbVector<double> &root_freqs = getRootFrequencies();
    
    std::map<std::vector<unsigned>, double> sample_probs;
    double sample_probs_sum = 0.0;
    
    // calculate probabilities for each state
    if ( use_cladogenetic_events == true )
    {
        // iterate over each cladogenetic event possible
        // and initialize probabilities for each clado event
        for (it = eventMap.begin(); it != eventMap.end(); it++)
        {
            const std::vector<unsigned>& states = it->first;
            double speciation_rate = it->second;
            
            // we need to sample from the ancestor, left, and right states jointly,
            // so keep track of the probability of each clado event
            double prob = root_freqs[states[0]] * speciation_rate;
            sample_probs[ states ] = prob;
            sample_probs_sum += prob;
        }
    }
    else
    {
        for (size_t i = 0; i < num_states; i++)
        {
            std::vector<unsigned> states = boost::assign::list_of(i)(i)(i);
            sample_probs[ states ] = root_freqs[i] * speciation_rates[i];
            sample_probs_sum += root_freqs[i] * speciation_rates[i];
        }
    }
    
    // sample left and right character states from probs
    size_t l = 0, r = 0;
    
    if (sample_probs_sum == 0)
    {
        size_t u = rng->uniform01() * sample_probs.size();
        size_t v = 0;
        for (it = sample_probs.begin(); it != sample_probs.end(); it++)
        {
            if (u < v)
            {
                const std::vector<unsigned>& states = it->first;
                l = states[1];
                r = states[2];
                break;
            }
            v++;
        }
    }
    else
    {
        double u = rng->uniform01() * sample_probs_sum;
        
        for (it = sample_probs.begin(); it != sample_probs.end(); it++)
        {
            u -= it->second;
            if (u < 0.0)
            {
                const std::vector<unsigned>& states = it->first;
                l = states[1];
                r = states[2];
                break;
            }
        }
    }

    // make nodes for each daughter
    TopologyNode* left = new TopologyNode(1);
    left->setAge(t);
    root->addChild(left);
    left->setParent(root);
    left->setTimeInStates(std::vector<double>(num_states, 0.0));
    left->setNumberOfShiftEvents(0);
    lineages_in_state[l].push_back(1);
    nodes.push_back(left);

    TopologyNode* right = new TopologyNode(2);
    right->setAge(t);
    root->addChild(right);
    right->setParent(root);
    right->setTimeInStates(std::vector<double>(num_states, 0.0));
    right->setNumberOfShiftEvents(0);
    lineages_in_state[r].push_back(2);
    nodes.push_back(right);

    // simulate moving forward in time
    while (true) {

        // sum over all rates for all states (multiplied by num lineages in each state)
        double total_rate = 0;
        for (size_t i = 0; i < num_states; i++)
        {
            total_rate += total_rate_for_state[i] * lineages_in_state[i].size();
        }
        
        // draw the time to next event
        double dt = RbStatistics::Exponential::rv( total_rate, *rng );
        if (condition_on_num_tips == true)
        {
            t = t + dt;
        }
        else
        {
            t = t - dt;
        }

        if (t < 0 && condition_on_num_tips == false)
        {
            dt = dt - (0 - t);
            t = 0;
        }

        // extend all surviving branches to the new time
        size_t num_lineages = 0;
        for (size_t i = 0; i < num_states; i++)
        {
            for (size_t j = 0; j < lineages_in_state[i].size(); j++)
            {
                size_t idx = lineages_in_state[i][j];
                nodes[idx]->setAge(t);
                num_lineages++;
                std::vector<double> state_times = nodes[idx]->getTimeInStates();
                state_times[i] += dt;
                nodes[idx]->setTimeInStates(state_times);
            }
        }

        // stop and retry if we have too many surviving lineages
        if (num_lineages > max_num_lineages && condition_on_num_tips == false)
        {
            nodes.clear();
            delete tip_data;
            return false;
        }
        
        // stop and retry if we reached the max time
        if (t > max_time && condition_on_num_tips == true)
        {
            nodes.clear();
            delete tip_data;
            return false;
        }

        // stop if we reached the present when conditioning on root age
        if (t == 0 && condition_on_num_tips == false)
        {
            for (size_t i = 0; i < nodes.size(); i++)
            {
                if (nodes[i]->getAge() == t)
                {
                    std::stringstream ss;
                    ss << "sp" << i;
                    std::string name = ss.str();
                    nodes[i]->setName(name);
                }
            }

            // set CharacterData object for each tip state
            for (size_t i = 0; i < num_states; i++)
            {
                for (size_t j = 0; j < lineages_in_state[i].size(); j++)
                {
                    size_t this_node = lineages_in_state[i][j];
                    if (nodes[this_node]->isTip() == true)
                    {
                        DiscreteTaxonData<NaturalNumbersState> this_tip_data = DiscreteTaxonData<NaturalNumbersState>(nodes[this_node]->getName());
                        NaturalNumbersState state = NaturalNumbersState(i, num_states);
                        this_tip_data.addCharacter(state);
                        tip_data->addTaxonData(this_tip_data);
                    }
                }
                if (prune_extinct_lineages == false)
                {
                    for (size_t j = 0; j < extinct_lineages_in_state[i].size(); j++)
                    {
                        size_t this_node = extinct_lineages_in_state[i][j];
                        if (nodes[this_node]->isTip() == true)
                        {
                            DiscreteTaxonData<NaturalNumbersState> this_tip_data = DiscreteTaxonData<NaturalNumbersState>(nodes[this_node]->getName());
                            NaturalNumbersState state = NaturalNumbersState(i, num_states);
                            this_tip_data.addCharacter(state);
                            tip_data->addTaxonData(this_tip_data);
                        }
                    }
                }
            }
            break;
        }

        // determine the state for the event that occurred
        size_t event_state = 0;
        double u = rng->uniform01() * total_rate;
        for (size_t i = 0; i < num_states; i++)
        {
            u -= total_rate_for_state[i] * lineages_in_state[i].size();
            if (u < 0)
            {
                event_state = i;
                break;
            }
        }

        // determine the type of event
        std::string event_type = "";
        u = rng->uniform01() * total_rate_for_state[event_state];
        while (true) {
            u = u - extinction_rates[event_state];
            if (u < 0)
            {
                event_type = "extinction";
                break;
            }
            u = u - total_speciation_rates[event_state];
            if (u < 0)
            {
                event_type = "speciation";
                break;
            }
            u = u - total_anagenetic_rates[event_state];
            if (u < 0)
            {
                event_type = "anagenetic";
                break;
            }
        }

        // determine which lineage gets the event
        size_t event_index = 0;
        u = rng->uniform01() * static_cast<double>(lineages_in_state[event_state].size());
        event_index = lineages_in_state[event_state][floor(u)];

        if (event_type == "extinction")
        {
            extinct_lineages_in_state[event_state].push_back(event_index);
            lineages_in_state[event_state].erase(std::remove(lineages_in_state[event_state].begin(), lineages_in_state[event_state].end(), event_index), lineages_in_state[event_state].end());
            std::stringstream ss;
            ss << "ex" << event_index;
            std::string name = ss.str();
            nodes[event_index]->setName(name);
        }
        
        if (event_type == "anagenetic")
        {
            // remove this lineage from the current state
            lineages_in_state[event_state].erase(std::remove(lineages_in_state[event_state].begin(), lineages_in_state[event_state].end(), event_index), lineages_in_state[event_state].end());

            // draw a new state
            size_t new_state = 0;
            u = rng->uniform01() * total_anagenetic_rates[event_state];
            for (size_t i = 0; i < this->num_states; i++)
            {
                if (i != event_state)
                {
                    u -= rate_matrix->getRate( event_state, i, 0, getEventRate( t ) );
                    if (u < 0.0)
                    {
                        new_state = i;
                        break;
                    }
                }
            }
            lineages_in_state[new_state].push_back(event_index);
            
            // increment the counter for the shift events
            nodes[event_index]->setNumberOfShiftEvents( nodes[event_index]->getNumberOfShiftEvents() + 1 );
        }
        
        if (event_type == "speciation")
        {

            // stop if we reached the right number of lineages when conditioning on num tips
            if (num_lineages == exact_num_lineages && condition_on_num_tips == true)
            {
                // TODO trim off uniformly distributed time since last speciation event
                
                for (size_t i = 0; i < nodes.size(); i++)
                {
                    if (nodes[i]->getAge() == t)
                    {
                        std::stringstream ss;
                        ss << "sp" << i;
                        std::string name = ss.str();
                        nodes[i]->setName(name);
                    }
                }
                
                // shift node times back so root starts at 0.0
                for (size_t i = 0; i < nodes.size(); i++)
                {
                    nodes[i]->setAge( t - nodes[i]->getAge() );
                }
                
                // set CharacterData object for each tip state
                for (size_t i = 0; i < num_states; i++)
                {
                    for (size_t j = 0; j < lineages_in_state[i].size(); j++)
                    {
                        size_t this_node = lineages_in_state[i][j];
                        if (nodes[this_node]->isTip() == true)
                        {
                            DiscreteTaxonData<NaturalNumbersState> this_tip_data = DiscreteTaxonData<NaturalNumbersState>(nodes[this_node]->getName());
                            NaturalNumbersState state = NaturalNumbersState(i, num_states);
                            this_tip_data.addCharacter(state);
                            tip_data->addTaxonData(this_tip_data);
                        }
                    }
                    if (prune_extinct_lineages == false)
                    {
                        for (size_t j = 0; j < extinct_lineages_in_state[i].size(); j++)
                        {
                            size_t this_node = extinct_lineages_in_state[i][j];
                            if (nodes[this_node]->isTip() == true)
                            {
                                DiscreteTaxonData<NaturalNumbersState> this_tip_data = DiscreteTaxonData<NaturalNumbersState>(nodes[this_node]->getName());
                                NaturalNumbersState state = NaturalNumbersState(i, num_states);
                                this_tip_data.addCharacter(state);
                                tip_data->addTaxonData(this_tip_data);
                            }
                        }
                    }
                }
                break;
            }

            // gather the probabilities for each type of cladogenetic event
            std::map<std::vector<unsigned>, double> sample_probs;
            double sample_probs_sum = 0.0;
            if ( use_cladogenetic_events == true )
            {
                // iterate over each cladogenetic event possible
                for (it = eventMap.begin(); it != eventMap.end(); it++)
                {
                    const std::vector<unsigned>& states = it->first;
                    double speciation_rate = it->second;
                    if (states[0] == event_state)
                    {
                        // we need to sample from the ancestor, left, and right states jointly,
                        // so keep track of the probability of each clado event
                        double prob = speciation_rate;
                        sample_probs[ states ] = prob;
                        sample_probs_sum += prob;
                    }
                }
            }
            else
            {
                std::vector<unsigned> states = boost::assign::list_of(event_state)(event_state)(event_state);
                sample_probs[ states ] = speciation_rates[event_state];
                sample_probs_sum += speciation_rates[event_state];
            }
            
            // sample left and right character states from probs
            size_t l = 0, r = 0;
            
            if (sample_probs_sum == 0)
            {
                size_t u = rng->uniform01() * sample_probs.size();
                size_t v = 0;
                for (it = sample_probs.begin(); it != sample_probs.end(); it++)
                {
                    if (u < v)
                    {
                        const std::vector<unsigned>& states = it->first;
                        l = states[1];
                        r = states[2];
                        break;
                    }
                    v++;
                }
            }
            else
            {
                double u = rng->uniform01() * sample_probs_sum;
                
                for (it = sample_probs.begin(); it != sample_probs.end(); it++)
                {
                    u -= it->second;
                    if (u < 0.0)
                    {
                        const std::vector<unsigned>& states = it->first;
                        l = states[1];
                        r = states[2];
                        break;
                    }
                }
            }
            
            // make nodes for each daughter
            size_t index = nodes.size();
            TopologyNode* left = new TopologyNode(index);
            left->setAge(t);
            nodes[event_index]->addChild(left);
            left->setParent(nodes[event_index]);
            left->setTimeInStates(std::vector<double>(num_states, 0.0));
            left->setNumberOfShiftEvents( 0 );
            lineages_in_state[l].push_back(index);
            nodes.push_back(left);

            index = nodes.size();
            TopologyNode* right = new TopologyNode(index);
            right->setAge(t);
            nodes[event_index]->addChild(right);
            right->setParent(nodes[event_index]);
            right->setTimeInStates(std::vector<double>(num_states, 0.0));
            right->setNumberOfShiftEvents( 0 );
            lineages_in_state[r].push_back(index);
            nodes.push_back(right);
           
            // remove the parent node from our vector of current lineages
            lineages_in_state[event_state].erase(std::remove(lineages_in_state[event_state].begin(), lineages_in_state[event_state].end(), event_index), lineages_in_state[event_state].end());
        }
    }
   
    // make a tree object
    Tree *sim_tree = new Tree();
    sim_tree->setRoot(root, true);
    sim_tree->setRooted(true);
        
    // stop and retry if we have too few surviving lineages
    size_t num_lineages = 0;
    for (size_t i = 0; i < num_states; i++)
    {
        num_lineages += lineages_in_state[i].size();
    }
  
    // prune extinct lineage if necessary
    if (prune_extinct_lineages == true)
    {
        for (size_t i = 0; i < num_states; i++)
        {
            for (size_t j = 0; j < extinct_lineages_in_state[i].size(); j++)
            {
                size_t this_node = extinct_lineages_in_state[i][j];
                if (nodes[this_node]->isTip() == true)
                {
                    sim_tree->dropTipNodeWithName( nodes[this_node]->getName() );
                }
            }
        }
    }
    
    if (sim_tree->getNumberOfTips() < min_num_lineages && condition_on_num_tips == false)
    {
        delete tip_data;
        nodes.clear();
        delete sim_tree;
        return false;
    }
    
    if ( (sim_tree->getNumberOfTips() < 2 || sim_tree->getRoot().getNumberOfChildren() != 2 || sim_tree->getRoot().getAge() != process_age->getValue()) && condition_on_tree == true)
    {
        delete tip_data;
        nodes.clear();
        delete sim_tree;
        return false;
    }
    
    // update character history vectors
    resizeVectors(sim_tree->getNumberOfNodes());
    simmap = "";
    for (size_t i = 0; i < sim_tree->getNumberOfNodes(); i++)
    {
        double branch_total_speciation = 0.0;
        double branch_total_extinction = 0.0;
        for (size_t j = 0; j < num_states; j++)
        {
            time_in_states[j] += sim_tree->getNodes()[i]->getTimeInStates()[j];
            branch_total_speciation += sim_tree->getNodes()[i]->getTimeInStates()[j] * total_speciation_rates[j];
            branch_total_extinction += sim_tree->getNodes()[i]->getTimeInStates()[j] * extinction_rates[j];
        }
        if (sim_tree->getNodes()[i]->getBranchLength() > 0)
        {
            average_speciation[i] = branch_total_speciation/sim_tree->getNodes()[i]->getBranchLength();
            average_extinction[i] = branch_total_extinction/sim_tree->getNodes()[i]->getBranchLength();
            num_shift_events[i]   = sim_tree->getNodes()[i]->getNumberOfShiftEvents();
        }
    }
    
    // set the simulated values
    value->getTreeChangeEventHandler().removeListener( this );
    static_cast<TreeDiscreteCharacterData *>(this->value)->setTree( *sim_tree );
    delete sim_tree;
    value->getTreeChangeEventHandler().addListener( this );
    static_cast<TreeDiscreteCharacterData*>(this->value)->setCharacterData(tip_data);
    static_cast<TreeDiscreteCharacterData*>(this->value)->setTimeInStates(time_in_states);
    return true;
    
}


/**
 * Sorts times to run from present to past (0->inf) and orders par to match this.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::sortNonGlobalTimesAndParameters(std::vector<RbVector<double> >& par, std::vector<double> &times) const
{
    std::vector<double> times_sorted_ascending = times;
    std::vector<double> times_sorted_descending = times;

    sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );
    sort(times_sorted_descending.rbegin(), times_sorted_descending.rend() );

    // We want times in ascending order, so if they already are we're done here
    if ( times != times_sorted_ascending )
    {
        // If times are sorted in descending order, we just flip the parameter and time vectors
        if ( times == times_sorted_ascending )
        {
            std::reverse(times.begin(),times.end());
            std::reverse(par.begin(),par.end());
        }
        else
        {
            // Pair up the times and the parameter values so we can sort them together
            std::vector<std::pair<double, RbVector<double> > > times_par;
            for (size_t i=0; i<times.size(); ++i)
            {
                times_par.push_back(std::make_pair(times[i],par[i]));
            }

            std::sort(times_par.begin(),times_par.end());

            // Replace times with sorted times
            for (size_t i=0; i<times.size(); ++i)
            {
                times[i] = times_par[i].first;
                par[i] = times_par[i].second;
            }
        }
    }

    if ( times[0] < DBL_EPSILON )
    {
        throw RbException("User-specified interval times cannot include time = 0");
    }

}


/**
 * Sorts times to run from present to past (0->inf) and orders par to match this.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::sortNonGlobalTimesAndParameters(std::vector<double> &par, std::vector<double> &times) const
{
    std::vector<double> times_sorted_ascending = times;
    std::vector<double> times_sorted_descending = times;

    sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );
    sort(times_sorted_descending.rbegin(), times_sorted_descending.rend() );

    // We want times in ascending order, so if they already are we're done here
    if ( times != times_sorted_ascending )
    {
        // If times are sorted in descending order, we just flip the parameter and time vectors
        if ( times == times_sorted_ascending )
        {
            std::reverse(times.begin(),times.end());
            std::reverse(par.begin(),par.end());
        }
        else
        {
            // Pair up the times and the parameter values so we can sort them together
            std::vector<std::pair<double,double> > times_par;
            for (size_t i=0; i<times.size(); ++i)
            {
                times_par.push_back(std::make_pair(times[i],par[i]));
            }

            std::sort(times_par.begin(),times_par.end());

            // Replace times with sorted times
            for (size_t i=0; i<times.size(); ++i)
            {
                times[i] = times_par[i].first;
                par[i] = times_par[i].second;
            }
        }
    }

    if ( times[0] < DBL_EPSILON )
    {
        throw RbException("User-specified interval times cannot include time = 0");
    }

}



/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == process_age )
    {
        process_age = static_cast<const TypedDagNode<double>* >( newP );
    }
    if ( oldP == mu_const )
    {
        mu_const = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == mu_var )
    {
        mu_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
    if ( oldP == lambda_const )
    {
        lambda_const = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == lambda_var )
    {
        lambda_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
    if ( oldP == phi_const )
    {
        phi_const = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == phi_var )
    {
        phi_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
    if ( oldP == Q_const )
    {
        Q_const = static_cast<const TypedDagNode<RateGenerator>* >( newP );
    }
    if ( oldP == Q_var )
    {
        Q_var = static_cast<const TypedDagNode<RbVector<RateGenerator> >* >( newP );
    }
    if ( oldP == eta_const )
    {
        eta_const = static_cast<const TypedDagNode<double>* >( newP );
    }
    if ( oldP == eta_var )
    {
        eta_var = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == epoch_times_lambda )
    {
        epoch_times_lambda = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == epoch_times_mu )
    {
        epoch_times_mu = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == epoch_times_phi )
    {
        epoch_times_phi = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == epoch_times_gamma )
    {
        epoch_times_gamma = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == epoch_times_eta )
    {
        epoch_times_eta = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == epoch_times_Q )
    {
        epoch_times_Q = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == pi )
    {
        pi = static_cast<const TypedDagNode<Simplex>* >( newP );
    }
    if ( oldP == rho )
    {
        rho = static_cast<const TypedDagNode<double>* >( newP );
    }
    if ( oldP == rho_per_state )
    {
        rho_per_state = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if ( oldP == cladogenesis_matrix )
    {
        cladogenesis_matrix = static_cast<const TypedDagNode<CladogeneticSpeciationRateMatrix>* >( newP );
    }
    
}



/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    
    if ( affecter == process_age )
    {
        if ( use_origin == false)
        {
            value->getRoot().setAge( process_age->getValue() );
        }

        if ( dag_node != NULL )
        {
            dag_node->touchAffected();
        }
    }
    
    if ( affecter != this->dag_node )
    {
        
        for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
        {
            (*it) = true;
        }
        
        // flip the active likelihood pointers
        for (size_t index = 0; index < changed_nodes.size(); ++index)
        {
            if ( changed_nodes[index] == false )
            {
                active_likelihood[index] = (active_likelihood[index] == 0 ? 1 : 0);
                changed_nodes[index] = true;
            }
        }
    }
    
}


/**
 * Wrapper function for the ODE time stepper function.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::numericallyIntegrateProcess(std::vector< double > &likelihoods, double begin_age, double end_age, bool backward_time, bool extinction_only) const
{
    
    size_t index_epoch_begin = 0;
    size_t index_epoch_end = 0;
    
    if ( backward_time == true )
    {
        index_epoch_begin = computeEpochIndex( begin_age );
        index_epoch_end = computeEpochIndex( end_age );
    }

    double current_begin_age = begin_age;
        
    for ( size_t index_epoch=index_epoch_begin; index_epoch<=index_epoch_end; ++index_epoch )
    {
        
        double epoch_end = computeEpochEnd( index_epoch );
        double current_end_age = (end_age < epoch_end ? end_age : epoch_end );
        
        const RbVector<double> &extinction_rates = computeExtinctionRateAtTime(current_end_age);
        const RateGenerator &rg = getEventRateMatrix( current_begin_age );
        SSE_ODE ode = SSE_ODE(extinction_rates, &rg, getEventRate(current_end_age), backward_time, extinction_only);
        if ( use_cladogenetic_events == true )
        {
            cladogenesis_matrix->getValue(); // we must call getValue() to update the speciation and extinction rates in the event map
        
            // get cladogenesis event map (sparse speciation rate matrix)
            std::map<std::vector<unsigned>, double> event_map = cladogenesis_matrix->getValue().getEventMap();
        
            ode.setEventMap( event_map );
        }
        else
        {
            const RbVector<double> &speciation_rates = computeSpeciationRateAtTime(current_end_age);
            ode.setSpeciationRate( speciation_rates );
        }
    
        if ( phi_var != NULL || phi_const != NULL )
        {
            const RbVector<double> &fossilization_rates = computeFossilizationRateAtTime(current_end_age);
            ode.setSerialSamplingRate( fossilization_rates );
        }
    
        typedef boost::numeric::odeint::runge_kutta_dopri5< std::vector< double > > stepper_type;

        boost::numeric::odeint::integrate_adaptive( make_controlled( 1E-9, 1E-9, stepper_type() ) , ode , likelihoods , current_begin_age , current_end_age , dt );
//        boost::numeric::odeint::integrate_adaptive( stepper_type(), ode , likelihoods , current_begin_age , current_end_age , dt );
    
        // catch negative extinction probabilities that can result from
        // rounding errors in the ODE stepper
        for (size_t i = 0; i < 2 * num_states; ++i)
        {
        
            // Sebastian: The likelihoods here are probability densities (not log-transformed).
            // These are densities because they are multiplied by the probability density of the speciation event happening.
            likelihoods[i] = ( likelihoods[i] < 0.0 ? 0.0 : likelihoods[i] );
        }
        
        // catch too large extinction probabilities that can result from
        // rounding errors in the ODE stepper
        // for safety we set all likelihoods to nan if rounding errors happened
        bool rounding_error = false;
        for (size_t i = 0; i < num_states; ++i)
        {
            
            // Sebastian: The extinction probabilities here are probabilities (not log-transformed).
            // So they must be between 0 and 1.
            rounding_error |= ( likelihoods[i] > 1.0 );
            
        }
        
        if ( rounding_error == true )
        {
            for (size_t i = 0; i < (2*num_states); ++i)
            {
                
                // invalidate likelihoods
                likelihoods[i] = RbConstants::Double::nan;
                
            }
        }
        
        if ( index_epoch < index_epoch_end )
        {
            const RbVector<double>& surv_probs = computeSurvivalProbabilitiesAtTime( current_end_age );
            for (size_t i = 0; i < num_states; ++i)
            {
                
                // Sebastian: The extinction probabilities here are probabilities (not log-transformed).
                // So they must be between 0 and 1.
                likelihoods[i] = (1.0-surv_probs[i]) + surv_probs[i] * likelihoods[i];
                likelihoods[i+num_states] = surv_probs[i] * likelihoods[i+num_states];

            }
        }
        
        current_begin_age = current_end_age;
        
    }

    
}


/**
 * Resize various vectors depending on the current number of nodes.
 */
void EpisodicStateDependentSpeciationExtinctionFossilizationProcess::resizeVectors(size_t num_nodes)
{
    active_likelihood = std::vector<bool>(num_nodes, false);
    changed_nodes = std::vector<bool>(num_nodes, false);
    dirty_nodes = std::vector<bool>(num_nodes, true);
    node_partial_likelihoods = std::vector<std::vector<std::vector<double> > >(num_nodes, std::vector<std::vector<double> >(2,std::vector<double>(2*num_states,0)));
    scaling_factors = std::vector<std::vector<double> >(num_nodes, std::vector<double>(2,0.0) );
    average_speciation = std::vector<double>(num_nodes, 0.0);
    average_extinction = std::vector<double>(num_nodes, 0.0);
    num_shift_events = std::vector<std::int64_t>(num_nodes, 0.0);
    time_in_states = std::vector<double>(num_states, 0.0);
}
