#include <boost/assign/list_of.hpp>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstddef>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "DistributionLognormal.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "BDS_ODE.h"
#include "DistributionExponential.h"
#include "HomologousDiscreteCharacterData.h"
#include "FastBirthDeathShiftProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "AbstractDiscreteTaxonData.h"
#include "AbstractTaxonData.h"
#include "DiscreteCharacterState.h"
#include "DiscreteTaxonData.h"
#include "NaturalNumbersState.h"
#include "RbException.h"
#include "RbSettings.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "Simplex.h"
#include "Taxon.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TreeDiscreteCharacterData.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "boost/numeric/odeint.hpp" // IWYU pragma: keep


namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }


using namespace RevBayesCore;


/**
 * Constructor.
 *
 * The constructor connects the parameters of the birth-death process (DAG structure)
 * and initializes the probability density by computing the combinatorial constant of the tree structure.
 */
FastBirthDeathShiftProcess::FastBirthDeathShiftProcess(const TypedDagNode<double> *age,
                                                       const TypedDagNode<double> *sp,
                                                       const TypedDagNode<double> *ext,
                                                       const TypedDagNode<double> *sp_sd,
                                                       const TypedDagNode<double> *ext_sd,
                                                       const TypedDagNode<double> *r_sp,
                                                       const TypedDagNode<double> *r_ext,
                                                       size_t num_classes, 
                                                       const std::string &cdt,
                                                       bool uo,
                                                       size_t min_num_lineages,
                                                       size_t max_num_lineages,
                                                       size_t exact_num_lineages,
                                                       double max_t,
                                                       bool prune,
                                                       bool condition_on_tip_states,
                                                       bool condition_on_num_tips,
                                                       bool condition_on_tree) : TypedDistribution<Tree>( new TreeDiscreteCharacterData() ),
    condition( cdt ),
    active_likelihood( std::vector<bool>(5, 0) ),
    changed_nodes( std::vector<bool>(5, false) ),
    dirty_nodes( std::vector<bool>(5, true) ),
    node_partial_likelihoods( std::vector<std::vector<std::vector<double> > >(5, std::vector<std::vector<double> >(2,std::vector<double>(2*num_classes*num_classes,0))) ),
    num_states( num_classes*num_classes ),
    scaling_factors( std::vector<std::vector<double> >(5, std::vector<double>(2,0.0) ) ),
    use_origin( uo ),
    sample_character_history( false ),
    average_speciation( std::vector<double>(5, 0.0) ),
    average_extinction( std::vector<double>(5, 0.0) ),
    num_shift_events( std::vector<long>(5, 0.0) ),
    time_in_states( std::vector<double>(num_classes*num_classes, 0.0) ),    
    simmap( "" ),
    process_age( age ),
    num_rate_classes ( num_classes ),
    speciation_scale( sp ),
    extinction_scale( ext ),
    speciation_sd( sp_sd ),
    extinction_sd( ext_sd ),
    alpha( r_sp ),
    beta( r_ext ),
    lambda( std::vector<double>(num_classes * num_classes, 0.0) ),
    mu( std::vector<double>(num_classes * num_classes, 0.0) ),
    rho( new ConstantNode<double>("", new double(1.0)) ),
    min_num_lineages( min_num_lineages ),
    max_num_lineages( max_num_lineages ),
    exact_num_lineages( exact_num_lineages ),
    max_time( max_t ),
    prune_extinct_lineages( prune ),
    condition_on_tip_states( condition_on_tip_states ),
    condition_on_num_tips( condition_on_num_tips ),
    condition_on_tree( condition_on_tree )
{
    addParameter( sp );
    addParameter( ext );
    addParameter( sp_sd );
    addParameter( ext_sd );
    addParameter( alpha );
    addParameter( beta );
    addParameter( rho );
    addParameter( process_age );


    // set the new Q matrix
    boost::numeric::ublas::matrix<double> Qm(num_states, num_states);
    Qmatrix = Qm;
    updateQmatrix();
    update_rates();

    if ( min_num_lineages > max_num_lineages )
    {
        throw RbException("minNumLineages cannot be greater than maxNumLineages.");
    }
    
    // set the length of the time slices used by the ODE for numerical integration
    dt = process_age->getValue() / 500 * 10.0;

    value->getTreeChangeEventHandler().addListener( this );

}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
FastBirthDeathShiftProcess* FastBirthDeathShiftProcess::clone( void ) const
{
    FastBirthDeathShiftProcess* tmp = new FastBirthDeathShiftProcess( *this );
    tmp->getValue().getTreeChangeEventHandler().addListener(tmp);
    return tmp;
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
FastBirthDeathShiftProcess::~FastBirthDeathShiftProcess( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!

    // remove myself from the tree listeners
    value->getTreeChangeEventHandler().removeListener( this );
}


std::vector<double> FastBirthDeathShiftProcess::calculateTotalSpeciationRatePerState( void ) const
{
    std::vector<double> total_rates = std::vector<double>(num_states, 0);
    std::vector<double> speciation_rates;
    std::map<std::vector<unsigned>, double>::iterator it;

    for (size_t i = 0; i < num_states; i++)
    {
        total_rates[i] += lambda[i];
    }
    return total_rates;
}


std::vector<double> FastBirthDeathShiftProcess::calculateTotalAnageneticRatePerState( void ) const
{
    std::vector<double> total_rates = std::vector<double>(num_states, 0);
    for (size_t i = 0; i < num_states; i++)
    {
        for (size_t j = 0; j < num_states; j++)
        {
            if (i != j)
            {
                total_rates[i] += Qmatrix(i,j);
            }
        }
    }

    return total_rates;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FastBirthDeathShiftProcess::computeLnProbability( void )
{
    
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


void FastBirthDeathShiftProcess::computeNodeProbability(const RevBayesCore::TopologyNode &node, size_t node_index) const
{
   
    // check for recomputation
    if ( dirty_nodes[node_index] == true || sample_character_history == true )
    {
        // mark as computed
        dirty_nodes[node_index] = false;
        
        std::vector<double> &node_likelihood = node_partial_likelihoods[node_index][active_likelihood[node_index]];

        if ( node.isTip() ) // this is a tip node
        {
            TreeDiscreteCharacterData* tree = static_cast<TreeDiscreteCharacterData*>( this->value );

            std::vector<double> sampling;
            std::vector<double> extinction;
            if ( rho != NULL )
            {
                sampling   = std::vector<double>(num_states, rho->getValue());
                extinction = std::vector<double>(num_states, 1.0 - rho->getValue());
            }
            else
            {
                throw RbException("A species sampling probability needs to be set.");
            }

            if ( node.isFossil() )
            {
                throw(RbException("The model does not support fossil tips."));
            }
            
            for (size_t j = 0; j < num_states; ++j)
            {
                node_likelihood[j] = extinction[j];
                node_likelihood[num_states+j] = sampling[j];
            }

            scaling_factors[node_index][active_likelihood[node_index]] = 0.0;
        }
        else // this is an internal node
        {
            const TopologyNode          &left           = node.getChild(0);
            size_t                      left_index      = left.getIndex();
            computeNodeProbability( left, left_index );

            const TopologyNode          &right          = node.getChild(1);
            size_t                      right_index     = right.getIndex();
            computeNodeProbability( right, right_index );
            
            // get the likelihoods of descendant nodes
            const std::vector<double> &left_likelihoods  = node_partial_likelihoods[left_index][active_likelihood[left_index]];
            const std::vector<double> &right_likelihoods = node_partial_likelihoods[right_index][active_likelihood[right_index]];

            // merge descendant likelihoods
            for (size_t i=0; i<num_states; ++i)
            {
                node_likelihood[i] = left_likelihoods[i];

                node_likelihood[num_states + i] = left_likelihoods[num_states + i] * right_likelihoods[num_states + i];
                node_likelihood[num_states + i] *= lambda[i];
            }

            // propagate the scaling factor (not re-scaling yet)
            scaling_factors[node_index][active_likelihood[node_index]] = scaling_factors[left_index][active_likelihood[left_index]] + scaling_factors[right_index][active_likelihood[right_index]];
            
        }
      
        // find the time span for the ODE
        double begin_age;
        double end_age;

        if (node.isRoot()){
            begin_age = getRootAge();

            if (use_origin == true){
                end_age = getOriginAge();
            }else{
                end_age = begin_age;
            }
        }else{
            begin_age = node.getAge();
            end_age = node.getParent().getAge();
        }
        
        // numerically integrate over the branch 
        numericallyIntegrateProcess(node_likelihood, begin_age, end_age, true, false);
        
        // rescale the probability densities at the "end" of the branch
        if ( RbSettings::userSettings().getUseScaling() == true ) 
        {
            double max;
            max = rescaleProbabilities(node_likelihood);

            scaling_factors[node_index][active_likelihood[node_index]] += log(max);
        }
    }

}

double FastBirthDeathShiftProcess::rescaleProbabilities(std::vector<double>& probabilities) const
{
    double max = 0.0;
    for (size_t i=0; i<num_states; ++i)
    {
        if ( probabilities[num_states+i] > max )
        {
            max = probabilities[num_states+i];
        }
    }
    
    for (size_t i=0; i<num_states; ++i)
    {
        probabilities[num_states+i] /= max;
    }

    return max;
}


double FastBirthDeathShiftProcess::computeRootLikelihood( void ) const
{
    // get the likelihoods of descendant nodes
    const TopologyNode     &root            = value->getRoot();
    size_t                  node_index      = root.getIndex();
  
    computeNodeProbability( root, node_index);
    std::vector<double> &node_likelihood = node_partial_likelihoods[node_index][active_likelihood[node_index]];

    // sum the root likelihoods
    const RbVector<double> &freqs = getRootFrequencies();
    double prob = 0.0;

    for (size_t i = 0; i < num_states; ++i)
    {
        prob += freqs[i] * node_likelihood[num_states + i];
    }
    
    return log(prob) + scaling_factors[node_index][active_likelihood[node_index]];
}


void FastBirthDeathShiftProcess::drawStochasticCharacterMap(std::vector<std::string>& character_histories, bool set_amb_char_data)
{
    // first populate partial likelihood vectors along all the branches
    sample_character_history = true;
    computeLnProbability();

    size_t attempts = 0;
    bool success = false;
    while (!success)
    {
        if (attempts == 100000)
        {
            throw RbException("After 100000 attempts a character history could not be sampled with a non-zero probability. Try increasing nTimeSlices.");
        }

        for (size_t i = 0; i < num_states; i++) 
        {
            time_in_states[i] = 0.0;
        }

        // now begin the root-to-tip pass, drawing ancestral states for each time slice conditional on the start states

        // get the likelihoods of descendant nodes
        const TopologyNode          &root               = value->getRoot();
        size_t                       node_index         = root.getIndex();
        const TopologyNode          &left               = root.getChild(0);
        size_t                       left_index         = left.getIndex();
        const std::vector< double > &left_likelihoods   = node_partial_likelihoods[left_index][active_likelihood[left_index]];
        const TopologyNode          &right              = root.getChild(1);
        size_t                       right_index        = right.getIndex();
        const std::vector< double > &right_likelihoods  = node_partial_likelihoods[right_index][active_likelihood[right_index]];
        
        // get root frequencies
        const RbVector<double> &freqs = getRootFrequencies();
        
        std::map<std::vector<unsigned>, double> sample_probs;
        double sample_probs_sum = 0.0;
        std::map<std::vector<unsigned>, double>::iterator it;
        
        // calculate probabilities for each state
        for (size_t i = 0; i < num_states; i++)
        {
            double likelihood = left_likelihoods[num_states + i] * right_likelihoods[num_states + i] * lambda[i];
            std::vector<unsigned> states = boost::assign::list_of(i)(i)(i);
            sample_probs[ states ] = likelihood * freqs[i];
            sample_probs_sum += likelihood * freqs[i];
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
            bool success_l = recursivelyDrawStochasticCharacterMap(left, l, character_histories, set_amb_char_data);
            bool success_r = recursivelyDrawStochasticCharacterMap(right, r, character_histories, set_amb_char_data);
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


bool FastBirthDeathShiftProcess::recursivelyDrawStochasticCharacterMap(const TopologyNode &node, size_t start_state, std::vector<std::string>& character_histories, bool set_amb_char_data)
{
    size_t node_index = node.getIndex();
    //std::vector<double> speciation_rates = calculateTotalSpeciationRatePerState();
    //std::vector<double> extinction_rates = mu->getValue();
    
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
        total_speciation_rate += lambda[current_state];
        total_extinction_rate +=     mu[current_state];
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
        total_speciation_rate += lambda[new_state];
        total_extinction_rate +=     mu[new_state];
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
        for (size_t i = transition_times.size(); i > 0; i--)
        {
            simmap_string = simmap_string + StringUtilities::toString(transition_states[i - 1]) + "," + StringUtilities::toString(transition_times[i - 1]);
            if (i != 1)
            {
                simmap_string = simmap_string + ":";
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
        for (size_t i = 0; i < num_states; i++)
        {
            double prob = left_likelihoods[num_states + i] * right_likelihoods[num_states + i] * lambda[i];
            prob *= branch_conditional_probs[num_states + i];
            std::vector<unsigned> states = boost::assign::list_of(i)(i)(i);
            sample_probs[ states ] = prob;
            sample_probs_sum += prob;
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
        total_speciation_rate += lambda[a];
        total_extinction_rate += mu[a];
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
        bool success_l = recursivelyDrawStochasticCharacterMap(left, l, character_histories, set_amb_char_data);
        bool success_r = recursivelyDrawStochasticCharacterMap(right, r, character_histories, set_amb_char_data);
        return success_l && success_r;
    }
    return true;
}




void FastBirthDeathShiftProcess::fireTreeChangeEvent( const RevBayesCore::TopologyNode &n, const unsigned& m )
{
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
}


const RevBayesCore::AbstractHomologousDiscreteCharacterData& FastBirthDeathShiftProcess::getCharacterData() const
{
    return static_cast<TreeDiscreteCharacterData*>(this->value)->getCharacterData();
}


void FastBirthDeathShiftProcess::drawJointConditionalAncestralStates(std::vector<size_t>& startStates, std::vector<size_t>& endStates)
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
    
    // get root frequencies
    const RbVector<double> &freqs = getRootFrequencies();
    
    std::map<std::vector<unsigned>, double> sample_probs;
    double sample_probs_sum = 0.0;
    std::map<std::vector<unsigned>, double>::iterator it;
    
    // calculate probabilities for each state
    for (size_t i = 0; i < num_states; i++) {
            double likelihood = left_likelihoods[num_states + i] * right_likelihoods[num_states + i] * lambda[i];
            std::vector<unsigned> states = boost::assign::list_of(i)(i)(i);
            sample_probs[ states ] = likelihood * freqs[i];
            sample_probs_sum += likelihood * freqs[i];
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


void FastBirthDeathShiftProcess::recursivelyDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<size_t>& startStates, std::vector<size_t>& endStates)
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
        for (size_t i = 0; i < num_states; i++){
                double prob = left_likelihoods[num_states + i] * right_likelihoods[num_states + i] * lambda[i];
                prob *= branch_conditional_probs[num_states + i];
                std::vector<unsigned> states = boost::assign::list_of(i)(i)(i);
                sample_probs[ states ] = prob;
                sample_probs_sum += prob;
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


void FastBirthDeathShiftProcess::recursivelyFlagNodeDirty( const RevBayesCore::TopologyNode &n ) {

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


RevLanguage::RevPtr<RevLanguage::RevVariable> FastBirthDeathShiftProcess::executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found)
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
    return TypedDistribution<Tree>::executeProcedure( name, args, found );
}


void FastBirthDeathShiftProcess::executeMethod(const std::string &name, const std::vector<const DagNode *> &args, RbVector<long> &rv) const
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


void FastBirthDeathShiftProcess::executeMethod(const std::string &name, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
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
 * Get the affected nodes by a change of this node.
 * If the root age has changed than we need to call get affected again.
 */
void FastBirthDeathShiftProcess::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    if ( affecter == process_age )
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
}



double FastBirthDeathShiftProcess::getOriginAge( void ) const
{

    return process_age->getValue();
}


std::vector<double> FastBirthDeathShiftProcess::getAverageExtinctionRatePerBranch( void ) const
{
    return average_extinction;
}


std::vector<double> FastBirthDeathShiftProcess::getAverageSpeciationRatePerBranch( void ) const
{
    return average_speciation;
}


std::vector<long> FastBirthDeathShiftProcess::getNumberOfShiftEventsPerBranch( void ) const
{
    return num_shift_events;
}


std::vector<double> FastBirthDeathShiftProcess::getTimeInStates( void ) const
{
    return time_in_states;
}


/**
 * By default, the root age is assumed to be equal to the origin time.
 * This should be overridden if a distinct root age is needed
 */
double FastBirthDeathShiftProcess::getRootAge( void ) const
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
std::vector<double> FastBirthDeathShiftProcess::getRootFrequencies(void) const
{
    return std::vector<double>(num_states, 1.0/num_states);
}


/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void FastBirthDeathShiftProcess::keepSpecialization(const DagNode *affecter)
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


double FastBirthDeathShiftProcess::lnProbTreeShape(void) const
{
    // the birth death divergence times density is derived for a (ranked) unlabeled oriented tree
    // so we convert to a (ranked) labeled non-oriented tree probability by multiplying by 2^{n+m-1} / n!
    // where n is the number of extant tips, m is the number of extinct tips

    int num_taxa = (int)value->getNumberOfTips();
    int num_extinct = (int)value->getNumberOfExtinctTips();
    int num_sa = (int)value->getNumberOfSampledAncestors();

    return (num_taxa - num_sa - 1) * RbConstants::LN2 - RbMath::lnFactorial(num_taxa - num_extinct);
}


std::vector<double> FastBirthDeathShiftProcess::pExtinction(double start, double end) const
{
    
    
    std::vector<double> sampling_probability;
    if ( rho != NULL )
    {
        sampling_probability = std::vector<double>(num_states, rho->getValue());
    }
    else
    {
        throw RbException("A global sampling fraction needs to be set.");
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


double FastBirthDeathShiftProcess::pSurvival(double start, double end) const
{

    // delegate to specific function that manages survival between origin and root
    return pSurvival(start, end, use_origin == false);
}


double FastBirthDeathShiftProcess::pSurvival(double start, double end, bool speciation) const
{

    std::vector< double > initial_state = pExtinction(start,end);
    std::vector<double> speciation_rates = calculateTotalSpeciationRatePerState();

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
void FastBirthDeathShiftProcess::redrawValue( void )
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
void FastBirthDeathShiftProcess::restoreSpecialization(const DagNode *affecter)
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

void FastBirthDeathShiftProcess::setSampleCharacterHistory(bool sample_history)
{
    sample_character_history = sample_history;
}


void FastBirthDeathShiftProcess::setSamplingFraction(const TypedDagNode<double> *f)
{
    // remove the old parameter first
    this->removeParameter( rho );

    // set the value
    rho = f;
    
    // add the new parameter
    this->addParameter( rho );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}

/**
 * Set the current value.
 */
void FastBirthDeathShiftProcess::setValue(Tree *v, bool f )
{
    if (v->isBinary() == false)
    {
        throw RbException("The character-dependent birth death process is only implemented for binary trees.");
    }

    value->getTreeChangeEventHandler().removeListener( this );

    // delegate to super class
    //    TypedDistribution<Tree>::setValue(v, f);
    static_cast<TreeDiscreteCharacterData *>(this->value)->setTree( *v );

    resizeVectors(v->getNumberOfNodes());
    
    // clear memory
    delete v;
    
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


bool FastBirthDeathShiftProcess::simulateTreeConditionedOnTips( size_t attempts )
{

    if ( prune_extinct_lineages == false )
    {
        throw RbException("Simulations conditioned on the tip states are currently implemented only when pruneExtinctLineages is set to true.");
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
    std::vector<double> total_speciation_rates = calculateTotalSpeciationRatePerState();
    std::vector<double> total_anagenetic_rates = calculateTotalAnageneticRatePerState();
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
                mu[i] -= t/10; // TODO make this a user option
                if (mu[i] < 0)
                {
                    mu[i] = 0.0;
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
                r[i] += mu[i] + total_anagenetic_rates[i];
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
                prob_extinction[i] = mu[i] * (lineages_in_state[i].size() + 1) * exp(-1 * dt * total_rate_ext);

                for (size_t j = 0; j < num_states; ++j)
                {
                    if (i != j && lineages_in_state[j].size() > 0) 
                    {
                        prob_transition[i][j] = Qmatrix(i,j) * (lineages_in_state[i].size() + 1) * exp(-1 * dt * total_rate_ana[j]);
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
            branch_total_extinction += sim_tree->getNodes()[i]->getTimeInStates()[j] * mu[j];
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
bool FastBirthDeathShiftProcess::simulateTree( size_t attempts )
{

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
    std::vector<double> total_speciation_rates = calculateTotalSpeciationRatePerState();
    std::vector<double> total_anagenetic_rates = calculateTotalAnageneticRatePerState();
    std::vector<double> total_rate_for_state = std::vector<double>(num_states, 0);
    for (size_t i = 0; i < num_states; i++)
    {
        total_rate_for_state[i] = mu[i] + total_speciation_rates[i] + total_anagenetic_rates[i];
    }

    // get the speciation rates, extinction rates, and Q matrix
    std::map<std::vector<unsigned>, double>::iterator it;


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

    // get root frequencies
    const RbVector<double> &root_freqs = getRootFrequencies();
    
    std::map<std::vector<unsigned>, double> sample_probs;
    double sample_probs_sum = 0.0;
    
    // calculate probabilities for each state
    for (size_t i = 0; i < num_states; i++){
        std::vector<unsigned> states = boost::assign::list_of(i)(i)(i);
        sample_probs[ states ] = root_freqs[i] * lambda[i];
        sample_probs_sum += root_freqs[i] * lambda[i];
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
            u = u - mu[event_state];
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
                    u -= Qmatrix(i,0);
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

            std::vector<unsigned> states = boost::assign::list_of(event_state)(event_state)(event_state);
            sample_probs[ states ] = lambda[event_state];
            sample_probs_sum += lambda[event_state];
            
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
            branch_total_extinction += sim_tree->getNodes()[i]->getTimeInStates()[j] * mu[j];
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
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void FastBirthDeathShiftProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == process_age )
    {
        process_age = static_cast<const TypedDagNode<double>* >( newP );
    }
    if ( oldP == speciation_scale )
    {
        speciation_scale = static_cast<const TypedDagNode<double>* >( newP );
        update_rates();
    }
    if ( oldP == extinction_scale )
    {
        extinction_scale = static_cast<const TypedDagNode<double>* >( newP );
        update_rates();
    }
    if ( oldP == speciation_sd )
    {
        speciation_sd = static_cast<const TypedDagNode<double>* >( newP );
        update_rates();
    }
    if ( oldP == extinction_sd )
    {
        extinction_sd = static_cast<const TypedDagNode<double>* >( newP );
        update_rates();
    }
    if ( oldP == alpha )
    {
        alpha = static_cast<const TypedDagNode<double>* >( newP );
        updateQmatrix();
    }
    if ( oldP == beta )
    {
        beta = static_cast<const TypedDagNode<double>* >( newP );
        updateQmatrix();
    }
    if ( oldP == rho )
    {
        rho = static_cast<const TypedDagNode<double>* >( newP );
    }
}



/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void FastBirthDeathShiftProcess::touchSpecialization(const DagNode *affecter, bool touchAll)
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
void FastBirthDeathShiftProcess::numericallyIntegrateProcess(std::vector< double > &likelihoods, double begin_age, double end_age, bool backward_time, bool extinction_only) const
{
    const std::vector<double> &speciation_rates = lambda;
    const std::vector<double> &extinction_rates = mu;
    const double &alpha_ref = alpha->getValue();
    const double &beta_ref = beta->getValue();


    // construct the Q matrix
    //updateQmatrix();

    boost::numeric::ublas::matrix<double> &Qref = Qmatrix;
    //const size_t num_classes = sqrt(speciation_rates.size());


    //BDS_ODE ode = BDS_ODE(speciation_rates, extinction_rates, &getEventRateMatrix());
    BDS_ODE ode = BDS_ODE(speciation_rates, extinction_rates, num_rate_classes, alpha_ref, beta_ref);
   
    typedef boost::numeric::odeint::runge_kutta_dopri5< std::vector< double > > stepper_type;

    boost::numeric::odeint::integrate_adaptive( make_controlled( 1E-6, 1E-3, stepper_type() ) , ode , likelihoods , begin_age , end_age , dt );
    //boost::numeric::odeint::integrate_adaptive( stepper_type(), ode , likelihoods , begin_age , end_age , dt );
    
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
    
}

void FastBirthDeathShiftProcess::update_rates(){
    size_t num_categories = num_rate_classes * num_rate_classes;

    double lambda_scale = speciation_scale->getValue();
    double lambda_sd    = speciation_sd   ->getValue();

    double mu_scale = extinction_scale->getValue();
    double mu_sd    = extinction_sd   ->getValue();

    std::vector<double> v_lambda;
    std::vector<double> v_mu;

    for (size_t i=0; i < num_rate_classes; ++i){
        double p = (i+0.5)/num_rate_classes;

        double lambda_quantile = RbStatistics::Lognormal::quantile(log(lambda_scale), lambda_sd, p);
        double mu_quantile     = RbStatistics::Lognormal::quantile(log(mu_scale),     mu_sd,     p);

        v_lambda.push_back(lambda_quantile);
        v_mu.push_back(mu_quantile);
    }

    size_t q = 0;
    for (size_t i = 0; i < num_rate_classes; i++){
        for (size_t j = 0; j < num_rate_classes; j++){
            lambda[q] = v_lambda[i]; 
            mu[q] = v_mu[j]; 
            q++;
        }
    }
}


void FastBirthDeathShiftProcess::updateQmatrix(){

    const double &a = alpha ->getValue();
    const double &b = beta ->getValue();

    size_t n = sqrt(num_states);
    std::vector<size_t> A;
    std::vector<size_t> B;
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < n; j++){
            A.push_back(i);
            B.push_back(j);
        }
    }

    for (size_t i = 0; i < num_states; i++){
        for (size_t j = 0; j < num_states; j++){
            //std::cout << "i= " << i << "  j= " << j << std::endl;
            //std::cout << "Qith entry: \t " << Qmatrix(i,j) << std::endl;

            size_t n_changes = 0;
            if (A[i] != A[j]){
                n_changes += 1;
            }
            if (B[i] != B[j]){
                n_changes += 1;
            }

            if (n_changes == 1){
                if (A[i] != A[j]){
                    Qmatrix(i,j) = a / (n-1);
                }else{
                    Qmatrix(i,j) = b / (n-1);
                }
            }else{
                Qmatrix(i,j) = 0;
            }

        }
    }
    for (size_t i; i < num_states; i++){
        Qmatrix(i,i) = -(a+b);
    }
} 


/**
 * Resize various vectors depending on the current number of nodes.
 */
void FastBirthDeathShiftProcess::resizeVectors(size_t num_nodes)
{
    active_likelihood = std::vector<bool>(num_nodes, false);
    changed_nodes = std::vector<bool>(num_nodes, false);
    dirty_nodes = std::vector<bool>(num_nodes, true);
    node_partial_likelihoods = std::vector<std::vector<std::vector<double> > >(num_nodes, std::vector<std::vector<double> >(2,std::vector<double>(2*num_states,0)));
    scaling_factors = std::vector<std::vector<double> >(num_nodes, std::vector<double>(2,0.0) );
    average_speciation = std::vector<double>(num_nodes, 0.0);
    average_extinction = std::vector<double>(num_nodes, 0.0);
    num_shift_events = std::vector<long>(num_nodes, 0.0);
    time_in_states = std::vector<double>(num_states, 0.0);    
}
