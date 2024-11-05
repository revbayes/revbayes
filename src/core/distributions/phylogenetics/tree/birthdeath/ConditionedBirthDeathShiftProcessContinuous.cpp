#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "ConditionedBirthDeathShiftProcessContinuous.h"
#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "BranchHistory.h"
#include "BranchHistoryContinuous.h"
#include "CharacterEvent.h"
#include "CharacterEventCompare.h"
#include "CharacterEventContinuous.h"
#include "CharacterHistoryContinuous.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }

using namespace RevBayesCore;

ConditionedBirthDeathShiftProcessContinuous::ConditionedBirthDeathShiftProcessContinuous( const TypedDagNode<double> *a, const TypedDagNode<double> *root_sp, const TypedDagNode<double> *root_ex, TypedDistribution<double> *s, TypedDistribution<double> *e, const TypedDagNode<double > *ev, const TypedDagNode< double > *r, const std::string &cdt, const std::vector<Taxon> &n) : AbstractCharacterHistoryBirthDeathProcess(),
    root_age( a ),
    root_speciation( root_sp ),
    root_extinction( root_ex ),
    base_distribution_speciation( s ),
    base_distribution_extinction( e ),
    shift_rate( ev ),
    rho( r ),
    branch_histories( NULL, 1 ),
    condition( cdt ),
    taxa( n )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( root_age );
    addParameter( root_speciation );
    addParameter( root_extinction );
    addParameter( shift_rate );
    addParameter( rho );
    
    // add the parameters of the distribution
    if ( base_distribution_speciation != NULL )
    {
        const std::vector<const DagNode*>& sp_pars = base_distribution_speciation->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = sp_pars.begin(); it != sp_pars.end(); ++it)
        {
            this->addParameter( *it );
        }
    }
    
    if ( base_distribution_extinction != NULL )
    {
        const std::vector<const DagNode*>& ex_pars = base_distribution_extinction->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = ex_pars.begin(); it != ex_pars.end(); ++it)
        {
            this->addParameter( *it );
        }
    }

    
    num_taxa = taxa.size();
    
    // the combinatorial factor for the probability of a labelled history is
    // 2^{n-1} / ( n! * (n-1)! )
    // but since the probability of the divergence times contains the factor (n-1)! we simply store
    // 2^{n-1} / n!
    double lnFact = 0.0;
    for (size_t i = 2; i <= num_taxa; i++)
    {
        lnFact += std::log(i);
    }
    
    event_prior_only = false;
    
    if ( event_prior_only == false )
    {
        log_tree_topology_prob = (num_taxa - 1) * RbConstants::LN2 - lnFact ;
    }
    else
    {
        log_tree_topology_prob = 0.0;
    }
    
    simulateTree();

}


ConditionedBirthDeathShiftProcessContinuous::~ConditionedBirthDeathShiftProcessContinuous()
{
    
}


/**
 * Recursive call to attach ordered interior node times to the time tree psi. Call it initially with the
 * root of the tree.
 */
void ConditionedBirthDeathShiftProcessContinuous::attachTimes(Tree* psi, std::vector<TopologyNode *> &nodes, size_t index, const std::vector<double> &interiorNodeTimes, double originTime )
{
    
    if (index < num_taxa-1)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;
        
        // Randomly draw one node from the list of nodes
        size_t node_index = static_cast<size_t>( floor(rng->uniform01()*nodes.size()) );
        
        // Get the node from the list
        TopologyNode* parent = nodes.at(node_index);
        psi->getNode( parent->getIndex() ).setAge( originTime - interiorNodeTimes[index] );
        
        // Remove the randomly drawn node from the list
        nodes.erase(nodes.begin()+std::int64_t(node_index));
        
        // Add the left child if an interior node
        TopologyNode* leftChild = &parent->getChild(0);
        if ( !leftChild->isTip() )
        {
            nodes.push_back(leftChild);
        }
        
        // Add the right child if an interior node
        TopologyNode* rightChild = &parent->getChild(1);
        if ( !rightChild->isTip() )
        {
            nodes.push_back(rightChild);
        }
        
        // Recursive call to this function
        attachTimes(psi, nodes, index+1, interiorNodeTimes, originTime);
    }
    
}


/** Build random binary tree to size num_taxa. The result is a draw from the uniform distribution on histories. */
void ConditionedBirthDeathShiftProcessContinuous::buildRandomBinaryHistory(std::vector<TopologyNode*> &tips)
{
    
    if (tips.size() < num_taxa)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;
        
        // Randomly draw one node from the list of tips
        size_t index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );
        
        // Get the node from the list
        TopologyNode* parent = tips.at(index);
        
        // Remove the randomly drawn node from the list
        tips.erase(tips.begin()+std::int64_t(index));
        
        // Add a left child
        TopologyNode* leftChild = new TopologyNode();
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        tips.push_back(leftChild);
        
        // Add a right child
        TopologyNode* rightChild = new TopologyNode();
        parent->addChild(rightChild);
        rightChild->setParent(parent);
        tips.push_back(rightChild);
        
        // Recursive call to this function
        buildRandomBinaryHistory(tips);
    }
}


/* Clone function */
ConditionedBirthDeathShiftProcessContinuous* ConditionedBirthDeathShiftProcessContinuous::clone( void ) const
{
    
    return new ConditionedBirthDeathShiftProcessContinuous( *this );
}


double ConditionedBirthDeathShiftProcessContinuous::computeBranchProbability(double begin, double end, double sp_rate, double ex_rate, double sh_rate)
{
    
    double ln_prob = 0.0;
    
    double sampling = rho->getValue();
    
    // add the probability that there was no shiftevent
    ln_prob -= sh_rate*(end-begin);

    // we need to compute P(N(T)=1|N(start)=1) / P(N(T)=1|N(end)=1)
    // where                      (b-d)^2 * exp(-(b-d)(T-t))
    //         P(N(T)=1|N(t)=1) = --------------------------
    //                            (b-d*exp(-(b-d)(T-t))^2
    //
    // or with incomplete sampling:
    // where                      r(b-d)^2 * exp(-(b-d)(T-t))
    //         P(N(T)=1|N(t)=1) = --------------------------
    //                            (rb+(b*(1-r)-d)*exp(-(b-d)(T-t)))^2

//    ln_prob += 2*log(sp_rate-ex_rate) - ((sp_rate-ex_rate)*(present-begin));
//    ln_prob -= 2*log(sp_rate - ex_rate*exp((sp_rate-ex_rate)*(present-begin)));
//
//    ln_prob -= 2*log(sp_rate-ex_rate) - ((sp_rate-ex_rate)*(present-end));
//    ln_prob += 2*log(sp_rate - ex_rate*exp((sp_rate-ex_rate)*(present-end)));

    if ( event_prior_only == false )
    {
        ln_prob += log(sampling) + 2*log(sp_rate-ex_rate) - ((sp_rate-ex_rate)*(end));
        ln_prob -= 2*log(sampling*sp_rate + (sp_rate*(1-sampling) - ex_rate)*exp((ex_rate-sp_rate)*(end)));

        ln_prob -= log(sampling) + 2*log(sp_rate-ex_rate) - ((sp_rate-ex_rate)*(begin));
        ln_prob += 2*log(sampling*sp_rate + (sp_rate*(1-sampling) - ex_rate)*exp((ex_rate-sp_rate)*(begin)));
    }
    
    return ln_prob;
    
}


/* Compute probability */
double ConditionedBirthDeathShiftProcessContinuous::computeLnProbability( void )
{
    
    // Variable declarations and initialization
    double ln_prob = 0.0;
    double age = root_age->getValue();
    
    // we need to check that the root age matches
    if ( age != value->getRoot().getAge() )
    {
        return RbConstants::Double::neginf;
    }
    
    
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
            else if ( (the_node.getAge() - (*it)->getParent().getAge()) > 1E-6 && the_node.isSampledAncestorTip() == true )
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
            else if ( the_node.getBranchLength() > 1E-6 )
            {
                return RbConstants::Double::neginf;
            }
            
        }
        
    }
    
    // compute the probability at the root
    ln_prob = computeRootLikelihood();
    
    return ln_prob + log_tree_topology_prob;
}


double ConditionedBirthDeathShiftProcessContinuous::computeNodeProbability(const RevBayesCore::TopologyNode &node, size_t node_index)
{
    double ln_prob_node = 0.0;
    
    // check for recomputation
    if ( true )
    {
        
        double value_tipwards_birth  = computeStartValue( node_index, 0 );
        
        if ( node.isTip() )
        {
            // this is a tip node
            double sampling_probability = rho->getValue();
            
            ln_prob_node += log( sampling_probability );
        }
        else
        {
            // this is an internal node
            const TopologyNode &left = node.getChild(0);
            size_t left_index = left.getIndex();
            double ln_prob_left = computeNodeProbability( left, left_index );
            const TopologyNode &right = node.getChild(1);
            size_t right_index = right.getIndex();
            double ln_prob_right = computeNodeProbability( right, right_index );
            
            if ( event_prior_only == false )
            {
                ln_prob_node += ln_prob_left + ln_prob_right + log(value_tipwards_birth);
            }
            else
            {
                ln_prob_node += ln_prob_left + ln_prob_right;
            }

        }
        
        // remember that we go back in time (rootwards)
        double begin_time = node.getAge();
        double branch_length = node.getBranchLength();
        double end_time = begin_time + branch_length;
        
        const BranchHistory& bh = branch_histories[ node_index ];
        const std::multiset<CharacterEvent*,CharacterEventCompare>& hist = bh.getHistory();
        for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=hist.begin(); it!=hist.end(); ++it)
        {
            
            // we need to set the current rate category
            double current_value_birth = computeStateValue( node.getIndex(), 0, begin_time );
            double current_value_death = computeStateValue( node.getIndex(), 1, begin_time );
            
            if ( base_distribution_speciation == NULL )
            {
                if ( current_value_birth != root_speciation->getValue() )
                {
                    ln_prob_node = RbConstants::Double::neginf;
                }
            }
            else
            {
                base_distribution_speciation->setValue( new double(current_value_birth) );
                ln_prob_node += base_distribution_speciation->computeLnProbability();
            }
            
            if ( base_distribution_extinction == NULL )
            {
                if ( current_value_death != root_extinction->getValue() )
                {
                    ln_prob_node = RbConstants::Double::neginf;
                }
            }
            else
            {
                base_distribution_extinction->setValue( new double(current_value_death) );
                ln_prob_node += base_distribution_extinction->computeLnProbability();
            }
            
            CharacterEvent* event = *it;
            double event_time = event->getAge();

            ln_prob_node += computeBranchProbability(begin_time, event_time, current_value_birth, current_value_death, shift_rate->getValue() );
            ln_prob_node += log( shift_rate->getValue() );

            begin_time = event_time;

        }

        double value_rootwards_birth = computeStartValue( node.getParent().getIndex(), 0 );
        double value_rootwards_death = computeStartValue( node.getParent().getIndex(), 1 );

        ln_prob_node += computeBranchProbability(begin_time, end_time, value_rootwards_birth, value_rootwards_death, shift_rate->getValue());
        
    }
    
    return ln_prob_node;
}


double ConditionedBirthDeathShiftProcessContinuous::computeStartValue(size_t i, size_t j) const
{
    
    size_t node_index = i;
    while ( value->getNode(node_index).isRoot() == false && branch_histories[node_index].getNumberEvents() == 0)
    {
        node_index = value->getNode(node_index).getParent().getIndex();
    }
    
    
    if ( value->getNode(node_index).isRoot() == false )
    {
        const BranchHistory &bh = branch_histories[ node_index ];
        const std::multiset<CharacterEvent*, CharacterEventCompare> &h = bh.getHistory();
        
        CharacterEventContinuous *event = static_cast<CharacterEventContinuous*>(*h.begin());
        return event->getState(j);
        
    }
    else
    {
        if ( j == 0 )
        {
            return root_speciation->getValue();
        }
        else
        {
            return root_extinction->getValue();
        }
        
    }
    
}




double ConditionedBirthDeathShiftProcessContinuous::computeStateValue(size_t i, size_t j, double time) const
{
    
    size_t node_index = i;
    double event_value = 0;
    bool found = false;
    do
    {
        if ( value->getNode(node_index).isRoot() == true )
        {

            if ( j == 0 )
            {
                event_value = root_speciation->getValue();
            }
            else
            {
                event_value = root_extinction->getValue();
            }
            found = true;
        }
        else if ( branch_histories[node_index].getNumberEvents() > 0 )
        {
            
            const BranchHistory &bh = branch_histories[ node_index ];
            const std::multiset<CharacterEvent*, CharacterEventCompare> &h = bh.getHistory();
            for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=h.begin(); it!=h.end(); ++it)
            {
                CharacterEventContinuous* event = static_cast<CharacterEventContinuous*>(*it);
                double event_time = event->getAge();
                if ( event_time > time )
                {
                    found = true;
                    event_value = event->getState(j);
                    break;
                }
            }
            
            node_index = value->getNode(node_index).getParent().getIndex();
            
        }
        else
        {
            node_index = value->getNode(node_index).getParent().getIndex();
            
        }
        
        
    } while ( found == false );
    
    
    return event_value;
}


double ConditionedBirthDeathShiftProcessContinuous::computeRootLikelihood( void )
{
    
    const TopologyNode &root = value->getRoot();
    
    // fill the like
    const TopologyNode &left = root.getChild(0);
    size_t left_index = left.getIndex();
    double ln_prob_left = computeNodeProbability( left, left_index );
    const TopologyNode &right = root.getChild(1);
    size_t right_index = right.getIndex();
    double ln_prob_right = computeNodeProbability( right, right_index );
    
    
    // now compute the likelihoods of this internal node
    double ln_prob = ln_prob_left + ln_prob_right;
    
    if ( condition == "survival" )
    {
//        ln_prob -= 2*log( 1.0-leftStates[ root_state->getValue()-1 ] );
    }
    
    return ln_prob;
}


void ConditionedBirthDeathShiftProcessContinuous::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<std::int64_t> &rv) const
{
    
    if ( n == "numberEvents" )
    {
        size_t num_branches = branch_histories.getNumberBranches();
        rv.clear();
        rv.resize( num_branches );
        
        for (size_t i = 0; i < num_branches; ++i)
        {
            rv[i] = int(branch_histories[i].getNumberEvents());
        }
        
    }
    else
    {
        throw RbException("The birth-death-shift process does not have a member method called '" + n + "'.");
    }
    
}


void ConditionedBirthDeathShiftProcessContinuous::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
    
    if ( n == "averageSpeciationRate" )
    {
        size_t num_branches = branch_histories.getNumberBranches();
        rv.clear();
        rv.resize( num_branches );
        
        for (size_t i = 0; i < num_branches; ++i)
        {
            const TopologyNode &node = this->value->getNode( i );
            const BranchHistory& bh = branch_histories[ i ];
            const std::multiset<CharacterEvent*,CharacterEventCompare>& hist = bh.getHistory();
            double value_rootwards = computeStartValue( node.getParent().getIndex(), 0 );
            
            double rate = 0;
            double begin_time = 0.0;
            double branch_length = node.getBranchLength();
            for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=hist.begin(); it!=hist.end(); ++it)
            {
                CharacterEventContinuous* event = static_cast<CharacterEventContinuous*>(*it);
                double end_time = event->getAge() - node.getAge();
                double time_interval = (end_time - begin_time) / branch_length;
                
                // we need to set the current rate caterogy
                double current_rate = event->getState(0);
                
                rate += time_interval * current_rate;
                
                begin_time = end_time;
            }
            rate += (branch_length-begin_time)/branch_length * value_rootwards;
            
            rv[i] = rate;
            
        }
        
    }
    else if ( n == "averageExtinctionRate" )
    {
        size_t num_branches = branch_histories.getNumberBranches();
        rv.clear();
        rv.resize( num_branches );
        
        for (size_t i = 0; i < num_branches; ++i)
        {
            const TopologyNode &node = this->value->getNode( i );
            const BranchHistory& bh = branch_histories[ i ];
            const std::multiset<CharacterEvent*,CharacterEventCompare>& hist = bh.getHistory();
            double value_rootwards = computeStartValue( node.getParent().getIndex(), 1 );
            
            double rate = 0;
            double begin_time = 0.0;
            double branch_length = node.getBranchLength();
            for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=hist.begin(); it!=hist.end(); ++it)
            {
                CharacterEventContinuous* event = static_cast<CharacterEventContinuous*>(*it);
                double end_time = event->getAge() - node.getAge();
                double time_interval = (end_time - begin_time) / branch_length;
                
                // we need to set the current rate caterogy
                double current_rate = event->getState(1);
                
                rate += time_interval * current_rate;
                
                begin_time = end_time;
            }
            rate += (branch_length-begin_time)/branch_length * value_rootwards;
            
            rv[i] = rate;
            
        }
        
    }
    else if ( n == "probabilitySpeciationRateIncrease" )
    {
        size_t num_branches = branch_histories.getNumberBranches();
        rv.clear();
        rv.resize( num_branches );
        
        for (size_t i = 0; i < num_branches; ++i)
        {
            const TopologyNode &node = this->value->getNode( i );
            const BranchHistory& bh = branch_histories[ i ];
            // const std::multiset<CharacterEvent*,CharacterEventCompare>& hist = bh.getHistory();
            bh.getHistory();
            double value_rootwards = computeStartValue( node.getParent().getIndex(), 0 );
            double value_tipwards  = computeStartValue( node.getIndex(), 0 );
            
            rv[i] = ( value_rootwards < value_tipwards ? 1.0 : 0.0 );
            
        }
        
    }
    else if ( n == "probabilitySpeciationRateDecrease" )
    {
        size_t num_branches = branch_histories.getNumberBranches();
        rv.clear();
        rv.resize( num_branches );
        
        for (size_t i = 0; i < num_branches; ++i)
        {
            const TopologyNode &node = this->value->getNode( i );
            const BranchHistory& bh = branch_histories[ i ];
//            const std::multiset<CharacterEvent*,CharacterEventCompare>& hist = bh.getHistory();
            bh.getHistory();
            double value_rootwards = computeStartValue( node.getParent().getIndex(), 0 );
            double value_tipwards  = computeStartValue( node.getIndex(), 0 );
            
            rv[i] = ( value_rootwards > value_tipwards ? 1.0 : 0.0 );
            
        }
        
    }
    else if ( n == "probabilityExtinctionRateIncrease" )
    {
        size_t num_branches = branch_histories.getNumberBranches();
        rv.clear();
        rv.resize( num_branches );
        
        for (size_t i = 0; i < num_branches; ++i)
        {
            const TopologyNode &node = this->value->getNode( i );
            const BranchHistory& bh = branch_histories[ i ];
//            const std::multiset<CharacterEvent*,CharacterEventCompare>& hist = bh.getHistory();
            bh.getHistory();
            double value_rootwards = computeStartValue( node.getParent().getIndex(), 1 );
            double value_tipwards  = computeStartValue( node.getIndex(), 1 );
            
            rv[i] = ( value_rootwards < value_tipwards ? 1.0 : 0.0 );
            
        }
        
    }
    else if ( n == "probabilityExtinctionRateDecrease" )
    {
        size_t num_branches = branch_histories.getNumberBranches();
        rv.clear();
        rv.resize( num_branches );
        
        for (size_t i = 0; i < num_branches; ++i)
        {
            const TopologyNode &node = this->value->getNode( i );
            const BranchHistory& bh = branch_histories[ i ];
//            const std::multiset<CharacterEvent*,CharacterEventCompare>& hist = bh.getHistory();
            bh.getHistory();
            double value_rootwards = computeStartValue( node.getParent().getIndex(), 1 );
            double value_tipwards  = computeStartValue( node.getIndex(), 1 );
            
            rv[i] = ( value_rootwards > value_tipwards ? 1.0 : 0.0 );
            
        }
        
    }
    else
    {
        throw RbException("The heterogeneous rate birth-death process does not have a member method called '" + n + "'.");
    }
    
}


/**
 * Get the character history object.
 */
CharacterHistoryContinuous& ConditionedBirthDeathShiftProcessContinuous::getCharacterHistory( void )
{
    
    return branch_histories;
}

/**
 * Get the character history object.
 */
const CharacterHistoryContinuous& ConditionedBirthDeathShiftProcessContinuous::getCharacterHistory( void ) const
{
    
    return branch_histories;
}


TypedDistribution<double>* ConditionedBirthDeathShiftProcessContinuous::getExtinctionRateDistibution(void) const
{
    return base_distribution_extinction;
}


TypedDistribution<double>* ConditionedBirthDeathShiftProcessContinuous::getSpeciationRateDistibution(void) const
{
    return base_distribution_speciation;
}


double ConditionedBirthDeathShiftProcessContinuous::getRootExtinctionRate(void) const
{
    return root_extinction->getValue();
}


double ConditionedBirthDeathShiftProcessContinuous::getRootSpeciationRate(void) const
{
    return root_speciation->getValue();
}



void ConditionedBirthDeathShiftProcessContinuous::initializeBranchHistories(const TopologyNode &node, size_t nIdx, size_t idx, double val)
{
    
    if ( node.isRoot() == false )
    {
        CharacterHistoryContinuous& ch = branch_histories;
        const BranchHistoryContinuous& bh = ch[nIdx];
        const std::vector<CharacterEvent*>& parents = bh.getParentCharacters();
        const std::vector<CharacterEvent*>& children = bh.getParentCharacters();
    
        if ( parents.size() > 0 )
        {
            CharacterEvent *ce = parents[0];
            static_cast<CharacterEventContinuous*>( ce )->resize( 2 );
            static_cast<CharacterEventContinuous*>( ce )->setState( val, idx );
        }
        if ( children.size() > 0 )
        {
            CharacterEvent *ce = children[0];
            static_cast<CharacterEventContinuous*>( ce )->resize( 2 );
            static_cast<CharacterEventContinuous*>( ce )->setState( val, idx );
        }
    }
    
    if ( node.isTip() )
    {
        
    }
    else
    {
        // this is an internal node
        const TopologyNode &left = node.getChild(0);
        size_t left_index = left.getIndex();
        initializeBranchHistories( left, left_index, idx, val );
        const TopologyNode &right = node.getChild(1);
        size_t right_index = right.getIndex();
        initializeBranchHistories( right, right_index, idx, val );
        
        
    }
    
}


bool ConditionedBirthDeathShiftProcessContinuous::isExtinctionRateConstant( void ) const
{
    return base_distribution_extinction == NULL;
}


bool ConditionedBirthDeathShiftProcessContinuous::isSpeciationRateConstant( void ) const
{
    return base_distribution_speciation == NULL;
}


void ConditionedBirthDeathShiftProcessContinuous::redrawValue( void )
{
    simulateTree();
}


void ConditionedBirthDeathShiftProcessContinuous::setValue(Tree *v, bool force)
{
    
    // delegate to the parent class
    TypedDistribution< Tree >::setValue(v, force);
    
    branch_histories.setTree( value );

    initializeBranchHistories( value->getRoot(), value->getRoot().getIndex(), 0, root_speciation->getValue() );
    initializeBranchHistories( value->getRoot(), value->getRoot().getIndex(), 1, root_extinction->getValue() );

}


/** Simulate the tree conditioned on the time of origin */
void ConditionedBirthDeathShiftProcessContinuous::simulateTree( void )
{
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // Create the time tree object (topology + times)
    Tree *psi = new Tree();
    
    // Root the topology by setting the appropriate flag
    psi->setRooted( true );
    
    // Create the root node and a vector of nodes
    TopologyNode* root = new TopologyNode();
    std::vector<TopologyNode* > nodes;
    nodes.push_back(root);
    
    // Draw a random tree history
    buildRandomBinaryHistory(nodes);
    
    // Set the tip names
    for (size_t i=0; i<num_taxa; i++)
    {
        size_t index = size_t( floor(rng->uniform01() * nodes.size()) );
        
        // Get the node from the list
        TopologyNode* node = nodes.at(index);
        
        // Remove the randomly drawn node from the list
        nodes.erase(nodes.begin()+std::int64_t(index) );
        
        // Set taxon
        node->setTaxon( taxa[i] );
    }
    
    // Initialize the topology by setting the root
    psi->setRoot(root, true);
    
    // Now simulate the speciation times counted from originTime
    std::vector<double> intNodeTimes;
    double              t_or = root_age->getValue();
    intNodeTimes.push_back( 0.0 );  // For time of mrca
    for ( size_t i=0; i<num_taxa-2; i++ )
    {
        double t = rng->uniform01() * t_or;
        intNodeTimes.push_back( t );
    }
    
    // Sort the speciation times from 0.0 (root node) to the largest value
    std::sort( intNodeTimes.begin(), intNodeTimes.end() );
    
    // Attach times
    nodes.clear();
    nodes.push_back( root );
    attachTimes(psi, nodes, 0, intNodeTimes, t_or);
    for (size_t i = 0; i < num_taxa; ++i)
    {
        TopologyNode& node = psi->getTipNode(i);
        psi->getNode( node.getIndex() ).setAge( 0.0 );
    }
    
    
    // Finally store the new value
    value = psi;
    
    branch_histories.setTree( psi );
    
    initializeBranchHistories( value->getRoot(), value->getRoot().getIndex(), 0, root_speciation->getValue() );
    initializeBranchHistories( value->getRoot(), value->getRoot().getIndex(), 1, root_extinction->getValue() );

}

void ConditionedBirthDeathShiftProcessContinuous::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    if ( affecter == root_age && this->dag_node != NULL )
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
}

/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void ConditionedBirthDeathShiftProcessContinuous::keepSpecialization(const DagNode *affecter)
{
    
    if ( affecter == root_age && this->dag_node != NULL )
    {
        dag_node->keepAffected();
    }
    
}

/**
 * Restore the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void ConditionedBirthDeathShiftProcessContinuous::restoreSpecialization(const DagNode *affecter)
{
    
    if ( affecter == root_age && this->dag_node != NULL )
    {
        value->getNode( value->getRoot().getIndex() ).setAge( root_age->getValue() );
        dag_node->restoreAffected();
    }
    
    if ( affecter == root_speciation && isSpeciationRateConstant() == true )
    {
        CharacterHistoryContinuous& ch = branch_histories;
        size_t num_branches = ch.getNumberBranches();
        
        double lambda = root_speciation->getValue();
        
        for (size_t i=0; i<num_branches; ++i)
        {
            BranchHistoryContinuous& bh = ch[i];
            for ( size_t j=0; j<bh.getNumberEvents(); ++j )
            {
                bh.getEvent(j)->setState(lambda,0);
            }
            
        }
        
    }
    
    if ( affecter == root_extinction && isExtinctionRateConstant() == true )
    {
        CharacterHistoryContinuous& ch = branch_histories;
        size_t num_branches = ch.getNumberBranches();
        
        double mu = root_extinction->getValue();
        
        for (size_t i=0; i<num_branches; ++i)
        {
            BranchHistoryContinuous& bh = ch[i];
            for ( size_t j=0; j<bh.getNumberEvents(); ++j )
            {
                bh.getEvent(j)->setState(mu,1);
            }
            
        }
        
    }
    
}


void ConditionedBirthDeathShiftProcessContinuous::runMCMC( void )
{
    
    
    
}

/** Swap a parameter of the distribution */
void ConditionedBirthDeathShiftProcessContinuous::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == root_age)
    {
        root_age = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == root_speciation)
    {
        root_speciation = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == root_extinction)
    {
        root_extinction = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == shift_rate)
    {
        shift_rate = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if ( base_distribution_speciation != NULL )
    {
        const std::vector<const DagNode*>& sp_pars = base_distribution_speciation->getParameters();
        bool is_speciation_par = false;
        for (std::vector<const DagNode*>::const_iterator it = sp_pars.begin(); it != sp_pars.end(); ++it)
        {
            if ( *it == oldP )
            {
                is_speciation_par = true;
                break;
            }
        }
        if ( is_speciation_par == true )
        {
            base_distribution_speciation->swapParameter(oldP,newP);
        }
        
    }
    
    if ( base_distribution_extinction != NULL )
    {
        
        const std::vector<const DagNode*>& ex_pars = base_distribution_extinction->getParameters();
        bool is_extinction_par = false;
        for (std::vector<const DagNode*>::const_iterator it = ex_pars.begin(); it != ex_pars.end(); ++it)
        {
            if ( *it == oldP )
            {
                is_extinction_par = true;
                break;
            }
        }
        if ( is_extinction_par == true )
        {
            base_distribution_extinction->swapParameter(oldP,newP);
        }
        
    }
    
}

/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void ConditionedBirthDeathShiftProcessContinuous::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    
    if ( affecter == root_age )
    {
        value->getNode( value->getRoot().getIndex() ).setAge( root_age->getValue() );
        
        if ( this->dag_node != NULL )
        {
            dag_node->touchAffected();
        }
        
    }
    
    if ( affecter == root_speciation && isSpeciationRateConstant() == true )
    {
        CharacterHistoryContinuous& ch = branch_histories;
        size_t num_branches = ch.getNumberBranches();
        
        double lambda = root_speciation->getValue();
        
        for (size_t i=0; i<num_branches; ++i)
        {
            BranchHistoryContinuous& bh = ch[i];
            for ( size_t j=0; j<bh.getNumberEvents(); ++j )
            {
                bh.getEvent(j)->setState(lambda,0);
            }
            
        }
        
    }
    
    if ( affecter == root_extinction && isExtinctionRateConstant() == true )
    {
        CharacterHistoryContinuous& ch = branch_histories;
        size_t num_branches = ch.getNumberBranches();
        
        double mu = root_extinction->getValue();
        
        for (size_t i=0; i<num_branches; ++i)
        {
            BranchHistoryContinuous& bh = ch[i];
            for ( size_t j=0; j<bh.getNumberEvents(); ++j )
            {
                bh.getEvent(j)->setState(mu,1);
            }
            
        }
        
    }
    
}
