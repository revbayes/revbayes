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
#include "PhyloCharacterEventDistribution.h"
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

PhyloCharacterEventDistribution::PhyloCharacterEventDistribution( const TypedDagNode<RbVector<double> >* root, const std::vector<TypedDistribution<double>* >& bd, const TypedDagNode<Tree>* t, const TypedDagNode<double> *s, const std::vector< std::string >& n) : TypedDistribution< CharacterHistoryDiscrete >( new CharacterHistoryDiscrete() ),
    root_values( root ),
    base_distribution( bd ),
    shift_rate( s ),
    tree( t ),
    names( n )
{
    
    num_values_per_event = names.size();
    
    if ( num_values_per_event != base_distribution.size() )
    {
        throw RbException("You need to provide the same number of base distributions as names in PhyloCharacterEventDistribution.");
    }

    if ( num_values_per_event != root_values->getValue().size() )
    {
        throw RbException("You need to provide the same number of root values as names in PhyloCharacterEventDistribution.");
    }

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    
    addParameter( shift_rate );
    addParameter( tree );
    
    // add all the root values
    addParameter( root_values );
    
    // add all the base distribution parameters
    for (size_t i=0; i<num_values_per_event; ++i)
    {
        // some base distributions can be NULL if this value should remain constant
        if ( base_distribution[i] != NULL )
        {
            const std::vector<const DagNode*>& sp_pars = base_distribution[i]->getParameters();
            for (std::vector<const DagNode*>::const_iterator it = sp_pars.begin(); it != sp_pars.end(); ++it)
            {
                this->addParameter( *it );
            }
        }
    }
    
    // connect the tree to the character history object
    value->setTree( &tree->getValue() );
    initializeBranchHistories( tree->getValue().getRoot(), tree->getValue().getRoot().getIndex() );
    
    // we also now initialize our value vectors with one element, which is the root category
    event_values = RbVector< RbVector<double> >( num_values_per_event, RbVector<double>(1,0.0) );
    for (size_t i=0; i<num_values_per_event; ++i)
    {
        event_values[i][0] = root_values->getValue()[i];
    }

}


PhyloCharacterEventDistribution::~PhyloCharacterEventDistribution()
{
    
}


/* Clone function */
PhyloCharacterEventDistribution* PhyloCharacterEventDistribution::clone( void ) const
{
    
    return new PhyloCharacterEventDistribution( *this );
}


/* Compute probability */
double PhyloCharacterEventDistribution::computeLnProbability( void )
{
    
    // Variable declarations and initialization
    double ln_prob = 0.0;
    
    // compute the probability at the root
    const TopologyNode &root = tree->getValue().getRoot();
    
    // fill the like
    const TopologyNode &left = root.getChild(0);
    size_t left_index = left.getIndex();
    double ln_prob_left = computeNodeProbability( left, left_index );
    const TopologyNode &right = root.getChild(1);
    size_t right_index = right.getIndex();
    double ln_prob_right = computeNodeProbability( right, right_index );
    
    
    // now compute the likelihoods of this internal node
    ln_prob = ln_prob_left + ln_prob_right;
    
    return ln_prob;
}


double PhyloCharacterEventDistribution::computeNodeProbability(const RevBayesCore::TopologyNode &node, size_t node_index)
{
    double ln_prob_node = 0.0;
    
    CharacterHistoryDiscrete& branch_histories = *value;

    // check for recomputation
    if ( true )
    {
                
        if ( node.isTip() == false )
        {
            
            // this is an internal node
            const TopologyNode &left = node.getChild(0);
            size_t left_index = left.getIndex();
            double ln_prob_left = computeNodeProbability( left, left_index );
            const TopologyNode &right = node.getChild(1);
            size_t right_index = right.getIndex();
            double ln_prob_right = computeNodeProbability( right, right_index );
            
            ln_prob_node += ln_prob_left + ln_prob_right;
            
        }
        
        // remember that we go back in time (rootwards)
        double begin_time = node.getAge();
        double branch_length = node.getBranchLength();
        double end_time = begin_time + branch_length;
        
        const BranchHistory& bh = branch_histories[ node_index ];
        const std::multiset<CharacterEvent*,CharacterEventCompare>& hist = bh.getHistory();
        for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=hist.begin(); it!=hist.end(); ++it)
        {
            
            for ( size_t i=0; i<num_values_per_event; ++i )
            {
                // we need to set the current rate category
                double current_value = computeStateValue( node.getIndex(), i, begin_time );
                TypedDistribution<double>* this_base_distribution = base_distribution[i];
                if ( this_base_distribution == NULL )
                {
                    if ( current_value != root_values->getValue()[i] )
                    {
                        ln_prob_node = RbConstants::Double::neginf;
                    }
                }
                else
                {
                    this_base_distribution->setValue( new double(current_value) );
                    ln_prob_node += this_base_distribution->computeLnProbability();
                }
            }
            
            CharacterEvent* event = *it;
            double event_time = event->getAge();
            
            // add the probability that there was shift event
            ln_prob_node -= shift_rate->getValue()*(event_time-begin_time);
            ln_prob_node += log( shift_rate->getValue() );

            begin_time = event_time;

        }
        
        // add the probability that there was no shift event in the remaining time
        ln_prob_node -= shift_rate->getValue()*(end_time-begin_time);
        
    }
    
    return ln_prob_node;
}


double PhyloCharacterEventDistribution::computeStartValue(size_t i, size_t value_index) const
{
    CharacterHistoryDiscrete& branch_histories = *value;

    size_t node_index = i;
    while ( tree->getValue().getNode(node_index).isRoot() == false && branch_histories[node_index].getNumberEvents() == 0)
    {
        node_index = tree->getValue().getNode(node_index).getParent().getIndex();
    }
    
    
    if ( tree->getValue().getNode(node_index).isRoot() == false )
    {
        const BranchHistory &bh = branch_histories[ node_index ];
        const std::multiset<CharacterEvent*, CharacterEventCompare> &h = bh.getHistory();
        
        CharacterEventDiscrete *event = static_cast<CharacterEventDiscrete*>(*h.begin());
        size_t event_index = event->getState();
        return event_values[value_index][event_index];
        
    }
    else
    {
        return root_values->getValue()[value_index];

        // this should be the same as
        // return event_values[value_index][0];
    }
    
}




double PhyloCharacterEventDistribution::computeStateValue(size_t i, size_t value_index, double time) const
{
    CharacterHistoryDiscrete& branch_histories = *value;

    size_t node_index = i;
    double event_value = 0;
    bool found = false;
    do
    {
        if ( tree->getValue().getNode(node_index).isRoot() == true )
        {

            event_value = root_values->getValue()[value_index];
            found = true;
        }
        else if ( branch_histories[node_index].getNumberEvents() > 0 )
        {
            
            const BranchHistory &bh = branch_histories[ node_index ];
            const std::multiset<CharacterEvent*, CharacterEventCompare> &h = bh.getHistory();
            for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=h.begin(); it!=h.end(); ++it)
            {
                CharacterEventDiscrete* event = static_cast<CharacterEventDiscrete*>(*it);
                double event_time = event->getAge();
                if ( event_time > time )
                {
                    found = true;
                    size_t event_index = event->getState();
                    event_value = event_values[value_index][event_index];
                    break;
                }
            }
            
            node_index = tree->getValue().getNode(node_index).getParent().getIndex();
            
        }
        else
        {
            node_index = tree->getValue().getNode(node_index).getParent().getIndex();
            
        }
        
        
    } while ( found == false );
    
    
    return event_value;
}


void PhyloCharacterEventDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<long> &rv) const
{
    CharacterHistoryDiscrete& branch_histories = *value;

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
        throw RbException() << "The PhyloCharacterEvent process does not have a member method called '" << n << "'.";
    }
    
}


void PhyloCharacterEventDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
    
    if ( n == "getValues" )
    {
        const std::string &val_name = static_cast<const TypedDagNode<std::string> * >( args[0] )->getValue();
        
        // need to find the index
        size_t index = 0;
        bool found = false;
        for (size_t i=0; i<num_values_per_event; ++i)
        {
            if ( names[i] == val_name )
            {
                found = true;
                index = i;
            }
        }
        
        if ( found == false )
        {
            throw RbException() << "Couldn't find values for name '" << val_name << "'.";
        }
            
        rv = event_values[index];
    }
    else if ( n == "getRealValues" )
    {
        const std::string &val_name = static_cast<const TypedDagNode<std::string> * >( args[0] )->getValue();
        
        // need to find the index
        size_t index = 0;
        bool found = false;
        for (size_t i=0; i<num_values_per_event; ++i)
        {
            if ( names[i] == val_name )
            {
                found = true;
                index = i;
            }
        }
        
        if ( found == false )
        {
            throw RbException() << "Couldn't find values for name '" << val_name << "'.";
        }
            
        rv = event_values[index];
    }
    else if ( n == "getRealPosValues" )
    {
        const std::string &val_name = static_cast<const TypedDagNode<std::string> * >( args[0] )->getValue();
        
        // need to find the index
        size_t index = 0;
        bool found = false;
        for (size_t i=0; i<num_values_per_event; ++i)
        {
            if ( names[i] == val_name )
            {
                found = true;
                index = i;
            }
        }
        
        if ( found == false )
        {
            throw RbException() << "Couldn't find values for name '" << val_name << "'.";
        }
            
        rv = event_values[index];
    }
    else
    {
        throw RbException() << "The Phylo-Character-Event Process does not have a member method called '" << n << "'.";
    }
    
}


/**
 * Get the character history object.
 */
CharacterHistoryDiscrete& PhyloCharacterEventDistribution::getCharacterHistory( void )
{
    CharacterHistoryDiscrete& branch_histories = *value;

    return branch_histories;
}

/**
 * Get the character history object.
 */
const CharacterHistoryDiscrete& PhyloCharacterEventDistribution::getCharacterHistory( void ) const
{
    CharacterHistoryDiscrete& branch_histories = *value;

    return branch_histories;
}


//TypedDistribution<double>* PhyloCharacterEventDistribution::getSpeciationRateDistibution(void) const
//{
//    return base_distribution_speciation;
//}
//
//
//double PhyloCharacterEventDistribution::getRootExtinctionRate(void) const
//{
//    return root_extinction->getValue();
//}



void PhyloCharacterEventDistribution::initializeBranchHistories(const TopologyNode &node, size_t nIdx)
{
    CharacterHistoryDiscrete& branch_histories = *value;

    if ( node.isRoot() == false )
    {
        CharacterHistoryDiscrete& ch = branch_histories;
        const BranchHistoryDiscrete& bh = ch[nIdx];
        const std::vector<CharacterEvent*>& parents = bh.getParentCharacters();
        const std::vector<CharacterEvent*>& children = bh.getParentCharacters();
    
        if ( parents.size() > 0 )
        {
            CharacterEvent *ce = parents[0];
            static_cast<CharacterEventDiscrete*>( ce )->setState( 0 );
        }
        if ( children.size() > 0 )
        {
            CharacterEvent *ce = children[0];
            static_cast<CharacterEventDiscrete*>( ce )->setState( 0 );
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
        initializeBranchHistories( left, left_index );
        const TopologyNode &right = node.getChild(1);
        size_t right_index = right.getIndex();
        initializeBranchHistories( right, right_index );
        
        
    }
    
}


void PhyloCharacterEventDistribution::redrawValue( void )
{
    // this code should be written (TODO: @Priscilla)
}


void PhyloCharacterEventDistribution::setValue(CharacterHistoryDiscrete *v, bool force)
{
    
    // delegate to the parent class
    TypedDistribution< CharacterHistoryDiscrete >::setValue(v, force);
    
    CharacterHistoryDiscrete& branch_histories = *value;
    branch_histories.setTree( &tree->getValue() );

    initializeBranchHistories( tree->getValue().getRoot(), tree->getValue().getRoot().getIndex() );

}


void PhyloCharacterEventDistribution::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    // if the tree has changed, so will have the character history object (which includes the tree)
    // thus, the value of this distribution has changed and we keep flaggin downstream DAG nodes
    if ( affecter == tree && this->dag_node != NULL )
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
    // if one of the root values has changed (marked as dirty), then we need to pass down this message
    if ( affecter == root_values && this->dag_node != NULL )
    {
        dag_node->initiateGetAffectedNodes( affected );
    }

    
}

/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void PhyloCharacterEventDistribution::keepSpecialization(const DagNode *affecter)
{
    
    // if the tree was marked to be kept (not dirty anymore), so should the character history object as well (which includes the tree)
    // thus, we keep flagging downstream DAG nodes
    if ( affecter == tree && this->dag_node != NULL )
    {
        dag_node->keepAffected();
    }
    
    // if one of the root values has changed and is now called to be kept (accepted), then we need to pass down this message
    if ( affecter == root_values && this->dag_node != NULL )
    {
        dag_node->keepAffected();
    }

    
}

/**
 * Restore the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void PhyloCharacterEventDistribution::restoreSpecialization(const DagNode *affecter)
{
    
    // if the tree was marked to be restored, so should the character history object as well (which includes the tree)
    // thus, we keep sending this message to downstream DAG nodes
    if ( affecter == tree && this->dag_node != NULL )
    {
        dag_node->restoreAffected();
    }
    
    // if one of the root values has changed, then we need to change the first element in the value vectors too
    
    if ( affecter == root_values )
    {
        for (size_t i=0; i<num_values_per_event; ++i)
        {
            event_values[i][0] = root_values->getValue()[i];
        }
    }
    
}


/** Swap a parameter of the distribution */
void PhyloCharacterEventDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == tree)
    {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
        
        // connect the tree to the character history object
        value->setTree( &tree->getValue() );
        initializeBranchHistories( tree->getValue().getRoot(), tree->getValue().getRoot().getIndex() );
    }
    
    if (oldP == shift_rate)
    {
        shift_rate = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    // if one of the root values has changed
    if ( oldP == root_values )
    {
        root_values = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
    // if one of the base distribution parameters has changed
    for (size_t i=0; i<num_values_per_event; ++i)
    {
        if ( base_distribution[i] != NULL )
        {
            const std::vector<const DagNode*>& this_dist_pars = base_distribution[i]->getParameters();
            bool is_this_dist_par = false;
            for (std::vector<const DagNode*>::const_iterator it = this_dist_pars.begin(); it != this_dist_pars.end(); ++it)
            {
                if ( *it == oldP )
                {
                    is_this_dist_par = true;
                    break;
                }
            }
            if ( is_this_dist_par == true )
            {
                base_distribution[i]->swapParameter(oldP,newP);
            }
        }
    }
    
}

/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void PhyloCharacterEventDistribution::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    
    // if the tree has changed and touched/fladded as dirty, so should the character history object as well (which includes the tree)
    // thus, we keep sending this message to downstream DAG nodes
    if ( affecter == tree && this->dag_node != NULL )
    {
        dag_node->restoreAffected();
    }
    
    // if one of the root values has changed, then we need to change the first element in the value vectors too
    if ( affecter == root_values )
    {
        for (size_t i=0; i<num_values_per_event; ++i)
        {
            event_values[i][0] = root_values->getValue()[i];
        }
    }
    
}
