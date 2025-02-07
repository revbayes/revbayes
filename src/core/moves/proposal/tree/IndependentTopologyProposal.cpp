#include "IndependentTopologyProposal.h"

#include <algorithm>
#include <cstddef>
#include <string>

#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

IndependentTopologyProposal::IndependentTopologyProposal( StochasticNode<Tree> *v, TypedDistribution<Tree>* d,  StochasticNode< RbVector<double> > *r, Clade og ) : Proposal(),
    outgroup( og ),
    proposal_distribution( d ),
    rates_node( r ),
    variable( v ),
    stored_ln_num_rankings( 0.0 ),
    init( false )
{

    addNode( rates_node );
    addNode( variable );

    outgroup.resetTaxonBitset( variable->getValue().getTaxonBitSetMap() );

    recursivelyRank( variable->getValue().getRoot(), stored_ln_num_rankings );

}


IndependentTopologyProposal::IndependentTopologyProposal( StochasticNode<Tree> *v, TypedDistribution<Tree>* d,  std::vector<StochasticNode<double>* >r, Clade og ) : Proposal(),
    outgroup( og ),
    proposal_distribution( d ),
    rates_node( NULL ),
    rates_vector( r ),
    variable( v ),
    stored_ln_num_rankings( 0.0 ),
    init( false )
{
    
    addNode( variable );

    for(size_t i = 0; i < rates_vector.size(); i++)
    {
        addNode( rates_vector[i] );
    }

    outgroup.resetTaxonBitset( variable->getValue().getTaxonBitSetMap() );

    recursivelyRank( variable->getValue().getRoot(), stored_ln_num_rankings );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
IndependentTopologyProposal* IndependentTopologyProposal::clone( void ) const
{
    
    return new IndependentTopologyProposal( *this );
}


/** Clean the proposal */
void IndependentTopologyProposal::cleanProposal( void )
{
    stored_ln_num_rankings = proposal_ln_num_rankings;
}

const std::string& IndependentTopologyProposal::getProposalName( void ) const
{
    static std::string name = "IndependentTopologyMove";
    
    return name;
}


double IndependentTopologyProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/** Perform the move */
double IndependentTopologyProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    double hr = 0.0;

    stored_tree = variable->getValue();

    // tree proposal ratio numerator
    hr += proposal_distribution->computeLnProbability();

    // draw a new tree
    proposal_distribution->redrawValue();

    // tree proposal ratio denominator
    hr -= proposal_distribution->computeLnProbability();

    Tree proposal_tree = proposal_distribution->getValue();

    // if necessary, root the new tree topology
    if( proposal_tree.isRooted() == false )
    {
        TopologyNode *node;

        bool make_bifurcating = false;
        // get the outgroup node
        if( outgroup.size() > 0 )
        {
            proposal_tree.reroot( outgroup, make_bifurcating, true );
            node = proposal_tree.getRoot().getNode( outgroup, true );
        }
        // or root the tree randomly
        else
        {
            node = &proposal_tree.getNode( size_t( rng->uniform01() * (proposal_tree.getNumberOfNodes() - 1) ) );
            proposal_tree.reroot(*node, make_bifurcating, true);
        }

        // get the old root
        TopologyNode& parent = node->getParent();

        // create a new root
        TopologyNode* root = new TopologyNode( );

        // attach the root between the node and its parent
        parent.removeChild(node);
        parent.addChild( root );
        root->setParent(&parent);

        root->addChild(node);
        node->setParent(root);

        // randomly split the branch length descending from the root
        double brlen = proposal_tree.getBranchLengthForNode(*node);
        double u = rng->uniform01();

        node->setBranchLength( brlen * u );
        root->setBranchLength( brlen * (1.0 - u) );

        // jacobian of root branch transformation
        hr += brlen;

        // reroot the tree
        proposal_tree.reroot(*node, make_bifurcating, true);
        proposal_tree.setRooted(true);

        // restore tip indices
        for (size_t i=0; i<proposal_tree.getNumberOfTips(); i++)
        {
            proposal_tree.getTipNodeWithName(variable->getValue().getTipNode(i).getName()).setIndex(i);
        }
        proposal_tree.orderNodesByIndex();
    }

    // first get the proposed branch lengths
    std::vector<double> branch_lengths;
    const std::vector<TopologyNode*>& proposal_nodes = proposal_tree.getNodes();
    for(size_t i = 0; i < proposal_nodes.size(); i++)
    {
        if( proposal_nodes[i]->isRoot() == false )
        {
            branch_lengths.push_back( proposal_tree.getBranchLengthForNode(*proposal_nodes[i]) );
        }
    }

    // set the tip ages and sampled ancestor flags
    for(size_t i = 0; i < proposal_tree.getNumberOfTips(); i++)
    {
        proposal_nodes[i]->setAge( variable->getValue().getTipNode(i).getAge() );

        if( variable->getValue().getTipNode(i).isSampledAncestorTip() )
        {
            proposal_nodes[i]->getParent().setAge( proposal_nodes[i]->getAge() );
            proposal_nodes[i]->setSampledAncestor( true );
        }
    }

    // get the sorted node ages
    const std::vector<TopologyNode*>& stored_nodes = variable->getValue().getNodes();
    std::vector<double> ages;
    for(size_t i = 0; i < stored_nodes.size(); i++)
    {
        if( stored_nodes[i]->isInternal() == true && stored_nodes[i]->isSampledAncestorTip() == false )
        {
            ages.push_back( stored_nodes[i]->getAge() );
        }
    }
    std::sort(ages.begin(), ages.end());


    // choose a new ranking uniformly at random
    proposal_ln_num_rankings = 0;
    std::vector<size_t> ranking = recursivelyRank( proposal_tree.getRoot(), proposal_ln_num_rankings );

    // ranking proposal ratio
    hr += proposal_ln_num_rankings - stored_ln_num_rankings;


    // update the ranked node ages
    for(size_t i = 0; i < ranking.size(); i++)
    {
        proposal_nodes[ranking[i]]->setAge( ages[i] );
    }

    variable->getValue() = proposal_tree;

    // update the substitution rates, if applicable
    if( rates_node != NULL )
    {
        stored_rates = rates_node->getValue();

        std::vector<double> rates( stored_rates.size(), 0.0);

        for(size_t i = 0; i < stored_rates.size(); i++)
        {
            rates[i] = branch_lengths[i] / proposal_tree.getBranchLengthForNode(i);
        }

        rates_node->getValue() = rates;
    }
    else if( rates_vector.empty() == false )
    {
        stored_rates.clear();

        for(size_t i = 0; i < rates_vector.size(); i++)
        {
            stored_rates.push_back(rates_vector[i]->getValue());

            rates_vector[i]->getValue() = branch_lengths[i] / proposal_tree.getBranchLengthForNode(i);
        }
    }

    // use the first proposal to initialize the values
    if( init == false )
    {
        hr = RbConstants::Double::inf;
        init = true;
    }

    return  hr;
}


/** Prepare the proposal */
void IndependentTopologyProposal::prepareProposal( void )
{
    
    // nothing to do
}



/** Print the parameter summary */
void IndependentTopologyProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // nothing to do
}


/** Recursively choose a random subtree node ranking */
std::vector<size_t> IndependentTopologyProposal::recursivelyRank( const TopologyNode& n, double& ln_num_rankings )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    std::vector<size_t> ranking;

    for(size_t i = 0; i < n.getNumberOfChildren(); i++)
    {
        if( n.getChild(i).isTip() )
        {
            continue;
        }

        // get the child ranking
        std::vector<size_t> child_ranking = recursivelyRank( n.getChild(i), ln_num_rankings );

        if( child_ranking.empty() == true )
        {
            continue;
        }

        // choose new ranks for the child nodes
        std::vector<size_t> pos;
        for(size_t j = 0; j < child_ranking.size(); j++)
        {
            size_t p = rng->uniform01()*(ranking.size() + 1.0);
            pos.push_back(p);
        }

        std::sort(pos.begin(), pos.end());

        // update the ranking
        for(size_t j = 0; j < pos.size() ; j++)
        {
            ranking.insert(ranking.begin() + pos[j], child_ranking[j]);

            for(size_t k = j+1; k < pos.size(); k++)
            {
                pos[k]++;
            }
        }

        // update the number of rankings
        ln_num_rankings += RbMath::lnChoose(ranking.size(), child_ranking.size());
    }

    if( n.isSampledAncestorTipOrParent() == false )
    {
        // put the current node at the highest rank
        ranking.push_back( n.getIndex() );
    }

    return ranking;
}



void IndependentTopologyProposal::undoProposal( void )
{
    
    // reset to the old tree
    variable->setValue( stored_tree.clone() );
    
    if( rates_node != NULL )
    {
        rates_node->getValue() = stored_rates;
    }
    else if( rates_vector.empty() == false )
    {
        for(size_t i = 0; i < rates_vector.size(); i++)
        {
            rates_vector[i]->getValue() = stored_rates[i];
        }
    }
}


void IndependentTopologyProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/** Tune the proposal */
void IndependentTopologyProposal::tune(double r)
{
    
    // nothing to do
}


void IndependentTopologyProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if( oldN == variable )
    {
        variable = static_cast<StochasticNode<Tree>* >( newN );
    }
    else if( oldN == rates_node )
    {
        rates_node = static_cast<StochasticNode<RbVector<double> >* >( newN );
    }
    else if( rates_vector.empty() == false )
    {
        for(size_t i = 0; i < rates_vector.size(); i++)
        {
            if(oldN == rates_vector[i])
            {
                rates_vector[i] = static_cast<StochasticNode<double>* >( newN );
            }
        }
    }
    
}

