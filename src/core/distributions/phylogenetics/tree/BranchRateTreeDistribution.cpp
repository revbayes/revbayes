#include "BranchRateTreeDistribution.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <map>
#include <string>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "TreeChangeEventMessage.h"

using namespace RevBayesCore;

BranchRateTreeDistribution::BranchRateTreeDistribution(const TypedDagNode<Tree>* tt, TypedDistribution<double>* brp) : TypedDistribution<Tree>( new Tree() ),
    branch_rate_prior( brp ),
    time_tree( tt ),
    num_taxa( tt->getValue().getNumberOfTips() )
{

    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = branch_rate_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }


    simulateTree();

}


BranchRateTreeDistribution::BranchRateTreeDistribution(const BranchRateTreeDistribution &d) : TypedDistribution<Tree>( d ),
    branch_rate_prior( d.branch_rate_prior->clone() ),
    time_tree( d.time_tree ),
    num_taxa( d.num_taxa )
{

    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = branch_rate_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }


}



BranchRateTreeDistribution::~BranchRateTreeDistribution()
{

    delete branch_rate_prior;
    // the tree will be deleted automatically by the base class

}


BranchRateTreeDistribution& BranchRateTreeDistribution::operator=(const BranchRateTreeDistribution &d)
{

    if ( this != &d )
    {
        TypedDistribution<Tree>::operator=( d );

        // remove the old branch-rate-prior parameters
        const std::vector<const DagNode*>& old_pars = branch_rate_prior->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = old_pars.begin(); it != old_pars.end(); ++it)
        {
            this->removeParameter( *it );
        }
        delete branch_rate_prior;

        branch_rate_prior           = d.branch_rate_prior->clone();
        time_tree                   = d.time_tree;
        num_taxa                    = d.num_taxa;

        // add the parameters of the distribution
        const std::vector<const DagNode*>& pars = branch_rate_prior->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
        {
            this->addParameter( *it );
        }

    }

    return *this;
}


BranchRateTreeDistribution* BranchRateTreeDistribution::clone( void ) const
{

    return new BranchRateTreeDistribution( *this );
}


RbBitSet BranchRateTreeDistribution::collectSplits(const TopologyNode& node, RbBitSet& intaxa, std::vector<RbBitSet>& splits) const
{

    std::vector<RbBitSet> child_splits;

    RbBitSet taxa( num_taxa );

    if ( node.isTip() )
    {
        node.getTaxa(taxa);
    }
    else
    {
        for (size_t i = 0; i < node.getNumberOfChildren(); i++)
        {
            const TopologyNode &child_node = node.getChild(i);

            child_splits.push_back( collectSplits(child_node, taxa, splits) );
        }
    }

    intaxa |= taxa;

    splits[node.getIndex()] = taxa;

    RbBitSet taxa_rev = taxa;
    taxa_rev.flip();
    splits[node.getIndex()] = taxa_rev;

    return taxa;
}


RbBitSet BranchRateTreeDistribution::collectTreeSample(const TopologyNode& n, RbBitSet& intaxa, std::map<RbBitSet, double>& split_branch_lengths)
{
    double bl = n.getBranchLength();

    std::vector<RbBitSet> child_splits;

    RbBitSet taxa( num_taxa );

    if ( n.isTip() )
    {
        n.getTaxa(taxa);
    }
    else
    {
        for (size_t i = 0; i < n.getNumberOfChildren(); i++)
        {
            const TopologyNode &child_node = n.getChild(i);

            child_splits.push_back( collectTreeSample(child_node, taxa, split_branch_lengths) );
        }
    }

    intaxa |= taxa;

    split_branch_lengths[taxa] = bl;

    RbBitSet taxa_rev = taxa;
    taxa_rev.flip();
    split_branch_lengths[taxa_rev] = bl;


    return taxa;
}


double BranchRateTreeDistribution::computeLnProbability( void )
{

    double ln_prob = 0.0;

    // make the time tree unrooted
    const Tree &time_tree_copy = time_tree->getValue();
    Tree *time_tree_unrooted = time_tree_copy.clone();
    time_tree_unrooted->unroot();

    // get our branch length tree
    const Tree &branch_length_tree = *value;

    // compare if the time tree and branch length tree topologies match
    const std::map<std::string, size_t> &time_tree_taxon_bitmap = time_tree_unrooted->getTaxonBitSetMap();
    const std::map<std::string, size_t> &branch_length_tree_taxon_bitmap = branch_length_tree.getTaxonBitSetMap();

    if ( time_tree_taxon_bitmap != branch_length_tree_taxon_bitmap )
    {
        std::cerr << "Ooohhh" << std::endl;
    }
    
    // @Allison: Check that the topologies are identical.
    // If they are not identical, then you need to return -Inf for the probability
    // the cheap way of checking topologies uses the newick strings (without branches)

    // compute the branch rates as r = bl / t
    const std::vector<TopologyNode*> &time_tree_nodes = time_tree_copy.getNodes();
//    const std::vector<TopologyNode*> &branch_length_tree_nodes = branch_length_tree.getNodes();
    
//    std::vector<double> rates(time_tree_nodes.size()-1,0.0);
    

    // get the clades for this tree
    RbBitSet b( branch_length_tree.getNumberOfTips(), false );
    std::map<RbBitSet, double> split_to_branch_lengths;
    collectTreeSample(branch_length_tree.getRoot(), b, split_to_branch_lengths);

    for (size_t i=0; i<time_tree_nodes.size(); ++i)
    {
            
        TopologyNode* the_time_node = time_tree_nodes[i];
        // check if the node is the root node
        if ( the_time_node->isRoot() == true )
        {
            continue;
        }

        RbBitSet this_split = RbBitSet(num_taxa);
        the_time_node->getTaxa( this_split );

        std::map<RbBitSet, double>::const_iterator it_branch_length = split_to_branch_lengths.find( this_split );
        if ( it_branch_length == split_to_branch_lengths.end() )
        {
            throw RbException("Problem in ultrametric tree distribution. Couldn't find branch length ...");
        }
        double branch_exp_num_events = it_branch_length->second;
        
        // check if the node is a descendant of the root
        if ( the_time_node->getParent().isRoot() == true )
        {
            // do something with the root branch, i.e., use a root branch fraction
            // @Allison: Add the root branch fraction variable in here
            branch_exp_num_events *= 0.5;
            
        }

        double branch_time = the_time_node->getBranchLength();

        double branch_rate = branch_exp_num_events / branch_time;

//        rates[i] = branch_rate;

        // compute the probability of each rate
        branch_rate_prior->setValue( new double(branch_rate) );
        ln_prob += branch_rate_prior->computeLnProbability();

    }

    return ln_prob;
}



void BranchRateTreeDistribution::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{

    if ( m == TreeChangeEventMessage::DEFAULT || m == TreeChangeEventMessage::TOPOLOGY )
    {

//        dirty_topology = true;
    }


}

void BranchRateTreeDistribution::redrawValue( void )
{
    simulateTree();
}


void BranchRateTreeDistribution::setValue(RevBayesCore::Tree *v, bool force)
{

    // delegate to super class
    TypedDistribution<Tree>::setValue( v, force );

    // Sebastian: check if anything special needs to be done
}


void BranchRateTreeDistribution::simulateTree( void )
{
    delete value;

    // make our own copy of this time tree as an unrooted tree
    const Tree &time_tree_copy = time_tree->getValue();
    Tree *value = time_tree_copy.clone();
    value->unroot();

    // draw new branch rates
    const std::vector<TopologyNode*> &tree_nodes = value->getNodes();

    for (size_t i=0; i<tree_nodes.size(); ++i)
    {
        TopologyNode* the_node = tree_nodes[i];
        double branch_time = the_node->getBranchLength();

        branch_rate_prior->redrawValue();
        double new_branch_rate = branch_rate_prior->getValue();
        double new_branch_length = new_branch_rate * branch_time;

        the_node->setBranchLength( new_branch_length );
    }

}


/** Swap a parameter of the distribution */
void BranchRateTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{

    if ( branch_rate_prior != NULL )
    {
        branch_rate_prior->swapParameter(oldP,newP);
    }
}
