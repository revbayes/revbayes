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
#include "TypedDagNode.h"
#include "TreeChangeEventMessage.h"

namespace RevBayesCore { class DagNode; }
using namespace RevBayesCore;

BranchRateTreeDistribution::BranchRateTreeDistribution(const TypedDagNode<Tree>* tt, TypedDistribution<double>* brp, TypedDagNode<double> *rbf) : TypedDistribution<Tree>( new Tree() ),
    branch_rate_prior( brp ),
    time_tree( tt ),
    root_branch_fraction( rbf ),
    num_taxa( tt->getValue().getNumberOfTips() )
{

    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = branch_rate_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( root_branch_fraction );


    simulateTree();

}


BranchRateTreeDistribution::BranchRateTreeDistribution(const BranchRateTreeDistribution &d) : TypedDistribution<Tree>( d ),
    branch_rate_prior( d.branch_rate_prior->clone() ),
    time_tree( d.time_tree ),
    root_branch_fraction( d.root_branch_fraction ),
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
        root_branch_fraction        = d.root_branch_fraction;
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
    Tree *branch_length_tree_copy = branch_length_tree.clone();

    // compare if the time tree and branch length tree topologies match
    const std::map<std::string, size_t> &time_tree_taxon_bitmap = time_tree_unrooted->getTaxonBitSetMap();
    const std::map<std::string, size_t> &branch_length_tree_taxon_bitmap = branch_length_tree.getTaxonBitSetMap();

    if ( time_tree_taxon_bitmap != branch_length_tree_taxon_bitmap )
    {
        std::cerr << "Ooohhh" << std::endl;
    }

    // Check that the topologies are identical
    //std::cout << time_tree_unrooted->getPlainNewickRepresentation() << std::endl;
    //std::cout << branch_length_tree_copy->getPlainNewickRepresentation() << std::endl;
    std::string time_tree_newick = time_tree_unrooted->getPlainNewickRepresentation();
    std::string branch_length_tree_newick = branch_length_tree_copy->getPlainNewickRepresentation();

    if ( time_tree_newick != branch_length_tree_newick )
    {
        return RbConstants::Double::neginf;
    }

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
            throw RbException("Problem in branch rate tree distribution. Couldn't find branch length ...");
        }
        // get branch length
        double branch_exp_num_events = it_branch_length->second;
        // get branch time
        double branch_time = the_time_node->getBranchLength();

        // check if the node is a descendant of the root
        if ( the_time_node->getParent().isRoot() == true )
        {

            // there is a Jacobian term because we use the root branch (i.e., the sum of the two branches subtending the root) as a single variable
            // however, the probability density is defined on the two branches subtending the root.
            // the Jacobian is J = root_branch = (left_branch + right_branch)
            ln_prob += log( branch_exp_num_events ) / 2.0; // we need to divide by 2 because we will go twice into this if statement

            // do something with the root branch, i.e., use a root branch fraction
            double frac = 1.0;
            if ( root_branch_fraction != NULL )
            {
                if ( the_time_node == &(the_time_node->getParent().getChild(0)) )
                {
                    frac = root_branch_fraction->getValue();
                }
                else
                {
                    frac = 1.0 - root_branch_fraction->getValue();
                }
            }
            else
            {
                double sum = the_time_node->getParent().getChild(0).getBranchLength() + the_time_node->getParent().getChild(1).getBranchLength();
                frac = branch_time / sum;
            }
            branch_exp_num_events *= frac;
        }


        double branch_rate = branch_exp_num_events / branch_time;

//        rates[i] = branch_rate;

        // compute the probability of each rate
        branch_rate_prior->setValue( new double(branch_rate) );
        ln_prob += branch_rate_prior->computeLnProbability();

    }
    
    delete time_tree_unrooted;
    delete branch_length_tree_copy;

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

    // make our own copy of this time tree
    // our strategy is to loop through all nodes of the time tree, update
    // the branch times with new branch lengths, and then finally unroot
    // the tree to create our new branch length tree
    const Tree &time_tree_copy = time_tree->getValue();
    Tree *value = time_tree_copy.clone();
    
    // we need to change the tree setting so that the branch lengths are used and not the ages.
    // internally, our tree nodes store both an age and branch length variable, but only one should be used.
    value->getRoot().setUseAges( false, true );

    // get time tree nodes
    const std::vector<TopologyNode*> &time_tree_nodes = value->getNodes();

    // loop through time tree nodes and draw new rates to get new branch lengths
    for (size_t i=0; i<time_tree_nodes.size(); ++i)
    {
        TopologyNode* the_node = time_tree_nodes[i];

        // get the branch time
        double branch_time = the_node->getBranchLength();

        // draw new branch rate
        branch_rate_prior->redrawValue();
        double new_branch_rate = branch_rate_prior->getValue();

        // get new branch length
        double new_branch_length = new_branch_rate * branch_time;

        // set new branch length
        the_node->setBranchLength( new_branch_length );
    }

    // now unroot to get the final branch length tree
    value->unroot();

}


/** Swap a parameter of the distribution */
void BranchRateTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{

    if (oldP == root_branch_fraction )
    {
        root_branch_fraction = static_cast<const TypedDagNode<double>* >( newP );
    }

    //if ( branch_rate_prior != NULL )
    try
    {
        branch_rate_prior->swapParameter(oldP,newP);
    }
    catch (RbException &e)
    {

    }
}
