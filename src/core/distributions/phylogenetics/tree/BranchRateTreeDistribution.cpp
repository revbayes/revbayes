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
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"
#include "TreeChangeEventMessage.h"

namespace RevBayesCore { class DagNode; }
using namespace RevBayesCore;

BranchRateTreeDistribution::BranchRateTreeDistribution(const TypedDagNode<Tree>* tt, TypedDistribution<double>* brp, TypedDagNode<double> *rbf) : TypedDistribution<Tree>( new Tree() ),
    branch_rate_prior( brp ),
    time_tree( tt ),
    root_branch_fraction( rbf ),
    num_taxa( tt->getValue().getNumberOfTips() ),
    time_tree_unrooted( NULL ),
    stored_time_tree_unrooted( NULL ),
    touched_branch_length_tree( true ),
    touched_time_tree( true ),
    was_touched_branch_length_tree( false ),
    was_touched_time_tree( false ),
    newick_time_tree( ),
    newick_branch_length_tree( ),
    splits( ),
    split_to_branch_lengths( ),
    stored_newick_time_tree( ),
    stored_newick_branch_length_tree( ),
    stored_splits( ),
    stored_split_to_branch_lengths( )
{

    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = branch_rate_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        const DagNode* the_node = *it;
        this->addParameter( the_node );
    }

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( root_branch_fraction );
    this->addParameter( time_tree );


    simulateTree();

}


BranchRateTreeDistribution::BranchRateTreeDistribution(const BranchRateTreeDistribution &d) : TypedDistribution<Tree>( d ),
    branch_rate_prior( d.branch_rate_prior->clone() ),
    time_tree( d.time_tree ),
    root_branch_fraction( d.root_branch_fraction ),
    num_taxa( d.num_taxa ),
    time_tree_unrooted( NULL ),
    stored_time_tree_unrooted( NULL ),
    touched_branch_length_tree( true ),
    touched_time_tree( true ),
    was_touched_branch_length_tree( false ),
    was_touched_time_tree( false ),
    newick_time_tree( d.newick_time_tree ),
    newick_branch_length_tree( d.newick_branch_length_tree ),
    splits( d.splits ),
    split_to_branch_lengths( d.split_to_branch_lengths ),
    stored_newick_time_tree( d.stored_newick_time_tree ),
    stored_newick_branch_length_tree( d.stored_newick_branch_length_tree ),
    stored_splits( d.stored_splits ),
    stored_split_to_branch_lengths( d.stored_split_to_branch_lengths )
{
    
    
    if ( d.time_tree_unrooted != NULL ) time_tree_unrooted = d.time_tree_unrooted->clone();
    if ( d.stored_time_tree_unrooted != NULL ) stored_time_tree_unrooted = d.stored_time_tree_unrooted->clone();

    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = branch_rate_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        const DagNode* the_node = *it;
        this->addParameter( the_node );
    }


}



BranchRateTreeDistribution::~BranchRateTreeDistribution()
{

    delete branch_rate_prior;
    
    delete time_tree_unrooted;
    delete stored_time_tree_unrooted;
    
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
            const DagNode* the_node = *it;
            this->addParameter( the_node );
        }
        
        delete time_tree_unrooted;
        delete stored_time_tree_unrooted;
        time_tree_unrooted          = NULL;
        stored_time_tree_unrooted   = NULL;
        if ( d.time_tree_unrooted != NULL ) time_tree_unrooted = d.time_tree_unrooted->clone();
        if ( d.stored_time_tree_unrooted != NULL ) stored_time_tree_unrooted = d.stored_time_tree_unrooted->clone();
        
        touched_branch_length_tree          = true;
        touched_time_tree                   = true;
        was_touched_branch_length_tree      = false;
        was_touched_time_tree               = false;
        newick_time_tree                    = d.newick_time_tree;
        newick_branch_length_tree           = d.newick_branch_length_tree;
        splits                              = d.splits;
        split_to_branch_lengths             = d.split_to_branch_lengths;
        stored_newick_time_tree             = d.stored_newick_time_tree;
        stored_newick_branch_length_tree    = d.stored_newick_branch_length_tree;
        stored_splits                       = d.stored_splits;
        stored_split_to_branch_lengths      = d.stored_split_to_branch_lengths;

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


RbBitSet BranchRateTreeDistribution::collectTreeSample(const Tree& tree, const TopologyNode& n, RbBitSet& intaxa, std::map<RbBitSet, double>& split_branch_lengths)
{
    double bl = tree.getBranchLengthForNode(n);

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

            child_splits.push_back( collectTreeSample(tree, child_node, taxa, split_branch_lengths) );
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
    if ( touched_time_tree == true )
    {
        delete time_tree_unrooted;
        
        time_tree_unrooted = time_tree_copy.clone();
        time_tree_unrooted->unroot();
    }

    // get our branch length tree
    const Tree &branch_length_tree = *value;
    
    // we need to reroot the timetree
    if ( touched_time_tree == true || touched_branch_length_tree == true )
    {
        
        Clade outgroup = branch_length_tree.getRoot().getChild(0).getClade();
        
        // check first if the outgroup is contained in the tree
        bool strict = true;
        bool contains = time_tree_unrooted->getRoot().containsClade(outgroup, strict);

        // if the outgroup is not contained in the current time tree,
        // then the topology must mismatch and thus the probability is 0.0
        if ( contains == false )
        {
            return RbConstants::Double::neginf;
        }

        bool make_bifurcating = false;
        bool reindex = true;
        time_tree_unrooted->reroot(outgroup, make_bifurcating, reindex);
        
    }

    // compare if the time tree and branch length tree topologies match
    const std::map<std::string, size_t> &time_tree_taxon_bitmap = time_tree_unrooted->getTaxonBitSetMap();
    const std::map<std::string, size_t> &branch_length_tree_taxon_bitmap = branch_length_tree.getTaxonBitSetMap();

    // Check that the topologies are identical
    if ( touched_time_tree == true )
    {
        newick_time_tree = time_tree_unrooted->getPlainNewickRepresentation();
    }
    if ( touched_branch_length_tree == true )
    {
        newick_branch_length_tree = branch_length_tree.getPlainNewickRepresentation();
    }
    
    if ( newick_time_tree != newick_branch_length_tree )
    {
        return RbConstants::Double::neginf;
    }

    // compute the branch rates as r = bl / t
    const std::vector<TopologyNode*> &time_tree_nodes = time_tree_copy.getNodes();

    if ( touched_branch_length_tree == true )
    {
        // get the clades for this tree
        RbBitSet b( branch_length_tree.getNumberOfTips(), false );
        split_to_branch_lengths.clear();
        collectTreeSample(branch_length_tree, branch_length_tree.getRoot(), b, split_to_branch_lengths);
    }
    if ( touched_time_tree == true )
    {
        splits.clear();
        splits.resize( time_tree_nodes.size() );
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

            splits[i] = this_split;
        }
    }
    
    for (size_t i=0; i<time_tree_nodes.size(); ++i)
    {

        TopologyNode* the_time_node = time_tree_nodes[i];
        // check if the node is the root node
        if ( the_time_node->isRoot() == true )
        {
            continue;
        }

        const RbBitSet& this_split = splits[i];

        std::map<RbBitSet, double>::const_iterator it_branch_length = split_to_branch_lengths.find( this_split );
        if ( it_branch_length == split_to_branch_lengths.end() )
        {
            throw RbException("Problem in branch rate tree distribution. Couldn't find branch length ...");
        }
        // get branch length
        double branch_exp_num_events = it_branch_length->second;
        // get branch time
        double branch_time = time_tree_copy.getBranchLengthForNode(*the_time_node);

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
                double sum = time_tree_copy.getBranchLengthForNode(the_time_node->getParent().getChild(0)) + time_tree_copy.getBranchLengthForNode(the_time_node->getParent().getChild(1));
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
    
    touched_branch_length_tree  = false;
    touched_time_tree           = false;

    return ln_prob;
}



void BranchRateTreeDistribution::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{

    if ( m == TreeChangeEventMessage::DEFAULT || m == TreeChangeEventMessage::TOPOLOGY )
    {

//        dirty_topology = true;
    }


}


void BranchRateTreeDistribution::keepSpecialization(const DagNode* affecter)
{
    
    if ( was_touched_time_tree == true )
    {
        was_touched_time_tree       = false;
        touched_time_tree           = false;
        
        delete stored_time_tree_unrooted;
        stored_time_tree_unrooted   = NULL;
    }
    
    if ( was_touched_branch_length_tree == true )
    {
        was_touched_branch_length_tree  = false;
        touched_branch_length_tree      = false;
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
    value = time_tree_copy.clone();
    
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
        double branch_time = value->getBranchLengthForNode(*the_node);

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

void BranchRateTreeDistribution::restoreSpecialization(const DagNode *restorer)
{
    
    if ( was_touched_time_tree == true )
    {
        delete time_tree_unrooted;
        was_touched_time_tree       = false;
        touched_time_tree           = false;
        newick_time_tree            = stored_newick_time_tree;
        time_tree_unrooted          = stored_time_tree_unrooted;
        stored_time_tree_unrooted   = NULL;
        splits                      = stored_splits;
    }
    
    if ( was_touched_branch_length_tree == true )
    {
        was_touched_branch_length_tree  = false;
        touched_branch_length_tree      = false;
        newick_branch_length_tree       = stored_newick_branch_length_tree;
        split_to_branch_lengths         = stored_split_to_branch_lengths;
    }
    
}



/** Swap a parameter of the distribution */
void BranchRateTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{

    if (oldP == root_branch_fraction )
    {
        root_branch_fraction = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if ( oldP == time_tree )
    {
        time_tree = static_cast<const TypedDagNode<Tree>* >( newP );        
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


void BranchRateTreeDistribution::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    
    if ( toucher == time_tree )
    {
        touched_time_tree = true;

        if ( was_touched_time_tree == false )
        {
            was_touched_time_tree       = true;
            stored_newick_time_tree     = newick_time_tree;
            stored_time_tree_unrooted   = time_tree_unrooted;
            stored_splits               = splits;
        }
        else
        {
            delete time_tree_unrooted;
        }
        time_tree_unrooted = NULL;
        
    }
    
    if ( toucher == this->dag_node )
    {
        touched_branch_length_tree = true;

        if ( was_touched_branch_length_tree == false )
        {
            was_touched_branch_length_tree      = true;
            stored_newick_branch_length_tree    = newick_branch_length_tree;
            stored_split_to_branch_lengths      = split_to_branch_lengths;
        }
        
    }
    
}
