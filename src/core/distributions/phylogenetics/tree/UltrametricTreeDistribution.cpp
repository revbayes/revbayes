#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "TreeUtilities.h"
#include "UltrametricTreeDistribution.h"
#include "RbBitSet.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StringUtilities.h"
#include "Trace.h"
#include "TraceTree.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }


#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;

UltrametricTreeDistribution::UltrametricTreeDistribution( TypedDistribution<Tree>* tp, TypedDistribution<double>* rp, TypedDagNode<double> *ra, TypedDagNode<double> *rbf, const TraceTree &tree_trace, Trace<double>* dens, MEAN m) : TypedDistribution<Tree>( new Tree() ),
    tree_prior( tp ),
    rate_prior( rp ),
    root_age( ra ),
    root_branch_fraction( rbf ),
    sample_prior_density( dens ),
    num_samples( 0 ),
    sample_block_start( 0 ),
    sample_block_end( num_samples ),
    sample_block_size( num_samples ),
    ln_probs(num_samples, 0.0),
    num_taxa( 0 ),
    mean_method( m )
{

    size_t burnin = tree_trace.getBurnin();
    long n_samples = tree_trace.size();
    for (int i=burnin; i<n_samples; ++i)
    {
        trees.push_back( tree_trace.objectAt( i ) );
    }

    if ( trees.size() > 0 )
    {
        num_taxa = trees[0].getNumberOfTips();
    }

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( root_age );
    this->addParameter( root_branch_fraction );


    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars_tp = tree_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars_tp.begin(); it != pars_tp.end(); ++it)
    {
        this->addParameter( *it );
    }

    const std::vector<const DagNode*>& pars_rp = rate_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars_rp.begin(); it != pars_rp.end(); ++it)
    {
        this->addParameter( *it );
    }


    num_samples = trees.size();

    if ( sample_prior_density != NULL && tree_trace.size() != sample_prior_density->size() && tree_trace.getBurnin() != sample_prior_density->getBurnin() )
    {
        throw RbException() << "The size of the tree trace (" << tree_trace.size() << ") and the prior density trace (" << sample_prior_density->size() << ") does not match!";
    }


    sample_block_start = 0;
    sample_block_end   = num_samples;
#ifdef RB_MPI
    sample_block_start = size_t(floor( (double(this->pid - this->active_PID)   / this->num_processes ) * num_samples) );
    sample_block_end   = size_t(floor( (double(this->pid+1 - this->active_PID) / this->num_processes ) * num_samples) );
#endif
    sample_block_size  = sample_block_end - sample_block_start;

#ifdef RB_MPI
    // now we need to populate how is
    pid_per_sample = std::vector<size_t>( num_samples, 0 );
    for (size_t i = 0; i < num_samples; ++i)
    {
        pid_per_sample[i] = size_t(floor( double(i) / num_samples * this->num_processes ) ) + this->active_PID;
    }
#endif

    ln_probs = std::vector<double>(num_samples, 0.0);

    simulateTree();

    prepareTreeSamples( trees );

}


UltrametricTreeDistribution::UltrametricTreeDistribution( const UltrametricTreeDistribution &d ) : TypedDistribution<Tree>( d ),
    tree_prior( d.tree_prior->clone() ),
    rate_prior( d.rate_prior->clone() ),
    root_age( d.root_age ),
    root_branch_fraction( d.root_branch_fraction ),
    trees( d.trees ),
    sample_prior_density( NULL ),
    trees_newick( d.trees_newick ),
    topology_indices( d.topology_indices ),
    tree_branch_lengths( d.tree_branch_lengths ),
    outgroup( d.outgroup ),
    num_samples( d.num_samples ),
    sample_block_start( d.sample_block_start ),
    sample_block_end( d.sample_block_end ),
    sample_block_size( d.sample_block_size ),
    ln_probs( d.ln_probs ),
    num_taxa( d.num_taxa ),
    mean_method( d.mean_method )
{
    if ( d.sample_prior_density != NULL )
    {
        sample_prior_density = d.sample_prior_density->clone();
    }
}


UltrametricTreeDistribution::~UltrametricTreeDistribution()
{

    delete tree_prior;
    delete rate_prior;
    delete sample_prior_density;

}





UltrametricTreeDistribution& UltrametricTreeDistribution::operator=( const UltrametricTreeDistribution &d )
{

    if ( this != &d )
    {
        TypedDistribution<Tree>::operator=( d );

        // free memory
        delete tree_prior;
        delete rate_prior;
        delete sample_prior_density;

        tree_prior              = d.tree_prior->clone();
        rate_prior              = d.rate_prior->clone();
        root_age                = d.root_age;
        root_branch_fraction    = d.root_branch_fraction;
        trees                   = d.trees;
        sample_prior_density    = NULL;
        trees_newick            = d.trees_newick;
        topology_indices        = d.topology_indices;
        tree_branch_lengths     = d.tree_branch_lengths;
        outgroup                = d.outgroup;
        num_samples             = d.num_samples;
        sample_block_start      = d.sample_block_start;
        sample_block_end        = d.sample_block_end;
        sample_block_size       = d.sample_block_size;
        ln_probs                = d.ln_probs;
        num_taxa                = d.num_taxa;
        mean_method             = d.mean_method;

        if ( d.sample_prior_density != NULL )
        {
            sample_prior_density = d.sample_prior_density->clone();
        }

    }

    return *this;
}


/* Clone function */
UltrametricTreeDistribution* UltrametricTreeDistribution::clone( void ) const
{

    return new UltrametricTreeDistribution( *this );
}


Split UltrametricTreeDistribution::collectSplits(const TopologyNode& node, RbBitSet& intaxa, std::vector<Split>& splits) const
{

    std::vector<Split> child_splits;

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

    Split parent_split(taxa);
    splits[node.getIndex()] = parent_split;

    RbBitSet taxa_rev = taxa;
    taxa_rev.flip();
    Split parent_split_rev(taxa_rev);
    splits[node.getIndex()] = parent_split_rev;

    return parent_split;
}


Split UltrametricTreeDistribution::collectTreeSample(const Tree& tree, const TopologyNode& n, RbBitSet& intaxa, std::map<Split, double>& split_branch_lengths)
{
    double bl = tree.getBranchLengthForNode(n);

    std::vector<Split> child_splits;

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

    Split parent_split(taxa);
    split_branch_lengths[parent_split] = bl;

    RbBitSet taxa_rev = taxa;
    taxa_rev.flip();
    Split parent_split_rev(taxa_rev);
    split_branch_lengths[parent_split_rev] = bl;


    return parent_split;
}

void UltrametricTreeDistribution::computeBranchRates(const Tree &my_time_tree, const std::string &newick, const std::vector<Split> &splits, size_t index, std::vector<double> &rates) const
{

    if ( newick != trees_newick[index] )
    {
        return;
    }

    const std::vector<TopologyNode*> &time_nodes = my_time_tree.getNodes();
    for (size_t i=0; i<time_nodes.size(); ++i)
    {
        TopologyNode* the_time_node = time_nodes[i];
        if ( the_time_node->isRoot() == false )
        {

            double branch_time = my_time_tree.getBranchLengthForNode(*the_time_node);

            const std::map<Split, double> &s = tree_branch_lengths[index];
            if ( s.find(splits[i]) == s.end() )
            {
                rates[i] = RbConstants::Double::nan;
            }

            const std::map<Split, double> &branch_lengths = tree_branch_lengths[index];
            const Split &this_split = splits[i];
            std::map<Split, double>::const_iterator it_branch_length = branch_lengths.find( this_split );
            if ( it_branch_length == branch_lengths.end() )
            {
                throw RbException("Problem in ultrametric tree distribution. Couldn't find branch length ...");
            }
            double branch_exp_num_events = it_branch_length->second;
            if ( the_time_node->getParent().isRoot() == true )
            {
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
                    double sum = my_time_tree.getBranchLengthForNode( the_time_node->getParent().getChild(0) ) + my_time_tree.getBranchLengthForNode( the_time_node->getParent().getChild(1) );
                    frac = branch_time / sum;
                }
                branch_exp_num_events *= frac;
            }
            double branch_rate = branch_exp_num_events / branch_time;

            rates[i] = branch_rate;
        }

    }

}

double UltrametricTreeDistribution::computeBranchRateLnProbability(const Tree &my_time_tree, const std::string &newick, const std::vector<Split> &splits, size_t index) const
{

    // we need to check if the "outgroup" is present first
    // by "outgroup" I mean the left subtree of the rooted tree.
//    const TopologyNode &outgroup = my_tree.getRoot().getChild(0);
//    if ( sampled_tree.containsClade(outgroup, true) == false )
//    {
//        return RbConstants::Double::neginf;
//    }

    if ( newick != trees_newick[index] )
    {
        return RbConstants::Double::neginf;
    }

    // initialize the probability
    double ln_prob =  0.0;

    const std::vector<TopologyNode*> &nodes = my_time_tree.getNodes();
    std::vector<double> rates(nodes.size()-1,0.0);
    computeBranchRates(my_time_tree, newick, splits, index, rates);

    // We need to a bool to tell us if we've included one of the root branches already
    bool root_included = false;

    for (size_t i=0; i<nodes.size(); ++i)
    {
        TopologyNode* the_node = nodes[i];
        if ( the_node->isRoot() == false )
        {
            double branch_rate = rates[i];

            // double branch_time = the_node->getBranchLength();
            // double events = branch_rate * branch_time;

            rate_prior->setValue( new double(branch_rate) );

            // Check if we need only one branch rate from the root branches
            if ( the_node->getParent().isRoot() == true )
            {
                if ( root_branch_fraction != NULL )
                {
                    ln_prob += rate_prior->computeLnProbability();
                }
                else
                {
                    if ( root_included == false )
                    {
                        ln_prob += rate_prior->computeLnProbability();
                        root_included = true;
                    }
                }
            }
            else
            {
                ln_prob += rate_prior->computeLnProbability();
            }
        }

        if ( RbMath::isFinite( ln_prob ) == false )
        {
            break;
        }
    }

    return ln_prob;
}



/* Compute probability */
double UltrametricTreeDistribution::computeLnProbability( void )
{

    // create a temporary copy of the this tree
    const Tree &my_time_tree = *value;

    // make the tree non-rooted
    Tree *my_tree_unrooted = my_time_tree.clone();
    my_tree_unrooted->unroot();

    bool make_bifurcating = false;
    bool reindex = true;
    my_tree_unrooted->reroot( outgroup, make_bifurcating, reindex);

    std::string my_tree_newick = my_tree_unrooted->getPlainNewickRepresentation();

    RbBitSet b( my_time_tree.getNumberOfTips(), false );
    std::vector<Split> my_splits = std::vector<Split>(my_time_tree.getNumberOfNodes(), Split(b) );
    collectSplits(my_time_tree.getRoot(), b, my_splits);

    const std::map<std::string, size_t> &my_taxon_bitmap = my_time_tree.getTaxonBitSetMap();
    std::vector<double> probs    = std::vector<double>(num_samples, 0.0);
    std::vector<double> ln_probs = std::vector<double>(num_samples, 0.0);

//    size_t num_observed = 0;

    size_t BURNIN = 0;
    if ( sample_prior_density != NULL )
    {
        BURNIN = sample_prior_density->getBurnin();
    }

    // add the ln-probs for each sample
    for (size_t i = 0; i < num_samples; ++i)
    {
        if ( i >= sample_block_start && i < sample_block_end )
        {
            const std::map<std::string, size_t> &your_taxon_bitmap = trees[i].getTaxonBitSetMap();

            ln_probs[i] = computeBranchRateLnProbability( my_time_tree, my_tree_newick, my_splits, i );

            if ( sample_prior_density != NULL && RbMath::isFinite( ln_probs[i] ) )
            {
                ln_probs[i] -= sample_prior_density->getValues()[i+BURNIN];
            }

        }

    }

    delete my_tree_unrooted;

#ifdef RB_MPI
    for (size_t i = 0; i < num_samples; ++i)
    {

        if ( this->pid == pid_per_sample[i] )
        {

            // send the likelihood from the helpers to the master
            if ( this->process_active == false )
            {
                // send from the workers the log-likelihood to the master
                MPI_Send(&ln_probs[i], 1, MPI_DOUBLE, this->active_PID, 0, MPI_COMM_WORLD);
            }

        }
        // receive the likelihoods from the helpers
        else if ( this->process_active == true )
        {
            MPI_Status status;
            MPI_Recv(&ln_probs[i], 1, MPI_DOUBLE, pid_per_sample[i], 0, MPI_COMM_WORLD, &status);
        }

    }
#endif


    double max = 0;
    // add the ln-probs for each sample
    for (size_t i = 0; i < num_samples; ++i)
    {

#ifdef RB_MPI
        if ( this->process_active == true )
        {
#endif
            if ( i == 0 || max < ln_probs[i] )
            {
                max = ln_probs[i];
            }

#ifdef RB_MPI
        }
#endif

    }


    // Variable declarations and initialization
    double ln_prob = 0.0;
    double prob    = 0;

#ifdef RB_MPI
    if ( this->process_active == true )
    {
#endif

        // now normalize
        for (size_t i = 0; i < num_samples; ++i)
        {
            probs[i] = exp( ln_probs[i] - max );
            if ( mean_method == ARITHMETIC )
            {
                prob += probs[i];
            }
            else if ( mean_method == HARMONIC )
            {
                prob += 1/probs[i];
            }
        }

        if ( mean_method == ARITHMETIC )
        {
            ln_prob = std::log( prob ) + max - std::log( num_samples );
        }
        else if ( mean_method == HARMONIC )
        {
            ln_prob = std::log( num_samples ) + max - std::log( prob );
        }

#ifdef RB_MPI

        for (size_t i=this->active_PID+1; i<this->active_PID+this->num_processes; ++i)
        {
            MPI_Send(&ln_prob, 1, MPI_DOUBLE, int(i), 0, MPI_COMM_WORLD);
        }

    }
    else
    {

        MPI_Status status;
        MPI_Recv(&ln_prob, 1, MPI_DOUBLE, this->active_PID, 0, MPI_COMM_WORLD, &status);

    }
#endif

    // finally we add the tree prior
    tree_prior->setValue( value->clone() );
    ln_prob += tree_prior->computeLnProbability();

    return ln_prob;

}


void UltrametricTreeDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{

    if ( n == "getSampleProbabilities" )
    {

        bool log_transorm = static_cast<const TypedDagNode<Boolean>* >( args[0] )->getValue();

        rv.clear();
        rv.resize(num_samples);
        for (size_t i = 0; i < num_samples; ++i)
        {
            rv[i] = (log_transorm ? ln_probs[i] : exp(ln_probs[i]));
        }
    }
    else if ( n == "getBranchRates" )
    {

        int index = static_cast<const TypedDagNode<int>* >( args[0] )->getValue();

        Tree &this_tree = *value;
        size_t num_branches = this_tree.getNumberOfNodes() - 1;

        rv.clear();
        rv.resize( num_branches );

        // make the tree non-rooted
        Tree *my_tree_unrooted = this_tree.clone();
        my_tree_unrooted->unroot();

        bool make_bifurcating = false;
        bool reindex = true;
        my_tree_unrooted->reroot( outgroup, make_bifurcating, reindex);

        std::string my_tree_newick = my_tree_unrooted->getPlainNewickRepresentation();

        RbBitSet b( this_tree.getNumberOfTips(), false );
        std::vector<Split> my_splits = std::vector<Split>(this_tree.getNumberOfNodes(), Split(b) );
        collectSplits(this_tree.getRoot(), b, my_splits);

        computeBranchRates(this_tree, my_tree_newick, my_splits, index, rv);

    }
    else
    {
        throw RbException() << "An ultrametric tree distribution does not have a member method called '" << n << "'.";
    }

}


void UltrametricTreeDistribution::redrawValue( void )
{
    simulateTree();
}


/** Simulate the tree conditioned on the time of origin */
void UltrametricTreeDistribution::simulateTree( void )
{

    delete value;

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t index = size_t( rng->uniform01() * trees.size() ) ;

    value = trees[index].clone();
    TopologyNode *new_root = new TopologyNode( value->getNumberOfNodes()+1 );

    TopologyNode &old_root = value->getRoot();
    TopologyNode &og = old_root.getChild(2);
    old_root.removeChild(&og);
    new_root->addChild(&old_root);
    new_root->addChild( &og );
    old_root.setParent( new_root );
    og.setParent( new_root );

    double midpoint = value->getBranchLengthForNode(og) / 2.0;
    old_root.setBranchLength( midpoint );
    og.setBranchLength( midpoint );

    value->setRoot( new_root, true);
    value->setRooted( true );

    TreeUtilities::makeUltrametric(*value);
    TreeUtilities::rescaleTree(value->getRoot(), root_age->getValue()/value->getRoot().getAge());

}

void UltrametricTreeDistribution::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{

    if ( affecter == root_age)
    {
        dag_node->initiateGetAffectedNodes( affected );
    }

}

/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void UltrametricTreeDistribution::keepSpecialization(const DagNode *affecter)
{

    if ( affecter == root_age )
    {
        dag_node->keepAffected();
    }

}


void UltrametricTreeDistribution::prepareTreeSamples(const std::vector<Tree> &trees)
{

    std::vector<std::string> tip_names = trees[0].getTipNames();
    std::sort(tip_names.begin(),tip_names.end());

    outgroup = tip_names[0];

    trees_newick = std::vector<std::string>(num_samples,"");
    tree_branch_lengths = std::vector<std::map<Split, double> >(num_samples, std::map<Split, double>() );
    topology_indices.clear();

    std::map<Split, long>       clade_counts;
    std::map<std::string, long> tree_counts;

    size_t burnin = 0;

    for (size_t i = burnin; i < num_samples; ++i)
    {

        Tree tree = trees[i];

        bool make_bifurcating = false;
        bool reindex = true;
        tree.reroot( outgroup, make_bifurcating, reindex);

        std::string newick = tree.getPlainNewickRepresentation();
        trees_newick[i] = newick;

//        if ( topology_indices.find( newick ) == topology_indices.end() )
//        {
//
//        }
        topology_indices[newick].push_back( i );

        tree_counts[newick]++;

        // get the clades for this tree
        RbBitSet b( tree.getNumberOfTips(), false );
        collectTreeSample(tree, tree.getRoot(), b, tree_branch_lengths[i]);
    }

}

/**
 * Restore the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void UltrametricTreeDistribution::restoreSpecialization(const DagNode *affecter)
{

    if ( affecter == root_age )
    {
        value->getNode( value->getRoot().getIndex() ).setAge( root_age->getValue() );
        dag_node->restoreAffected();
    }

}

/** Swap a parameter of the distribution */
void UltrametricTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == root_age )
    {
        root_age = static_cast<const TypedDagNode<double>* >( newP );
    }

    if (oldP == root_branch_fraction )
    {
        root_branch_fraction = static_cast<const TypedDagNode<double>* >( newP );
    }

    try
    {
        rate_prior->swapParameter(oldP, newP);
    }
    catch (RbException &e)
    {

    }

    try
    {
        tree_prior->swapParameter(oldP, newP);
    }
    catch (RbException &e)
    {

    }

}

/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void UltrametricTreeDistribution::touchSpecialization(const DagNode *affecter, bool touchAll)
{

    if ( affecter == root_age )
    {
        value->getNode( value->getRoot().getIndex() ).setAge( root_age->getValue() );
        dag_node->touchAffected();
    }

}
