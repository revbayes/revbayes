#include "TreeScaleProposal.h"

#include <cstddef>
#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "TreeUtilities.h"
#include "Cloneable.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
TreeScaleProposal::TreeScaleProposal( StochasticNode<Tree> *t, StochasticNode< RbVector<Tree> > *vec_n, StochasticNode<double> *r, double d ) : Proposal(),
    tree( t ),
    vector_variable( vec_n),
    rootAge( r ),
    delta( d )
{
    // tell the base class to add the node
    addNode( tree );
    addNode( vector_variable );
    addNode( rootAge );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void TreeScaleProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
TreeScaleProposal* TreeScaleProposal::clone( void ) const
{
    
    return new TreeScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& TreeScaleProposal::getProposalName( void ) const
{
    static std::string name = "TreeScale";
    
    return name;
}


double TreeScaleProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


namespace {

// Follow rescaleSubtree's binary traversal: tips and sampled-ancestor parents stop scaling
// of their entire subtrees. Only the internal ages above these boundaries enter the Jacobian.
size_t countScaledAges(const TopologyNode& node)
{
    if (node.isTip() or node.isSampledAncestorParent())
        return 0;

    return 1 + countScaledAges(node.getChild(0)) + countScaledAges(node.getChild(1));
}

}

// Scale eligible internal ages together, keeping an optional root-age scalar synchronized.
double TreeScaleProposal::doProposal( void )
{
    
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
        
    // get the current tree either from the single variable or the vector of trees
    Tree *tmp = NULL;
    if ( tree != NULL )
    {
        tmp = &tree->getValue();
    }
    else
    {
        tree_index = floor(rng->uniform01() * vector_variable->getValue().size());
        tmp = &(vector_variable->getValue()[tree_index]);
    }
    Tree& tau = *tmp;
    
    TopologyNode& node = tau.getRoot();
    
    // we need to work with the times
    double my_age      = node.getAge();

    // rootAge represents the root coordinate, not an independent origin or stem age.
    // Check before mutation: undo must not overwrite a different scalar with the root age.
    if (rootAge != NULL)
    {
        if (rootAge->getValue() != my_age)
            throw RbException("mvTreeScale: rootAge must equal the tree's root age; do not pass an origin age as rootAge.");
        storedRootAge = rootAge->getValue();
    }
    
    // now we store all necessary values
    storedAges = std::vector<double>(tau.getNumberOfNodes(), 0.0);
    TreeUtilities::getAges(node, storedAges);
    const size_t num_scaled_ages = countScaledAges(node);
    
    // draw new ages 
    double u = rng->uniform01();
    double scaling_factor = std::exp( delta * ( u - 0.5 ) );
    
    // rescale the subtrees
    TreeUtilities::rescaleSubtree(node, scaling_factor );
    
    if ( rootAge != NULL )
    {
        // A sampled-ancestor root is a fixed boundary, so use its actual age after scaling.
        rootAge->setValue( new double(node.getAge()) );
    }
    
    // The log-symmetric scale draw gives a Hastings factor s^k for k scaled coordinates.
    // A root-age scalar represents the same coordinate as the tree root, not an extra one.
    // Fixed boundaries and their descendants contribute nothing, including a no-op at the root.
    double lnHastingsratio = log( scaling_factor ) * num_scaled_ages;
    
    return lnHastingsratio;
}


/**
 *
 */
void TreeScaleProposal::prepareProposal( void )
{
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void TreeScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "delta = ";
    if (name_only == false)
    {
        o << delta;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void TreeScaleProposal::undoProposal( void )
{
        
    // get the current tree either from the single variable or the vector of trees
    Tree *tmp = NULL;
    if ( tree != NULL )
    {
        tmp = &tree->getValue();
    }
    else
    {
        tmp = &(vector_variable->getValue()[tree_index]);
    }
    Tree& tau = *tmp;
    
    TopologyNode& node = tau.getRoot();
    
    // undo the proposal
    TreeUtilities::setAges(node, storedAges );
    
    if ( rootAge != NULL )
    {
        rootAge->setValue( new double(storedRootAge) );
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void TreeScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == tree )
    {
        tree = static_cast<StochasticNode<Tree>* >(newN);
    }
    else if ( oldN == vector_variable )
    {
        vector_variable = static_cast<StochasticNode< RbVector<Tree> >* >(newN);
    }
    else if ( oldN == rootAge )
    {
        rootAge = static_cast<StochasticNode<double>* >(newN);
    }
    
}


void TreeScaleProposal::setProposalTuningParameter(double tp)
{
    delta = tp;
}


/**
 * Tune the Proposal to accept at the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.234.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void TreeScaleProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        delta *= (1.0 + ((rate - p)/(1.0 - p)) );
    }
    else
    {
        delta /= (2.0 - rate/p);
    }
    
}

