#include <cmath>
#include <iostream>
#include <cstddef>
#include <vector>

#include "DistributionUniform.h"
#include "TipTimeSlideUniformProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "TypedDagNode.h"
#include "Proposal.h"
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
TipTimeSlideUniformProposal::TipTimeSlideUniformProposal( StochasticNode<Tree> *n, TypedDagNode<double> *o ) : Proposal(),
    tree( n ),
    origin( o ),
    use_index( false ),
    failed( false )
{
    // tell the base class to add the node
    addNode( tree );
    addNode( origin );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void TipTimeSlideUniformProposal::cleanProposal( void )
{
    failed = false; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
TipTimeSlideUniformProposal* TipTimeSlideUniformProposal::clone( void ) const
{
    
    return new TipTimeSlideUniformProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& TipTimeSlideUniformProposal::getProposalName( void ) const
{
    static std::string name = "TipTimeSlideUniform";
    
    return name;
}


double TipTimeSlideUniformProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Uniform-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Uniform(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double TipTimeSlideUniformProposal::doProposal( void )
{
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    std::vector<size_t> tips;

    if( use_index == false )
    {
        for (size_t i = 0; i < tau.getNumberOfTips(); ++i)
        {
            TopologyNode* node = &tau.getNode(i);
            if ( node->isFossil() )
            {
                tips.push_back(i);
            }

        }

        if ( tips.empty() )
        {
            failed = true;
            return 0;
        }

        // pick a random fossil node
        double u = rng->uniform01();
        node_index = tips[ size_t( std::floor(tips.size() * u) ) ];
    }
    
    TopologyNode* node = &tau.getNode(node_index);

    TopologyNode& parent = node->getParent();

    // we need to work with the times
    double my_age  = node->getAge();
    double max_age = parent.getAge();
    double min_age = 0.0;

    // adjust min and max age if sampled ancestor
    if (node->isSampledAncestorTip())
    {
        TopologyNode *sibling = &parent.getChild( 0 );
        if ( sibling == node )
        {
            sibling = &parent.getChild( 1 );
        }

        min_age = sibling->getAge();

        if (parent.isRoot())
        {
            if (origin == NULL)
                throw(RbException("Attempting to move root sampled ancestor, but no origin time provided."));

            max_age = origin->getValue();
        }
        else
        {
            TopologyNode& grandParent = parent.getParent();

            max_age = grandParent.getAge();
        }
    }

    // adjust min and max age given taxon data
    Taxon& taxon = node->getTaxon();
    if ( taxon.getMaxAge() < max_age ) {
    	max_age = taxon.getMaxAge();
    }
    if ( taxon.getMinAge() > min_age ) {
    	min_age = taxon.getMinAge();
    }


	// now we store all necessary values
	storedNode = node;
	storedAge = my_age;

	// draw new ages and compute the hastings ratio at the same time
	double my_new_age = min_age + (max_age - min_age) * rng->uniform01();
    
    // set the age
    node->setAge( my_new_age );

    return 0.0;
    
}


/**
 *
 */
void TipTimeSlideUniformProposal::prepareProposal( void )
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
void TipTimeSlideUniformProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void TipTimeSlideUniformProposal::undoProposal( void )
{
    
    // undo the proposal
    if ( failed == false )
    {
        storedNode->setAge( storedAge );
    }
}


/**
 * Set a specific tip index to sample
 */
void TipTimeSlideUniformProposal::useIndex( size_t index )
{
    use_index = true;
    node_index = index;
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void TipTimeSlideUniformProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if (oldN == tree)
    {
        tree = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else if (oldN == origin)
    {
        origin = static_cast<TypedDagNode<double>* >(newN) ;
    }
    
}


void TipTimeSlideUniformProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void TipTimeSlideUniformProposal::tune( double rate )
{
    
}

