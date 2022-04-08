#include "ElementsSwapSimplexProposal.h"

#include <stddef.h>
#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
ElementsSwapSimplexProposal::ElementsSwapSimplexProposal( StochasticNode<Simplex> *n ) : Proposal(),
simplex ( n )
{
    // tell the base class to add the node
    addNode( simplex );
    
}


void ElementsSwapSimplexProposal::addIndex( size_t v )
{
    if (v < simplex->getValue().size())
    {
        indices.insert(v);
    }
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void ElementsSwapSimplexProposal::cleanProposal( void )
{
    RbOrderedSet<DagNode*> affected;
    simplex->initiatefindUniqueDescendants( affected );
    
    for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
    {
        if ((*it)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC && (*it)->getNumberOfElements() == stored_simplex.size())
        {
            (*it)->clearTouchedElementIndices();
        }
    }
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
ElementsSwapSimplexProposal* ElementsSwapSimplexProposal::clone( void ) const
{
    return new ElementsSwapSimplexProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& ElementsSwapSimplexProposal::getProposalName( void ) const
{
    static std::string name = "ElementsSwapSimplex";
    
    return name;
}


double ElementsSwapSimplexProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Beta(val[index] * lambda)
 * where lambda is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double ElementsSwapSimplexProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // now we store all necessary values
    stored_simplex = simplex->getValue();
    const RbVector<double>& simplex_current = simplex->getValue();
    RbVector<double> simplex_new = simplex_current;
    
    // randomly draw two indices
    size_t chosen_index1;
    size_t chosen_index2;
    if ( indices.size() >= 2 )
    {
        size_t index1 = size_t( floor(rng->uniform01() * double(indices.size())) );
        size_t index2 = size_t( floor(rng->uniform01() * double(indices.size())) );
        while (index1 == index2)
        {
            index2 = size_t( floor(rng->uniform01() * double(indices.size())) );
        }
        
        std::set<size_t>::iterator it1 = indices.begin();
        std::set<size_t>::iterator it2 = indices.begin();
        std::advance(it1, index1);
        std::advance(it2, index2);
        
        chosen_index1 = *it1;
        chosen_index2 = *it2;
    }
    else
    {
        chosen_index1 = size_t( floor(rng->uniform01() * double(simplex_current.size())) );
        chosen_index2 = size_t( floor(rng->uniform01() * double(simplex_current.size())) );
        while (chosen_index1 == chosen_index2)
        {
            chosen_index2 = size_t( floor(rng->uniform01() * double(simplex_current.size())) );
        }
    }
    
    // swap the values
    simplex_new[chosen_index2] = simplex_current[chosen_index1];
    simplex_new[chosen_index1] = simplex_current[chosen_index2];

    simplex->setValue( new Simplex(simplex_new), false);
    
    RbOrderedSet<DagNode*> affected;
    simplex->initiatefindUniqueDescendants( affected );
    
    for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
    {
        if ((*it)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC && (*it)->getNumberOfElements() == simplex_new.size())
        {
            (*it)->addTouchedElementIndex(chosen_index1);
            (*it)->addTouchedElementIndex(chosen_index2);
        }
    }

    return 0.0;
    
}


/**
 *
 */
void ElementsSwapSimplexProposal::prepareProposal( void )
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
void ElementsSwapSimplexProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void ElementsSwapSimplexProposal::undoProposal( void )
{
    
    // undo the proposal
    simplex->setValue( new Simplex(stored_simplex), false);
    
    RbOrderedSet<DagNode*> affected;
    simplex->initiatefindUniqueDescendants( affected );
    
    for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
    {
        if ((*it)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC && (*it)->getNumberOfElements() == stored_simplex.size())
        {
            (*it)->clearTouchedElementIndices();
        }
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void ElementsSwapSimplexProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    simplex = static_cast<StochasticNode<Simplex>* >(newN) ;
}


void ElementsSwapSimplexProposal::setProposalTuningParameter(double tp)
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
void ElementsSwapSimplexProposal::tune( double rate )
{
    
}


