#include "RandomCategoryWalkProposal.h"


#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
RandomCategoryWalkProposal::RandomCategoryWalkProposal( StochasticNode< RbVector<long> >* n) : Proposal(),
variable( n )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void RandomCategoryWalkProposal::cleanProposal( void )
{
    ; //variable->clearTouchedElementIndices();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RandomCategoryWalkProposal* RandomCategoryWalkProposal::clone( void ) const
{
    
    return new RandomCategoryWalkProposal( *this );
}

/**
 * Perform the proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and Slides the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double RandomCategoryWalkProposal::doProposal( void )
{
    
    // reset the failed flag
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // get values
    RbVector<long> &val = variable->getValue();
    
    // how many elements did we have in the vector
    size_t num_elements = val.size();
    
    // choose one value in the vector
    size_t index = (size_t)(rng->uniform01() * num_elements);
    
    // store the index
    chosen_index = index;
    
    // randomly choose if we move up or down
    double u = rng->uniform01();
    
    // check if this resulted in a valid move
    if ( (u < 0.5 && index == 0) || ( u > 0.5 && index == num_elements-1) || (val[index] == 0) )
    {
        failed = true;
        
        // return a hastings ratio such that the move will not be accepted
        return RbConstants::Double::neginf;
    }
    if ( u < 0.5 )
    {
        chosen_neighbour = index - 1;
    }
    else
    {
        chosen_neighbour = index + 1;
    }
    
    // now update the elements
    --val[index];
    ++val[chosen_neighbour];
    
    // return the Hastings ratio of this symmetric proposal
    return 0.0;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RandomCategoryWalkProposal::getProposalName( void ) const
{
    static std::string name = "RandomCategoryWalk";
    
    return name;
}


double RandomCategoryWalkProposal::getProposalTuningParameter( void ) const
{
    // dummy return value
    return 0.0;
}


/**
 * Prepare the proposal, e.g., pick the element that we want to change.
 * Here we do not need to do any preparation.
 */
void RandomCategoryWalkProposal::prepareProposal( void )
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
void RandomCategoryWalkProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // nothing to print
    
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void RandomCategoryWalkProposal::undoProposal( void )
{
    
    if ( failed == false )
    {
        // get values
        RbVector<long> &val = variable->getValue();
        
        // now update the elements
        ++val[chosen_index];
        --val[chosen_neighbour];        
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void RandomCategoryWalkProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( variable == oldN )
    {
        variable = static_cast<StochasticNode< RbVector<long> > *>(newN);
    }
}


void RandomCategoryWalkProposal::setProposalTuningParameter(double tp)
{
    // nothing to be done
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * There is no tuning possible
 */
void RandomCategoryWalkProposal::tune( double rate )
{
    // nothing to be done
}

