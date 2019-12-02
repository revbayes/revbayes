#include "DistributionUniform.h"
#include "RootTimeSlideUniformProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
RootTimeSlideUniformProposal::RootTimeSlideUniformProposal( StochasticNode<Tree> *n, StochasticNode<double> *o ) : Proposal(),
    variable( n ),
    origin( o )
{
    // tell the base class to add the node
    addNode( variable );
    addNode( origin );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void RootTimeSlideUniformProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RootTimeSlideUniformProposal* RootTimeSlideUniformProposal::clone( void ) const
{
    
    return new RootTimeSlideUniformProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RootTimeSlideUniformProposal::getProposalName( void ) const
{
    static std::string name = "RootTimeSlideUniform";
    
    return name;
}


double RootTimeSlideUniformProposal::getProposalTuningParameter( void ) const
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
double RootTimeSlideUniformProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node = &tau.getRoot();
    
    // we need to work with the times
    double my_age      = node->getAge();
    double child_Age   = node->getChild( 0 ).getAge();
    if ( child_Age < node->getChild( 1 ).getAge())
    {
        child_Age = node->getChild( 1 ).getAge();
    }
    
    // now we store all necessary values
    storedAge = my_age;
    
    // draw new ages and compute the hastings ratio at the same time
    double my_new_age = (origin->getValue() - child_Age) * rng->uniform01() + child_Age;
    
    // set the age
    node->setAge( my_new_age );
    
    return 0.0;
    
}


/**
 *
 */
void RootTimeSlideUniformProposal::prepareProposal( void )
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
void RootTimeSlideUniformProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void RootTimeSlideUniformProposal::undoProposal( void )
{
    
    // undo the proposal
    variable->getValue().getRoot().setAge( storedAge );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void RootTimeSlideUniformProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if (oldN == variable)
    {
        variable = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else if (oldN == origin)
    {
        origin = static_cast<StochasticNode<double>* >(newN) ;
    }
    
}


void RootTimeSlideUniformProposal::setProposalTuningParameter(double tp)
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
void RootTimeSlideUniformProposal::tune( double rate )
{
    
}

