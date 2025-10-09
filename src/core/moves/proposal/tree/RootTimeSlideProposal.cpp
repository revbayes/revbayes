#include <cmath>
#include <iostream>

#include "DistributionNormal.h"
#include "RootTimeSlideProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
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
RootTimeSlideProposal::RootTimeSlideProposal( StochasticNode<Tree> *n, double d ) : Proposal(),
variable( n ),
delta( d )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void RootTimeSlideProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RootTimeSlideProposal* RootTimeSlideProposal::clone( void ) const
{
    
    return new RootTimeSlideProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RootTimeSlideProposal::getProposalName( void ) const
{
    static std::string name = "RootTimeSlide";
    
    return name;
}


double RootTimeSlideProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 * This proposal expands or contracts the distance between the root node and its oldest child.
 * Working on the log scale, we slide dist_root_oldest according to a variance 1 Bactrian kernel.
 *
 * \return The hastings ratio.
 */
double RootTimeSlideProposal::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    // get the root node
    TopologyNode* node = &tau.getRoot();
    // cannot move the root if it's a SA
    if(node->isSampledAncestorParent()) return RbConstants::Double::neginf;
    
    // we need to work with the times
    double my_age      = node->getAge();
    double child_Age   = node->getChild( 0 ).getAge();
    if ( child_Age < node->getChild( 1 ).getAge())
    {
        child_Age = node->getChild( 1 ).getAge();
    }
    
    // now we store all necessary values
    stored_node = node;
    stored_age  = my_age;
    
    
    double u           = rng->uniform01();
    double my_new_age  = my_age + ( delta * ( u - 0.5 ) );
    
    // set the age
    tau.getNode(node->getIndex()).setAge( my_new_age );
    
    // compute the Hastings ratio
    double lnHastingsratio = 0.0;
    
    // check for a valid age
    if ( my_new_age < child_Age )
    {
        lnHastingsratio = RbConstants::Double::neginf;
    }
    
    return lnHastingsratio;
    
}


/**
 *
 */
void RootTimeSlideProposal::prepareProposal( void )
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
void RootTimeSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void RootTimeSlideProposal::undoProposal( void )
{
    
    // undo the proposal
    stored_node->setAge( stored_age );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void RootTimeSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void RootTimeSlideProposal::setProposalTuningParameter(double tp)
{
    delta = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void RootTimeSlideProposal::tune( double rate )
{
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        delta *= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else
    {
        delta /= (2.0 - rate/p);
    }
    
    delta = fmin(10000, delta);
    
}

