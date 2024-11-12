#include <cmath>
#include <iostream>

#include "DistributionNormal.h"
#include "RootTimeScaleBactrianProposal.h"
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
RootTimeScaleBactrianProposal::RootTimeScaleBactrianProposal( StochasticNode<Tree> *n ) : Proposal(),
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
void RootTimeScaleBactrianProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RootTimeScaleBactrianProposal* RootTimeScaleBactrianProposal::clone( void ) const
{
    
    return new RootTimeScaleBactrianProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RootTimeScaleBactrianProposal::getProposalName( void ) const
{
    static std::string name = "RootTimeScaleBactrian";
    
    return name;
}


double RootTimeScaleBactrianProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * This proposal expands or contracts the distance between the root node and its oldest child.
 * Working on the log scale, we slide dist_root_oldest according to a variance 1 Bactrian kernel.
 *
 * \return The hastings ratio.
 */
double RootTimeScaleBactrianProposal::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    // get the root node
    TopologyNode* node = &tau.getRoot();
    double my_age      = node->getAge();

    // store info for undoProposal before we exit
    storedNode = node;
    storedAge = my_age;

    // cannot move the root if it's a SA
    if(node->isSampledAncestorParent()) return RbConstants::Double::neginf;
    
    // we need to work with the times
    double child_Age   = node->getChild( 0 ).getAge();
    if ( child_Age < node->getChild( 1 ).getAge())
    {
        child_Age = node->getChild( 1 ).getAge();
    }
    
    // the distance between root and child
    double dist_root_oldest = my_age - child_Age;
    
    // See Yang and Rodriguez (2013) SI eqns 19 and 20
    // Currently hard-coding m = 0.95
    double u = rng->uniform01();
    double delta = RbStatistics::Normal::rv(*GLOBAL_RNG) * RbConstants::BACT_SD;
    
    if ( u < 0.5 ) {
        delta += 0.95;
    } else {
        delta -= 0.95;
    }
    
    double scaling_factor = std::exp( delta );

    double my_new_age = child_Age + dist_root_oldest * scaling_factor;
    
    // set the age
    tau.getNode(node->getIndex()).setAge( my_new_age );
    
    // compute the Hastings ratio
    double lnHastingsratio = log( scaling_factor );
    
    return lnHastingsratio;
    
}


/**
 *
 */
void RootTimeScaleBactrianProposal::prepareProposal( void )
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
void RootTimeScaleBactrianProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void RootTimeScaleBactrianProposal::undoProposal( void )
{
    
    // undo the proposal
    storedNode->setAge( storedAge );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void RootTimeScaleBactrianProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void RootTimeScaleBactrianProposal::setProposalTuningParameter(double tp)
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
void RootTimeScaleBactrianProposal::tune( double rate )
{
    
}

