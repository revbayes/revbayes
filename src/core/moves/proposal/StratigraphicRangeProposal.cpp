#include "StratigraphicRangeProposal.h"

#include <cstddef>
#include <string>
#include <vector>

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "FossilizedBirthDeathSpeciationProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"

using namespace RevBayesCore;


StratigraphicRangeProposal::StratigraphicRangeProposal( StochasticNode<Tree> *n ) : Proposal(),
    variable( n )
{
    // tell the base class to add the node
    addNode( variable );

    if ( dynamic_cast<FossilizedBirthDeathSpeciationProcess *>( &variable->getDistribution() ) == NULL )
    {
        throw RbException("mvStratigraphicRange requires a dnFBDSP. A dnFBDRP keeps its appearances in the value, where mvMatrixElementSlide samples them.");
    }

    rangeProcess().setHasResampleMove();
}


/** The distribution, resolved on each use so a swapped node cannot leave a stale pointer behind. */
FossilizedBirthDeathSpeciationProcess& StratigraphicRangeProposal::rangeProcess( void ) const
{
    return static_cast<FossilizedBirthDeathSpeciationProcess &>( variable->getDistribution() );
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void StratigraphicRangeProposal::cleanProposal( void )
{
    ; // do nothing
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
StratigraphicRangeProposal* StratigraphicRangeProposal::clone( void ) const
{
    return new StratigraphicRangeProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& StratigraphicRangeProposal::getProposalName( void ) const
{
    static std::string name = "StratigraphicRange";

    return name;
}


double StratigraphicRangeProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Redraw one taxon's two appearances from the bins that reported them. The support is fixed by the
 * data rather than by the current state, which is what makes this an independence proposal with a
 * Hastings ratio of 1.
 *
 * \return The hastings ratio.
 */
double StratigraphicRangeProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    FossilizedBirthDeathSpeciationProcess &dist = rangeProcess();

    // both bases declare getTaxa; the range table is indexed by the range process's ordering
    size_t i = rng->uniform01() * dist.AbstractFossilizedBirthDeathRangeProcess::getTaxa().size();

    // marks taxon i dirty itself: the tree does not carry these ages, so no touched element index
    // names them and the pull cannot see the write
    dist.resampleFirstLast(i);

    return 0.0;
}


void StratigraphicRangeProposal::prepareProposal( void )
{

}


/**
 * Print the summary of the Proposal.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void StratigraphicRangeProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

}


/**
 * Reject the Proposal.
 */
void StratigraphicRangeProposal::undoProposal( void )
{
    // the distribution restores the two ages from restoreSpecialization, which the DAG calls next
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void StratigraphicRangeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    variable = static_cast<StochasticNode<Tree>* >(newN);
}
