#include "ExtinctionReversibleJumpProposal.h"

#include <cstddef>
#include <string>
#include <vector>

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "TreeDistributionProperties.h"

using namespace RevBayesCore;


ExtinctionReversibleJumpProposal::ExtinctionReversibleJumpProposal( StochasticNode<MatrixReal> *n ) : Proposal(),
    matrix( n )
{
    addNode( matrix );

    if ( dynamic_cast<AbstractFossilizedBirthDeathRangeProcess *>( &matrix->getDistribution() ) == NULL )
    {
        throw RbException("mvExtinctionRJSwitch requires a fossilized birth-death range process.");
    }

    collectEligible();
}


ExtinctionReversibleJumpProposal::ExtinctionReversibleJumpProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n )
{
    addNode( tree );

    if ( dynamic_cast<AbstractFossilizedBirthDeathRangeProcess *>( &tree->getDistribution() ) == NULL )
    {
        throw RbException("mvExtinctionRJSwitch requires a fossilized birth-death range process.");
    }

    // only an extended tree has extinction times as tip ages, so only it has a point mass
    const TreeDistributionProperties *props = dynamic_cast<const TreeDistributionProperties *>( &tree->getDistribution() );
    if ( props == NULL || props->isExtended() == false )
    {
        throw RbException("mvExtinctionRJSwitch requires an extended tree, whose tips are extinction times.");
    }

    collectEligible();
}


/** A taxon reported extant is pinned at the present and has no second state to jump to. */
void ExtinctionReversibleJumpProposal::collectEligible( void )
{
    // fixed for the run, so the 1/n of picking one cancels between the two directions
    const std::vector<Taxon> &taxa = rangeProcess().getTaxa();
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        if ( taxa[i].isExtinct() == true ) eligible.push_back( i );
    }
}


/** The distribution, resolved on each use so a swapped node cannot leave a stale pointer behind. */
AbstractFossilizedBirthDeathRangeProcess& ExtinctionReversibleJumpProposal::rangeProcess( void ) const
{
    Distribution &d = ( matrix != NULL ? (Distribution&)matrix->getDistribution()
                                       : (Distribution&)tree->getDistribution() );

    return dynamic_cast<AbstractFossilizedBirthDeathRangeProcess &>( d );
}


void ExtinctionReversibleJumpProposal::cleanProposal( void )
{
    ; // do nothing
}


ExtinctionReversibleJumpProposal* ExtinctionReversibleJumpProposal::clone( void ) const
{
    return new ExtinctionReversibleJumpProposal( *this );
}


const std::string& ExtinctionReversibleJumpProposal::getProposalName( void ) const
{
    static std::string name = "ExtinctionReversibleJump";

    return name;
}


double ExtinctionReversibleJumpProposal::getProposalTuningParameter( void ) const
{
    // the support is fixed by the data, so there is nothing to tune
    return RbConstants::Double::nan;
}


/**
 * Jump one taxon's extinction time between the present and a time drawn uniformly below its
 * youngest appearance.
 *
 * One direction is deterministic and the other draws from that uniform, so the ratio is the width
 * of the support, up when leaving the present and down when returning to it. The width is read
 * from the youngest appearance, which this move does not touch, so the forward and reverse
 * proposals see the same one.
 *
 * \return The hastings ratio.
 */
double ExtinctionReversibleJumpProposal::doProposal( void )
{
    if ( eligible.empty() == true ) return 0.0;

    RandomNumberGenerator* rng = GLOBAL_RNG;

    stored_index = eligible[ (size_t)( rng->uniform01() * eligible.size() ) ];

    double present = rangeProcess().getPresent();
    // tau_K, the youngest appearance. A tree does not carry it, so read it from the table there
    double last    = ( matrix != NULL ? matrix->getValue()[stored_index][1]
                                      : rangeProcess().getLastAppearance( stored_index ) );
    double width   = last - present;

    if ( width <= 0.0 ) return RbConstants::Double::neginf;

    // the latent status is the model index and the death time follows it, never the reverse
    bool survived = rangeProcess().hasSurvived( stored_index );
    double d = survived ? present + rng->uniform01() * width : present;

    rangeProcess().switchSurvived( stored_index );

    if ( matrix != NULL )
    {
        MatrixReal& v = matrix->getValue();
        stored_death = v[stored_index][0];
        v[stored_index][0] = d;

        matrix->addTouchedElementIndex( stored_index * v.getNumberOfColumns() );
    }
    else
    {
        // on an extended tree the tip age is the extinction time
        TopologyNode& tip = tree->getValue().getNode( stored_index );
        stored_death = tip.getAge();
        tip.setAge( d );
    }

    return survived ? log( width ) : -log( width );
}


void ExtinctionReversibleJumpProposal::prepareProposal( void )
{

}


void ExtinctionReversibleJumpProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

}


void ExtinctionReversibleJumpProposal::undoProposal( void )
{
    // the distribution restores the status from restoreSpecialization, which the DAG calls next
    if ( matrix != NULL )
    {
        MatrixReal& v = matrix->getValue();

        v[stored_index][0] = stored_death;

        matrix->addTouchedElementIndex( stored_index * v.getNumberOfColumns() );
    }
    else
    {
        tree->getValue().getNode( stored_index ).setAge( stored_death );
    }
}


void ExtinctionReversibleJumpProposal::setProposalTuningParameter(double tp)
{

}


void ExtinctionReversibleJumpProposal::tune( double rate )
{

}


void ExtinctionReversibleJumpProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( matrix != NULL ) matrix = static_cast<StochasticNode<MatrixReal>* >(newN);
    else                  tree   = static_cast<StochasticNode<Tree>* >(newN);
}
