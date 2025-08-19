#include "GibbsDrawCharacterHistoryProposal.h"

#include <cstddef>
#include <vector>

#include "RbException.h"
#include "StateDependentSpeciationExtinctionProcess.h"
#include "Tree.h"
#include "Cloneable.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
GibbsDrawCharacterHistoryProposal::GibbsDrawCharacterHistoryProposal( StochasticNode<Tree> *n ) : Proposal(),
    variable( n )
{
    // tell the base class to add the node
    addNode( variable );
    
    distribution = dynamic_cast< StateDependentSpeciationExtinctionProcess* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("The GibbsDrawCharacterHistoryProposal is currently only implemented for CDBDP distributions.");
    }
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void GibbsDrawCharacterHistoryProposal::cleanProposal( void )
{
    
    
}


GibbsDrawCharacterHistoryProposal* GibbsDrawCharacterHistoryProposal::clone( void ) const
{

        return new GibbsDrawCharacterHistoryProposal( *this );
}

/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double GibbsDrawCharacterHistoryProposal::doProposal( void )
{

    size_t num_nodes = distribution->getValue().getNumberOfNodes();
    std::vector<std::string> character_histories(num_nodes);
    distribution->drawStochasticCharacterMap(character_histories);
    
    return 0.0;
}


const std::string& GibbsDrawCharacterHistoryProposal::getProposalName( void ) const
{
    static std::string name = "GibbsDrawCharacterHistoryProposal";
    return name;
}


double GibbsDrawCharacterHistoryProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


void GibbsDrawCharacterHistoryProposal::prepareProposal( void )
{
    
}


void GibbsDrawCharacterHistoryProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    
}


void GibbsDrawCharacterHistoryProposal::undoProposal( void )
{
    
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void GibbsDrawCharacterHistoryProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<Tree>* >(newN) ;
    
    distribution = dynamic_cast< StateDependentSpeciationExtinctionProcess* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("The GibbsDrawCharacterHistoryProposal is currently only implemented for CDBDP distributions.");
    }
}


void GibbsDrawCharacterHistoryProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


void GibbsDrawCharacterHistoryProposal::tune( double rate )
{
    
    
}

