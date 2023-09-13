#include <cstddef>
#include <iostream>
#include <string>

#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "HomeologPhaseProposal.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class DnaState; }




/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
RevBayesCore::HomeologPhaseProposal::HomeologPhaseProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *c, std::string t1, std::string t2 ) : Proposal(),
    ctmc(c),
    tip1(t1),
    tip2(t2)
{

    addNode( ctmc );


}



/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RevBayesCore::HomeologPhaseProposal* RevBayesCore::HomeologPhaseProposal::clone( void ) const
{
    return new HomeologPhaseProposal( *this );
}


void RevBayesCore::HomeologPhaseProposal::cleanProposal(void)
{
}


/**
 * Perform the Proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double RevBayesCore::HomeologPhaseProposal::doProposal( void )
{
    ctmc->getValue().switchHomeologPhase(tip1, tip2);
    AbstractPhyloCTMCSiteHomogeneous<DnaState> *ctmc_dist = NULL;
    ctmc_dist = static_cast<AbstractPhyloCTMCSiteHomogeneous<DnaState>* >( &ctmc->getDistribution() );
    ctmc_dist->swap_taxon_name_2_tip_index(tip1, tip2);
    return 0.0;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RevBayesCore::HomeologPhaseProposal::getProposalName( void ) const
{
    static std::string name = "HomeologPhaseProposal";

    return name;
}


double RevBayesCore::HomeologPhaseProposal::getProposalTuningParameter(void) const
{
    return 0;
}


void RevBayesCore::HomeologPhaseProposal::prepareProposal(void)
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
void RevBayesCore::HomeologPhaseProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void RevBayesCore::HomeologPhaseProposal::setProposalTuningParameter(double tp)
{
}


/**
 * Swap the current ctmc for a new one.
 *
 * \param[in]     oldN     The old ctmc that needs to be replaced.
 * \param[in]     newN     The new ctmc.
 */
void RevBayesCore::HomeologPhaseProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    if (oldN == ctmc)
    {
        ctmc = static_cast<StochasticNode<AbstractHomologousDiscreteCharacterData>* >(newN) ;
    }
}


void RevBayesCore::HomeologPhaseProposal::tune(double r)
{
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the ctmc/DAG-node to its original value.
 */
void RevBayesCore::HomeologPhaseProposal::undoProposal( void )
{
    ctmc->getValue().switchHomeologPhase(tip1, tip2);
    AbstractPhyloCTMCSiteHomogeneous<DnaState> *ctmc_dist = NULL;
    ctmc_dist = static_cast<AbstractPhyloCTMCSiteHomogeneous<DnaState>* >( &ctmc->getDistribution() );
    ctmc_dist->swap_taxon_name_2_tip_index(tip1, tip2);
}

