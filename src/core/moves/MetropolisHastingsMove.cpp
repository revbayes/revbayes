#include <cstddef>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "DagNode.h"
#include "DebugMove.h"
#include "MetropolisHastingsMove.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathLogic.h"
#include "AbstractMove.h"
#include "RbOrderedSet.h"
#include "RbException.h"
#include "RbSettings.h"  // for debugMCMC setting
#include <range/v3/all.hpp> // for ranges::views

namespace views = ranges::views;

using namespace RevBayesCore;



/**
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 *
 * \param[in]    w   The weight how often the proposal will be used (per iteration).
 * \param[in]    t   If auto tuning should be used.
 */
MetropolisHastingsMove::MetropolisHastingsMove( Proposal *p, double w, bool t ) : AbstractMove(p->getNodes(), w, t),
    num_accepted_current_period( 0 ),
    num_accepted_total( 0 ),
    proposal( p )
{

    proposal->setMove( this );

}


/**
 * Copy constructor.
 * We need to create a deep copy of the proposal here.
 *
 * \param[in]   m   The object to copy.
 *
 */
MetropolisHastingsMove::MetropolisHastingsMove(const MetropolisHastingsMove &m) : AbstractMove(m),
    num_accepted_current_period( m.num_accepted_current_period ),
    num_accepted_total( m.num_accepted_total ),
    proposal( m.proposal->clone() )
{

    proposal->setMove( this );

}


/**
 * Basic destructor doing nothing.
 */
MetropolisHastingsMove::~MetropolisHastingsMove( void )
{

    delete proposal;
}


/**
 * Overloaded assignment operator.
 * We need a deep copy of the operator.
 */
MetropolisHastingsMove& MetropolisHastingsMove::operator=(const RevBayesCore::MetropolisHastingsMove &m)
{

    if ( this != &m )
    {
        // delegate
        AbstractMove::operator=( m );

        // free memory
        delete proposal;

        num_accepted_current_period     = m.num_accepted_current_period;
        num_accepted_total              = m.num_accepted_total;
        proposal                        = m.proposal->clone();

        proposal->setMove( this );

    }

    return *this;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the MetropolisHastingsMove.
 */
MetropolisHastingsMove* MetropolisHastingsMove::clone( void ) const
{

    return new MetropolisHastingsMove( *this );
}


/**
 * Get moves' name of object
 *
 * \return The moves' name.
 */
const std::string& MetropolisHastingsMove::getMoveName( void ) const
{

    return proposal->getProposalName();
}


double MetropolisHastingsMove::getMoveTuningParameter( void ) const
{
    return proposal->getProposalTuningParameter();
}


/**
 * How often was the move accepted
 */
size_t MetropolisHastingsMove::getNumberAcceptedCurrentPeriod( void ) const
{

    return num_accepted_current_period;
}


/**
 * How often was the move accepted
 */
size_t MetropolisHastingsMove::getNumberAcceptedTotal( void ) const
{

    return num_accepted_total;
}


/**
 * Get the proposal of the move
 *
 * \return The proposal object.
 */
Proposal& MetropolisHastingsMove::getProposal( void )
{

    return *proposal;
}


void MetropolisHastingsMove::performHillClimbingMove( double lHeat, double pHeat )
{

    // Propose a new value
    proposal->prepareProposal();
    double ln_hastings_ratio = proposal->doProposal();


    const RbOrderedSet<DagNode*> &affectedNodes = getAffectedNodes();
    const std::vector<DagNode*> nodes = getDagNodes();

    // first we touch all the nodes
    // that will set the flags for recomputation
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        // get the pointer to the current node
        DagNode* the_node = nodes[i];
        the_node->touch();
    }

    double lnPriorRatio = 0.0;
    double lnLikelihoodRatio = 0.0;


    // compute the probability of the current value for each node
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        // get the pointer to the current node
        DagNode* the_node = nodes[i];

        if ( RbMath::isAComputableNumber(lnPriorRatio) && RbMath::isAComputableNumber(lnLikelihoodRatio) && RbMath::isAComputableNumber(ln_hastings_ratio) )
        {
            if ( the_node->isClamped() )
            {
                lnLikelihoodRatio += the_node->getLnProbabilityRatio();
            }
            else
            {
                lnPriorRatio += the_node->getLnProbabilityRatio();
            }

        }

    }

    // then we recompute the probability for all the affected nodes
    for (RbOrderedSet<DagNode*>::const_iterator it = affectedNodes.begin(); it != affectedNodes.end(); ++it)
    {
        DagNode *the_node = *it;

        if ( RbMath::isAComputableNumber(lnPriorRatio) && RbMath::isAComputableNumber(lnLikelihoodRatio) && RbMath::isAComputableNumber(ln_hastings_ratio) )
        {
            if ( the_node->isClamped() )
            {
                lnLikelihoodRatio += the_node->getLnProbabilityRatio();
            }
            else
            {
                lnPriorRatio += the_node->getLnProbabilityRatio();
            }
        }

    }

    // exponentiate with the chain heat
    double ln_posterior_ratio = pHeat * (lHeat * lnLikelihoodRatio + lnPriorRatio);

    if ( RbMath::isAComputableNumber(ln_posterior_ratio) == false || ln_posterior_ratio < 0.0 )
    {

        proposal->undoProposal();

        // call restore for each node
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            // get the pointer to the current node
            DagNode* the_node = nodes[i];
            the_node->restore();
        }
    }
    else
    {

        num_accepted_total++;
        num_accepted_current_period++;

        // call accept for each node
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            // get the pointer to the current node
            DagNode* the_node = nodes[i];
            the_node->keep();
        }

    }

}

void MetropolisHastingsMove::performMcmcMove( double prHeat, double lHeat, double pHeat )
{
    // These are the nodes we are (directly) modifying.
    const std::vector<DagNode*> nodes = getDagNodes();
    // These are the nodes that are (indirectly) affected.
    const RbOrderedSet<DagNode*> &affected_nodes = getAffectedNodes();

    int logMCMC = RbSettings::userSettings().getLogMCMC();
    int debugMCMC = RbSettings::userSettings().getDebugMCMC();

    // Compute PDFs for nodes and affected nodes if we are going to use them.
    std::map<const DagNode*, double> initialPdfs;
    if (logMCMC >= 3 or debugMCMC >= 1)
	initialPdfs = getNodePrs(nodes, affected_nodes);

    if (logMCMC >= 3)
    {
        std::cerr<<std::setprecision(11);
        for(auto& [node,pr]: initialPdfs)
            std::cerr<<"    BEFORE:   "<<node->getName()<<":  "<<pr<<"\n";
        std::cerr<<"\n";
    }

    /*
     * NOTE: When checking PDFs, ONLY touch/keep (directly modified) nodes.
     *       Do NOT touch/keep (indirectly) affected_nodes.
     *       If we touch/keep affected_nodes, we can hide problems by doing more recalculation.
     */

    if (debugMCMC >= 1)
    {
        // 1. Compute PDFs before proposal, before touch
        auto& untouched_before_proposal = initialPdfs;

        // 2. Touch nodes.
        for (auto node: nodes)
            node->touch();

        // 3. Compute PDFs before proposal, after touch
        auto touched_before_proposal = getNodePrs(nodes, affected_nodes);

        // 4. Keep nodes
        for (auto node: nodes)
            node->keep();

        // 5. Compare pdfs for each node
        compareNodePrs(proposal->getLongProposalName(), untouched_before_proposal, touched_before_proposal, "PDFs not up-to-date before proposal");
    }

    // Propose a new value
    proposal->prepareProposal();
    double ln_hastings_ratio = RbConstants::Double::neginf;
    try {
        ln_hastings_ratio = proposal->doProposal();
    }
    catch (const RbException &e)
    {
        if ( e.getExceptionType() != RbException::MATH_ERROR and e.getExceptionType() != RbException::SKIP_PROPOSAL)
        {
            throw e;
        }
    }

    // Identify nodes that proposal touches
    std::vector<DagNode*> touched_nodes = nodes; //proposal->identifyNodesToTouch();

    // first we touch all the nodes
    // that will set the flags for recomputation
    for (auto node: touched_nodes)
        node->touch();

    bool fail_probability = not std::isfinite(ln_hastings_ratio);

    double ln_prior_ratio = 0.0;
    double ln_likelihood_ratio = 0.0;

    bool zero_or_nan_to_finite = false;

    // compute the probability of the current value for each node
    for (auto node: views::concat(touched_nodes, affected_nodes))
    {
        if (fail_probability) break;

        if (not node->isStochastic()) continue;

        double ratio = 0;
        try {
            // There should be a previous lnProbability because the nodes have been touched.
            double prev = node->getPrevLnProbability();
            // Compute the current lnProbability.
            double current = node->getLnProbability();

            if (std::isfinite(current) and (std::isnan(prev) or prev == RbConstants::Double::neginf))
            {
                // If the log-probability changes from NaN or -Inf to something finite,
                // the ratio would be NaN or +Inf, and the move would be rejected.
                zero_or_nan_to_finite = true;
            }
            else
                ratio = current - prev;
        }
        catch (const RbException &e)
        {
            if ( e.getExceptionType() != RbException::MATH_ERROR )
                throw;

            ratio = RbConstants::Double::neginf;
        }

        if ( node->isClamped() )
            ln_likelihood_ratio += ratio;
        else
            ln_prior_ratio += ratio;

        if ( not std::isfinite(ratio)) fail_probability = true;
    }

    if (logMCMC >= 3)
    {
        for(auto& [node,pr]: getNodePrs(nodes, affected_nodes))
            std::cerr<<"    PROPOSED: "<<node->getName()<<":  "<<pr<<"\n";
        std::cerr<<"\n";
    }

    // exponentiate with the chain heat
    double ln_posterior_ratio = pHeat * (lHeat * ln_likelihood_ratio + prHeat * ln_prior_ratio);

    double ln_acceptance_ratio = ln_posterior_ratio + ln_hastings_ratio;

    bool rejected = false;

    if ( fail_probability )
    {
        // Reject moves where the posterior ratio is -Inf, +Inf, or NaN.
        // (Terms where the previous log-pdf is NaN or -Inf and the current log-pdf is finite are not included.)
        rejected = true;
    }
    else if (zero_or_nan_to_finite)
    {
        // If we get here and any term changed from NaN or -Inf -> finite then accept the move.
    }
    else if (ln_acceptance_ratio >= 0.0)
        ;
    else if (ln_acceptance_ratio < -300.0)
        rejected = true;
    else
    {
        double r = exp(ln_acceptance_ratio);
        // Accept or reject the move
        double u = GLOBAL_RNG->uniform01();
        if (u < r)
            ;
        else
            rejected = true;
    }

    if ( rejected )
    {
        proposal->undoProposal();

        // call restore for each node
        for (auto node: touched_nodes)
            node->restore();
    }
    else
    {
        if ( ln_posterior_ratio < -1000 )
            throw RbException() << "Accepted move '" << proposal->getProposalName() << "' with with posterior ratio of " << ln_posterior_ratio << " and Hastings ratio of " << ln_hastings_ratio << ".";

        num_accepted_total++;
        num_accepted_current_period++;

        // call accept for each node
        for (auto node: touched_nodes)
            node->keep();

        proposal->cleanProposal();
    }

    std::map<const DagNode*, double> finalPdfs;
    if (logMCMC >=3 or (debugMCMC >= 1 and rejected))
	finalPdfs = getNodePrs(nodes, affected_nodes);

    if (logMCMC >= 3)
    {
        for(auto& [node,pr]: finalPdfs)
            std::cerr<<"    FINAL:    "<<node->getName()<<":  "<<pr<<"\n";
        std::cerr<<"\n";
    }

    if (debugMCMC >=1 and rejected)
    {
        compareNodePrs(proposal->getLongProposalName(), initialPdfs, finalPdfs, "PDFs have changed after rejection and restore");
    }

    if (logMCMC >= 2)
    {
        std::cerr<<"    log(posterior_ratio) = "<<ln_posterior_ratio<<"  log(likelihood_ratio) = "<<ln_likelihood_ratio<<"   log(prior_ratio) = "<<ln_prior_ratio<<"\n";
        std::cerr<<"    log(acceptance_ratio) = "<<ln_acceptance_ratio<<"  log(hastings_ratio) = "<<ln_hastings_ratio<<"\n";
        std::cerr<<"  The move was " << (rejected ? "REJECTED." : "ACCEPTED.") << std::endl;
    }

    /*
    // This fixes the problem in #567.
    if (rejected)
    {
        for(auto node: nodes)
            node->touch();
        for(auto node: nodes)
            node->keep();
    }
    */

    /*
     * NOTE: Debug code probably shouldn't call touch/keep here:
     *
     *       Calling touch/keep after reject/restore can hide MCMC problems.
     *       Calling touch/keep after accept should be redundant.
     */
}


/**
 * Print the summary of the move.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void MetropolisHastingsMove::printSummary(std::ostream &o, bool current_period) const
{
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();

    o << std::fixed;
    o << std::setprecision(4);

    // print the name
    const std::string &n = getMoveName();
    size_t spaces = 40 - (n.length() > 40 ? 40 : n.length());
    o << n;
    for (size_t i = 0; i < spaces; ++i)
    {
        o << " ";
    }
    o << " ";

    // print the DagNode name
    const std::string &dn_name = (*nodes.begin())->getName();
    spaces = 20 - (dn_name.length() > 20 ? 20 : dn_name.length());
    o << dn_name;
    for (size_t i = 0; i < spaces; ++i)
    {
        o << " ";
    }
    o << " ";

    // print the weight
    int w_length = 4;
    if (weight > 0) w_length -= (int)log10(weight);
    for (int i = 0; i < w_length; ++i)
    {
        o << " ";
    }
    o << weight;
    o << " ";

    size_t num_tried = num_tried_total;
    size_t num_accepted = num_accepted_total;
    if (current_period == true)
    {
        num_tried = num_tried_current_period;
        num_accepted = num_accepted_current_period;
    }

    // print the number of tries
    int t_length = 9;
    if (num_tried > 0) t_length -= (int)log10(num_tried);
    for (int i = 0; i < t_length; ++i)
    {
        o << " ";
    }
    o << num_tried;
    o << " ";

    // print the number of accepted
    int a_length = 9;
    if (num_accepted > 0) a_length -= (int)log10(num_accepted);

    for (int i = 0; i < a_length; ++i)
    {
        o << " ";
    }
    o << num_accepted;
    o << " ";

    // print the acceptance ratio
    double ratio = num_accepted / (double)num_tried;
    if (num_tried == 0) ratio = 0;
    int r_length = 5;

    for (int i = 0; i < r_length; ++i)
    {
        o << " ";
    }
    o << ratio;
    o << " ";

    proposal->printParameterSummary( o, false );

    o << std::endl;

    o.setf(previousFlags);
    o.precision(previousPrecision);


}

/**
 * Reset the move counters. Here we only reset the counter for the number of accepted moves.
 *
 */
void MetropolisHastingsMove::resetMoveCounters( void )
{
    num_accepted_current_period = 0;
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new variable.
 */
void MetropolisHastingsMove::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    proposal->swapNode(oldN, newN);

}


void MetropolisHastingsMove::setMoveTuningParameter(double tp)
{
    proposal->setProposalTuningParameter(tp);
}


void MetropolisHastingsMove::setNumberAcceptedCurrentPeriod( size_t na )
{
    num_accepted_current_period = na;
}


void MetropolisHastingsMove::setNumberAcceptedTotal( size_t na )
{
    num_accepted_total = na;
}


/**
 * Tune the move to accept the desired acceptance ratio.
 * We only compute the acceptance ratio here and delegate the call to the proposal.
 */
void MetropolisHastingsMove::tune( void )
{

    if ( num_tried_current_period > 2 )
    {
        double rate = num_accepted_current_period / double(num_tried_current_period);

        proposal->tune( rate );
    }

}
