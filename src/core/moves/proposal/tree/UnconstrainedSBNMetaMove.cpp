#include "DagNode.h"
#include "UnconstrainedSBNMetaMove.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "RbMathLogic.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace RevBayesCore;


/**
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 *
 * \param[in]    w   The weight how often the proposal will be used (per iteration).
 * \param[in]    t   If auto tuning should be used.
 */
UnconstrainedSBNMetaMove::UnconstrainedSBNMetaMove( UnrootedUnconstrainedSBNProposal* u, Proposal* n, Proposal* s, std::vector<double> m, size_t n0, size_t c0, double w, bool t ) : AbstractMove(w, t),
    sbn( u ),
    nni( n ),
    spr( s ),
    num_accepted_current_period( 0 ),
    num_accepted_total( 0 ),
    relativeMoveWeights( m ),
    waitBeforeRecording( n0 ),
    waitBeforeFitting( c0 ),
    per_move_num_accepted_current_period( std::vector<unsigned int>(3,0) ),
    per_move_num_tried_current_period( std::vector<unsigned int>(3,0) ),
    trace(  )

{
std::cout << "initializing meta move" << std::endl;
    // Add all nodes to move
    sbn->setMove( this );
    nni->setMove( this );
    spr->setMove( this );

    const std::vector<DagNode*>& nodes1 = sbn->getNodes();
    for (size_t j=0; j<nodes1.size(); ++j)
    {
      addNode(nodes1[j]);
    }

    const std::vector<DagNode*>& nodes2 = sbn->getNodes();
    for (size_t j=0; j<nodes2.size(); ++j)
    {
      addNode(nodes2[j]);
    }

    const std::vector<DagNode*>& nodes3 = sbn->getNodes();
    for (size_t j=0; j<nodes3.size(); ++j)
    {
      addNode(nodes3[j]);
    }

    // remove all "core" nodes from affectedNodes so their probabilities are not double-counted
    for (size_t i = 0; i < affected_nodes.size(); ++i)
    {
        RbOrderedSet<DagNode*>::iterator it = affected_nodes.begin();
        std::advance(it, i);

        for (size_t j = 0; j < nodes.size(); ++j)
        {
            if ( nodes[j] == *it )
            {
                affected_nodes.erase(*it);
                --i;
                break;
            }

        }

    }
std::cout << "initialized meta move" << std::endl;

}


/**
 * Copy constructor.
 * We need to create a deep copy of the proposal here.
 *
 * \param[in]   m   The object to copy.
 *
 */
UnconstrainedSBNMetaMove::UnconstrainedSBNMetaMove(const UnconstrainedSBNMetaMove &m) : AbstractMove(m),
    sbn( m.sbn ),
    nni( m.nni ),
    spr( m.spr ),
    num_accepted_current_period( m.num_accepted_current_period ),
    num_accepted_total( m.num_accepted_total ),
    relativeMoveWeights( m.relativeMoveWeights ),
    waitBeforeRecording( m.waitBeforeRecording ),
    waitBeforeFitting( m.waitBeforeRecording ),
    per_move_num_accepted_current_period( m.per_move_num_accepted_current_period ),
    per_move_num_tried_current_period( m.per_move_num_tried_current_period ),
    trace( m.trace )

{
    // Add all nodes to move
    sbn->setMove( this );
    nni->setMove( this );
    spr->setMove( this );

}


/**
 * Basic destructor doing nothing.
 */
UnconstrainedSBNMetaMove::~UnconstrainedSBNMetaMove( void )
{

    delete sbn;
    delete nni;
    delete spr;
}


/**
 * Overloaded assignment operator.
 * We need a deep copy of the operator.
 */
UnconstrainedSBNMetaMove& UnconstrainedSBNMetaMove::operator=(const RevBayesCore::UnconstrainedSBNMetaMove &m)
{
    if ( this != &m )
    {
        // delegate
        AbstractMove::operator=( m );

        // free memory
        delete sbn;
        delete nni;
        delete spr;

        sbn                                  = m.sbn;
        nni                                  = m.nni;
        spr                                  = m.spr;
        num_accepted_current_period          = m.num_accepted_current_period;
        num_accepted_total                   = m.num_accepted_total;
        relativeMoveWeights                  = m.relativeMoveWeights;
        waitBeforeRecording                  = m.waitBeforeRecording;
        waitBeforeFitting           = m.waitBeforeFitting;
        per_move_num_accepted_current_period = m.per_move_num_accepted_current_period;
        per_move_num_tried_current_period    = m.per_move_num_tried_current_period;
        trace                                = m.trace;

        // Add all nodes to move
        sbn->setMove( this );
        nni->setMove( this );
        spr->setMove( this );

    }

    return *this;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the UnconstrainedSBNMetaMove.
 */
UnconstrainedSBNMetaMove* UnconstrainedSBNMetaMove::clone( void ) const
{

    return new UnconstrainedSBNMetaMove( *this );
}

/**
 * This function finds which of the component we should use for this proposal.
 * It uses the relative probabilities of the moves, and masks adaptive moves that we are not yet using.
 *
 * \return The index of the proposal to use for this move
 */
size_t UnconstrainedSBNMetaMove::getActiveProposalIndex( void ) const
{
  double u = GLOBAL_RNG->uniform01();

  double total_weight;

  // We only use the SBN move after we've been able to fit it to tree samples
  total_weight += num_tried_total > waitBeforeFitting ? relativeMoveWeights[0] : 0.0;
  total_weight += relativeMoveWeights[1];
  total_weight += relativeMoveWeights[2];

  u *= total_weight;

  size_t idx;

  if (num_tried_total > waitBeforeFitting)
  {
    if ( u < relativeMoveWeights[0] )
    {
      idx = 0;
    }
    else if ( u < relativeMoveWeights[0] + relativeMoveWeights[1] )
    {
      idx = 1;
    }
    else
    {
      idx = 2;
    }
  }
  else
  {
    if ( u < relativeMoveWeights[1] )
    {
      idx = 1;
    }
    else
    {
      idx = 2;
    }
  }

  return idx;
}

/**
 * Get moves' name of object
 *
 * \return The moves' name.
 */
const std::string& UnconstrainedSBNMetaMove::getMoveName( void ) const
{

    static std::string name = "MixtureMove(SBN+NNI+SPR)";

    return(name);
}


double UnconstrainedSBNMetaMove::getMoveTuningParameter( void ) const
{
    // There are many moves, no single tuning parameter exists
    return 0.0;
}


/**
 * How often was the move accepted
 */
size_t UnconstrainedSBNMetaMove::getNumberAcceptedCurrentPeriod( void ) const
{

    return num_accepted_current_period;
}


/**
 * How often was the move accepted
 */
size_t UnconstrainedSBNMetaMove::getNumberAcceptedTotal( void ) const
{

    return num_accepted_total;
}


// /**
//  * Get the proposal of the move
//  *
//  * \return The proposal object.
//  */
// Proposal& UnconstrainedSBNMetaMove::getProposal( void )
// {
//
//     return *proposal;
// }


void UnconstrainedSBNMetaMove::performHillClimbingMove( double lHeat, double pHeat )
{

    size_t idx = getActiveProposalIndex();
    Proposal* proposal;
    if (idx == 0)
    {
      proposal = sbn;
    }
    else if (idx == 1)
    {
      proposal = nni;
    }
    else
    {
      proposal = spr;
    }

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




void UnconstrainedSBNMetaMove::performMcmcMove( double prHeat, double lHeat, double pHeat )
{
    size_t idx = getActiveProposalIndex();

    ++per_move_num_tried_current_period[idx];

    Proposal* proposal;
    if (idx == 0)
    {
      proposal = sbn;
    }
    else if (idx == 1)
    {
      proposal = nni;
    }
    else
    {
      proposal = spr;
    }

    recordState();

    const RbOrderedSet<DagNode*> &affected_nodes = getAffectedNodes();
    const std::vector<DagNode*> nodes = getDagNodes();

    // Propose a new value
    proposal->prepareProposal();
    double ln_hastings_ratio = proposal->doProposal();


    // Identify nodes that proposal touches
    std::vector<DagNode*> touched_nodes = nodes; //proposal->identifyNodesToTouch();

    // first we touch all the nodes
    // that will set the flags for recomputation
    for (size_t i = 0; i < touched_nodes.size(); ++i)
    {

        // get the pointer to the current node
        DagNode* the_node = touched_nodes[i];

//        // should this node be touched?
//        std::vector<DagNode*>::iterator it = std::find( nodes_left_untouched.begin(), nodes_left_untouched.end(), the_node );
//        bool touch_node = (it == nodes_left_untouched.end());

        // flag for recomputation
//        if (touch_node) {
        the_node->touch();
//        }
    }

    double ln_prior_ratio = 0.0;
    double ln_likelihood_ratio = 0.0;


    // compute the probability of the current value for each node
    for (size_t i = 0; i < touched_nodes.size(); ++i)
    {
        // get the pointer to the current node
        DagNode* the_node = touched_nodes[i];

        if ( RbMath::isAComputableNumber(ln_prior_ratio) && RbMath::isAComputableNumber(ln_likelihood_ratio) && RbMath::isAComputableNumber(ln_hastings_ratio) )
        {
            if ( the_node->isClamped() )
            {
                ln_likelihood_ratio += the_node->getLnProbabilityRatio();
            }
            else
            {
                ln_prior_ratio += the_node->getLnProbabilityRatio();
            }

        }

    }

    // then we recompute the probability for all the affected nodes
    for (RbOrderedSet<DagNode*>::const_iterator it = affected_nodes.begin(); it != affected_nodes.end(); ++it)
    {
        DagNode *the_node = *it;

        if ( RbMath::isAComputableNumber(ln_prior_ratio) && RbMath::isAComputableNumber(ln_likelihood_ratio) && RbMath::isAComputableNumber(ln_hastings_ratio) )
        {
            if ( the_node->isClamped() )
            {
                ln_likelihood_ratio += the_node->getLnProbabilityRatio();
            }
            else
            {
                ln_prior_ratio += the_node->getLnProbabilityRatio();
            }
        }

    }

    // exponentiate with the chain heat
    double ln_posterior_ratio;
    ln_posterior_ratio = pHeat * (lHeat * ln_likelihood_ratio + prHeat * ln_prior_ratio);

    bool rejected = false;

	if ( RbMath::isAComputableNumber(ln_posterior_ratio) == false )
    {
        rejected = true;

        proposal->undoProposal();

        // call restore for each node
        for (size_t i = 0; i < touched_nodes.size(); ++i)
        {
            // get the pointer to the current node
            DagNode* the_node = touched_nodes[i];
            the_node->restore();
        }
	}
    else
    {

        // finally add the Hastings ratio
        double ln_acceptance_ratio = ln_posterior_ratio + ln_hastings_ratio;

        if (ln_acceptance_ratio >= 0.0)
        {


            num_accepted_total++;
            num_accepted_current_period++;

            // call accept for each node
            for (size_t i = 0; i < touched_nodes.size(); ++i)
            {
                // get the pointer to the current node
                DagNode* the_node = touched_nodes[i];
                the_node->keep();
            }

            proposal->cleanProposal();
        }
        else if (ln_acceptance_ratio < -300.0)
        {
            rejected = true;

            proposal->undoProposal();

            // call restore for each node
            for (size_t i = 0; i < touched_nodes.size(); ++i)
            {
                // get the pointer to the current node
                DagNode* the_node = touched_nodes[i];
                the_node->restore();
            }
        }
        else
        {
            double r = exp(ln_acceptance_ratio);
            // Accept or reject the move
            double u = GLOBAL_RNG->uniform01();
            if (u < r)
            {

                num_accepted_total++;
                num_accepted_current_period++;

                // call accept for each node
                for (size_t i = 0; i < touched_nodes.size(); ++i)
                {
                    // get the pointer to the current node
                    DagNode* the_node = touched_nodes[i];
                    the_node->keep();
                }

                proposal->cleanProposal();
            }
            else
            {
                rejected = true;

                proposal->undoProposal();

                // call restore for each node
                for (size_t i = 0; i < touched_nodes.size(); ++i)
                {
                    // get the pointer to the current node
                    DagNode* the_node = touched_nodes[i];
                    the_node->restore();
                }

            }

        }

    }

  if ( rejected == false)
  {
    ++per_move_num_accepted_current_period[idx];
  }
}


/**
 * Print the summary of the move.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void UnconstrainedSBNMetaMove::printSummary(std::ostream &o, bool current_period) const
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

    // proposal->printParameterSummary( o, false );

    o << std::endl;

    o.setf(previousFlags);
    o.precision(previousPrecision);


}

// Record state of chain (after checking if we should be recording yet)
void UnconstrainedSBNMetaMove::recordState( void )
{

}

/**
 * Reset the move counters. Here we only reset the counter for the number of accepted moves.
 *
 */
void UnconstrainedSBNMetaMove::resetMoveCounters( void )
{
    // Main counter for overall move
    num_accepted_current_period = 0;

    // Since our parent class does not know about our component moves, we reset both their number tried and number accepted
    per_move_num_accepted_current_period = std::vector<unsigned int>(3,0);
    per_move_num_tried_current_period = std::vector<unsigned int>(3,0);

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new variable.
 */
void UnconstrainedSBNMetaMove::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    sbn->swapNode(oldN, newN);
    nni->swapNode(oldN, newN);
    spr->swapNode(oldN, newN);
}


void UnconstrainedSBNMetaMove::setMoveTuningParameter(double tp)
{
    // There is no single tuning parameter, nothing to do here
}


void UnconstrainedSBNMetaMove::setNumberAcceptedCurrentPeriod( size_t na )
{
    num_accepted_current_period = na;
}


void UnconstrainedSBNMetaMove::setNumberAcceptedTotal( size_t na )
{
    num_accepted_total = na;
}


/**
 * Tune the move to accept the desired acceptance ratio.
 * We only compute the acceptance ratio here and delegate the call to the proposal.
 */
void UnconstrainedSBNMetaMove::tune( void )
{
    if (num_tried_total > waitBeforeFitting + waitBeforeRecording)
    {
      // Feed the SBN the trace

    }
}
