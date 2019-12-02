#include "DagNode.h"
#include "MetropolisHastingsMetaMove.h"
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
MetropolisHastingsMetaMove::MetropolisHastingsMetaMove( std::vector<Proposal*> *p, double w, bool t ) : AbstractMove(w, t),
    num_accepted_current_period( 0 ),
    num_accepted_total( 0 ),
    proposals( *p ),
    per_move_num_accepted_current_period( std::vector<unsigned int>(proposals.size(),0) ),
    per_move_num_tried_current_period( std::vector<unsigned int>(proposals.size(),0) )

{
    // Add all nodes to move
    for (size_t i=0; i<proposals.size(); ++i)
    {
      proposals.at(i)->setMove( this );
      const std::vector<DagNode*>& nodes = proposals.at(i)->getNodes();
      for (size_t j=0; j<nodes.size(); ++j)
      {
        addNode(nodes[j]);
      }
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

}


/**
 * Copy constructor.
 * We need to create a deep copy of the proposal here.
 *
 * \param[in]   m   The object to copy.
 *
 */
MetropolisHastingsMetaMove::MetropolisHastingsMetaMove(const MetropolisHastingsMetaMove &m) : AbstractMove(m),
    num_accepted_current_period( m.num_accepted_current_period ),
    num_accepted_total( m.num_accepted_total ),
    proposals( m.proposals ),
    per_move_num_accepted_current_period( m.per_move_num_accepted_current_period ),
    per_move_num_tried_current_period( m.per_move_num_tried_current_period )

{

    // Add all nodes to move
    for (size_t i=0; i<proposals.size(); ++i)
    {
      proposals.at(i)->setMove( this );
    }

}


/**
 * Basic destructor doing nothing.
 */
MetropolisHastingsMetaMove::~MetropolisHastingsMetaMove( void )
{

    proposals.clear();
    // delete &proposals;
}


/**
 * Overloaded assignment operator.
 * We need a deep copy of the operator.
 */
MetropolisHastingsMetaMove& MetropolisHastingsMetaMove::operator=(const RevBayesCore::MetropolisHastingsMetaMove &m)
{

    if ( this != &m )
    {
        // delegate
        AbstractMove::operator=( m );

        // free memory
        // delete &proposals;
        proposals.clear();

        num_accepted_current_period          = m.num_accepted_current_period;
        num_accepted_total                   = m.num_accepted_total;
        proposals                            = m.proposals;
        per_move_num_accepted_current_period = m.per_move_num_accepted_current_period;
        per_move_num_tried_current_period    = m.per_move_num_tried_current_period;

        // Add all nodes to move
        for (size_t i=0; i<proposals.size(); ++i)
        {
          proposals.at(i)->setMove( this );
        }

    }

    return *this;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the MetropolisHastingsMetaMove.
 */
MetropolisHastingsMetaMove* MetropolisHastingsMetaMove::clone( void ) const
{

    return new MetropolisHastingsMetaMove( *this );
}

size_t MetropolisHastingsMetaMove::getActiveProposalIndex( void ) const
{
  double u = GLOBAL_RNG->uniform01();
  size_t idx = std::floor(u * proposals.size());
  return idx;
}

/**
 * Get moves' name of object
 *
 * \return The moves' name.
 */
const std::string& MetropolisHastingsMetaMove::getMoveName( void ) const
{

    std::string s;

    s += "MixtureMove(";

    s += proposals.at(0)->getProposalName();

    if (proposals.size() > 1)
    {
      for (size_t i=1; i<proposals.size(); ++i)
      {
        s += "+";
        s += proposals.at(i)->getProposalName();
      }
    }
    s += ")";

    // static std::string name = "KernelMixture";

    static std::string name = s;

    return(name);
}


double MetropolisHastingsMetaMove::getMoveTuningParameter( void ) const
{
    // There are many moves, no single tuning parameter exists
    return 0.0;
}


/**
 * How often was the move accepted
 */
size_t MetropolisHastingsMetaMove::getNumberAcceptedCurrentPeriod( void ) const
{

    return num_accepted_current_period;
}


/**
 * How often was the move accepted
 */
size_t MetropolisHastingsMetaMove::getNumberAcceptedTotal( void ) const
{

    return num_accepted_total;
}


// /**
//  * Get the proposal of the move
//  *
//  * \return The proposal object.
//  */
// Proposal& MetropolisHastingsMetaMove::getProposal( void )
// {
//
//     return *proposal;
// }


void MetropolisHastingsMetaMove::performHillClimbingMove( double lHeat, double pHeat )
{

    size_t idx = getActiveProposalIndex();
    Proposal* proposal = proposals.at(idx);

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




void MetropolisHastingsMetaMove::performMcmcMove( double prHeat, double lHeat, double pHeat )
{
    size_t idx = getActiveProposalIndex();

    ++per_move_num_tried_current_period[idx];

    Proposal* proposal = proposals.at(idx);

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
void MetropolisHastingsMetaMove::printSummary(std::ostream &o, bool current_period) const
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

/**
 * Reset the move counters. Here we only reset the counter for the number of accepted moves.
 *
 */
void MetropolisHastingsMetaMove::resetMoveCounters( void )
{
    // Main counter for overall move
    num_accepted_current_period = 0;

    // Since our parent class does not know about our component moves, we reset both their number tried and number accepted
    per_move_num_accepted_current_period = std::vector<unsigned int>(proposals.size(),0);
    per_move_num_tried_current_period = std::vector<unsigned int>(proposals.size(),0);

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new variable.
 */
void MetropolisHastingsMetaMove::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    for (size_t i=0; i<proposals.size(); ++i)
    {
      proposals.at(i)->swapNode(oldN, newN);
    }

}


void MetropolisHastingsMetaMove::setMoveTuningParameter(double tp)
{
    // There is no single tuning parameter, nothing to do here
}


void MetropolisHastingsMetaMove::setNumberAcceptedCurrentPeriod( size_t na )
{
    num_accepted_current_period = na;
}


void MetropolisHastingsMetaMove::setNumberAcceptedTotal( size_t na )
{
    num_accepted_total = na;
}


/**
 * Tune the move to accept the desired acceptance ratio.
 * We only compute the acceptance ratio here and delegate the call to the proposal.
 */
void MetropolisHastingsMetaMove::tune( void )
{

    for (size_t i=0; i<proposals.size(); ++i)
    {

      if ( per_move_num_tried_current_period[i] > 2 )
      {
          double rate = double(per_move_num_accepted_current_period[i]) / double(per_move_num_tried_current_period[i]);

          proposals.at(i)->tune( rate );
      }
    }

}
